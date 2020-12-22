/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *<LicenseText>
 *
 * CitcomS by Louis Moresi, Shijie Zhong, Lijie Han, Eh Tan,
 * Clint Conrad, Michael Gurnis, and Eun-seo Choi.
 * Copyright (C) 1994-2005, California Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *</LicenseText>
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/*   Functions which solve the heat transport equations using Petrov-Galerkin
     streamline-upwind methods. The process is basically as described in Alex
     Brooks PhD thesis (Caltech) which refers back to Hughes, Liu and Brooks.  */

#include <sys/types.h>

#include "element_definitions.h"
#include "global_defs.h"
#include <math.h>
#include "advection_diffusion.h"
#include "parsing.h"
#include "initial_temperature.h"

static void set_diffusion_timestep(struct All_variables *E);
static void predictor(struct All_variables *E, double **field,
                      double **fielddot);
static void corrector(struct All_variables *E, double **field,
                      double **fielddot, double **Dfielddot);
static void pg_solver(struct All_variables *E,
                      double **T, double **Tdot, double **DTdot,
                      struct SOURCES *Q0,
                      double diff, int bc, unsigned int **FLAGS);
static void pg_shape_fn(struct All_variables *E, int el,
                        struct Shape_function *PG,
                        struct Shape_function_dx *GNx,
                        float VV[4][9], double rtf[4][9],
                        double diffusion, int m);
static void element_residual(struct All_variables *E, int el,
                             struct Shape_function *PG,
                             struct Shape_function_dx *GNx,
                             struct Shape_function_dA *dOmega,
                             float VV[4][9],
                             double **field, double **fielddot,
                             struct SOURCES *Q0,
                             double Eres[9], double rtf[4][9],
                             double diff, float **BC,
                             unsigned int **FLAGS, int m);
static void filter(struct All_variables *E);
static void process_heating(struct All_variables *E, int psc_pass);

/* ============================================
   Generic adv-diffusion for temperature field.
   ============================================ */


/***************************************************************/

void advection_diffusion_parameters(struct All_variables *E)
{

    /* Set intial values, defaults & read parameters*/
    int m=E->parallel.me;

    input_boolean("ADV",&(E->advection.ADVECTION),"on",m);
    input_boolean("filter_temp",&(E->advection.filter_temperature),"off",m);
    input_boolean("monitor_max_T",&(E->advection.monitor_max_T),"on",m);

    input_int("minstep",&(E->advection.min_timesteps),"1",m);
    input_int("maxstep",&(E->advection.max_timesteps),"1000",m);
    input_int("maxtotstep",&(E->advection.max_total_timesteps),"1000000",m);
    input_float("finetunedt",&(E->advection.fine_tune_dt),"0.9",m);
    input_float("fixed_timestep",&(E->advection.fixed_timestep),"0.0",m);
    input_float("adv_gamma",&(E->advection.gamma),"0.5",m);
    input_int("adv_sub_iterations",&(E->advection.temp_iterations),"2,1,nomax",m);

    input_float("inputdiffusivity",&(E->control.inputdiff),"1.0",m);


    return;
}


void advection_diffusion_allocate_memory(struct All_variables *E)
{
  int i,m;

  for(m=1;m<=E->sphere.caps_per_proc;m++)  {
    E->Tdot[m]= (double *)malloc((E->lmesh.nno+1)*sizeof(double));

    for(i=1;i<=E->lmesh.nno;i++)
      E->Tdot[m][i]=0.0;
    }

  return;
}


void PG_timestep_init(struct All_variables *E)
{

  set_diffusion_timestep(E);

  return;
}


void PG_timestep(struct All_variables *E)
{
    void std_timestep();
    void PG_timestep_solve();

    std_timestep(E);

    PG_timestep_solve(E);

    return;
}



/* =====================================================
   Obtain largest possible timestep (no melt considered)
   =====================================================  */


void std_timestep(struct All_variables *E)
{
    int i,d,n,nel,el,node,m;

    float global_fmin();
    void velo_from_element();

    float adv_timestep;
    float ts,uc1,uc2,uc3,uc,size,step,VV[4][9];

    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int nno=E->lmesh.nno;
    const int lev=E->mesh.levmax;
    const int ends=enodes[dims];
    const int sphere_key = 1;

    nel=E->lmesh.nel;

    if(E->advection.fixed_timestep != 0.0) {
      E->advection.timestep = E->advection.fixed_timestep;
      return;
    }

    adv_timestep = 1.0e8;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(el=1;el<=nel;el++) {

	velo_from_element(E,VV,m,el,sphere_key);

	uc=uc1=uc2=uc3=0.0;
	for(i=1;i<=ENODES3D;i++) {
	  uc1 += E->N.ppt[GNPINDEX(i,1)]*VV[1][i];
	  uc2 += E->N.ppt[GNPINDEX(i,1)]*VV[2][i];
	  uc3 += E->N.ppt[GNPINDEX(i,1)]*VV[3][i];
        }
	uc = fabs(uc1)/E->eco[m][el].size[1] + fabs(uc2)/E->eco[m][el].size[2] + fabs(uc3)/E->eco[m][el].size[3];

	step = (0.5/uc);
	adv_timestep = min(adv_timestep,step);
      }

    adv_timestep = E->advection.dt_reduced * adv_timestep;

    adv_timestep = 1.0e-32 + min(E->advection.fine_tune_dt*adv_timestep,
				 E->advection.diff_timestep);

    E->advection.timestep = global_fmin(E,adv_timestep);

/*     if (E->parallel.me==0) */
/*       fprintf(stderr, "adv_timestep=%g diff_timestep=%g\n",adv_timestep,E->advection.diff_timestep); */

    return;
}


void PG_timestep_solve(struct All_variables *E)
{

  double Tmaxd();
  void temperatures_conform_bcs();
  void lith_age_conform_tbc();
  void assimilate_lith_conform_bcs();
  int i,m,psc_pass,iredo;
  double time0,time1,T_interior1;
  double *DTdot[NCS], *T1[NCS], *Tdot1[NCS];

  E->advection.timesteps++;

  for(m=1;m<=E->sphere.caps_per_proc;m++)
    DTdot[m]= (double *)malloc((E->lmesh.nno+1)*sizeof(double));


  if(E->advection.monitor_max_T) {
     for(m=1;m<=E->sphere.caps_per_proc;m++)  {
         T1[m]= (double *)malloc((E->lmesh.nno+1)*sizeof(double));
         Tdot1[m]= (double *)malloc((E->lmesh.nno+1)*sizeof(double));
     }

     for(m=1;m<=E->sphere.caps_per_proc;m++)
         for (i=1;i<=E->lmesh.nno;i++)   {
             T1[m][i] = E->T[m][i];
             Tdot1[m][i] = E->Tdot[m][i];
         }

     /* get the max temperature for old T */
     T_interior1 = Tmaxd(E,E->T);
  }

  E->advection.dt_reduced = 1.0;
  E->advection.last_sub_iterations = 1;


  do {
    E->advection.timestep *= E->advection.dt_reduced;

    iredo = 0;
    if (E->advection.ADVECTION) {

      predictor(E,E->T,E->Tdot);

      for(psc_pass=0;psc_pass<E->advection.temp_iterations;psc_pass++)   {
        /* adiabatic, dissipative and latent heating*/
        if(E->control.disptn_number != 0)
          process_heating(E, psc_pass);

        /* XXX: replace inputdiff with refstate.thermal_conductivity */
	pg_solver(E,E->T,E->Tdot,DTdot,&(E->convection.heat_sources),E->control.inputdiff,1,E->node);
	corrector(E,E->T,E->Tdot,DTdot);
	temperatures_conform_bcs(E);
      }

      if(E->advection.monitor_max_T) {
          /* get the max temperature for new T */
          E->monitor.T_interior = Tmaxd(E,E->T);

          /* if the max temperature changes too much, restore the old
           * temperature field, calling the temperature solver using
           * half of the timestep size */
          if (E->monitor.T_interior/T_interior1 > E->monitor.T_maxvaried) {
              if(E->parallel.me==0) {
                  fprintf(stderr, "max T varied from %e to %e\n",
                          T_interior1, E->monitor.T_interior);
                  fprintf(E->fp, "max T varied from %e to %e\n",
                          T_interior1, E->monitor.T_interior);
              }
              for(m=1;m<=E->sphere.caps_per_proc;m++)
                  for (i=1;i<=E->lmesh.nno;i++)   {
                      E->T[m][i] = T1[m][i];
                      E->Tdot[m][i] = Tdot1[m][i];
                  }
              iredo = 1;
              E->advection.dt_reduced *= 0.5;
              E->advection.last_sub_iterations ++;
          }
      }
    }

  }  while ( iredo==1 && E->advection.last_sub_iterations <= 5);

  	E->advection.dt_before=E->advection.timestep;

  /*assign value to difference in T, by lhy*/
  for (i=1;i<=E->lmesh.nno;i++)   {
  	E->Tdiff[1][i] = E->T[1][i] - T1[1][i];
  }


  /* filter temperature to remove over-/under-shoot */
  if(E->advection.filter_temperature)
    filter(E);


  E->advection.total_timesteps++;
  E->monitor.elapsed_time += E->advection.timestep;

  if (E->advection.last_sub_iterations==5)
    E->control.keep_going = 0;

  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    free((void *) DTdot[m] );
  }

  if(E->advection.monitor_max_T) {
      for(m=1;m<=E->sphere.caps_per_proc;m++) {
          free((void *) T1[m] );
          free((void *) Tdot1[m] );
      }
  }

  if(E->control.lith_age) {
      if(E->parallel.me==0) fprintf(stderr,"PG_timestep_solve\n");
      lith_age_conform_tbc(E);
      assimilate_lith_conform_bcs(E);
  }

  return;
}


/***************************************************************/

static void set_diffusion_timestep(struct All_variables *E)
{
  float diff_timestep, ts;
  int m, el, d;

  float global_fmin();

  diff_timestep = 1.0e8;
  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(el=1;el<=E->lmesh.nel;el++)  {
      for(d=1;d<=E->mesh.nsd;d++)    {
	ts = E->eco[m][el].size[d] * E->eco[m][el].size[d];
	diff_timestep = min(diff_timestep,ts);
      }
    }

  diff_timestep = global_fmin(E,diff_timestep);
  E->advection.diff_timestep = 0.5 * diff_timestep;

  return;
}


/* ==============================
   predictor and corrector steps.
   ============================== */

static void predictor(struct All_variables *E, double **field,
                      double **fielddot)
{
  int node,m;
  double multiplier;

  multiplier = (1.0-E->advection.gamma) * E->advection.timestep;

  for (m=1;m<=E->sphere.caps_per_proc;m++)
    for(node=1;node<=E->lmesh.nno;node++)  {
      field[m][node] += multiplier * fielddot[m][node] ;
      fielddot[m][node] = 0.0;
    }

  return;
}


static void corrector(struct All_variables *E, double **field,
                      double **fielddot, double **Dfielddot)
{
  int node,m;
  double multiplier;

  multiplier = E->advection.gamma * E->advection.timestep;

  for (m=1;m<=E->sphere.caps_per_proc;m++)
    for(node=1;node<=E->lmesh.nno;node++) {
      field[m][node] += multiplier * Dfielddot[m][node];
      fielddot[m][node] +=  Dfielddot[m][node];
    }

  return;
}


/* ===================================================
   The solution step -- determine residual vector from
   advective-diffusive terms and solve for delta Tdot
   Two versions are available -- one for Cray-style
   vector optimizations etc and one optimized for
   workstations.
   =================================================== */


static void pg_solver(struct All_variables *E,
                      double **T, double **Tdot, double **DTdot,
                      struct SOURCES *Q0,
                      double diff, int bc, unsigned int **FLAGS)
{
    void get_rtf_at_vpts();
    void velo_from_element();

    int el,e,a,i,a1,m;
    double Eres[9],rtf[4][9];  /* correction to the (scalar) Tdot field */
    float VV[4][9];

    struct Shape_function PG;

    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int ends=enodes[dims];
    const int sphere_key = 1;
    const int lev=E->mesh.levmax;

    for (m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=1;i<=E->lmesh.nno;i++)
 	 DTdot[m][i] = 0.0;

    for (m=1;m<=E->sphere.caps_per_proc;m++)
       for(el=1;el<=E->lmesh.nel;el++)    {

          velo_from_element(E,VV,m,el,sphere_key);

          get_rtf_at_vpts(E, m, lev, el, rtf);

          /* XXX: replace diff with refstate.thermal_conductivity */
          pg_shape_fn(E, el, &PG, &(E->gNX[m][el]), VV,
                      rtf, diff, m);
          element_residual(E, el, &PG, &(E->gNX[m][el]), &(E->gDA[m][el]),
                           VV, T, Tdot,
                           Q0, Eres, rtf, diff, E->sphere.cap[m].TB,
                           FLAGS, m);

        for(a=1;a<=ends;a++) {
	    a1 = E->ien[m][el].node[a];
	    DTdot[m][a1] += Eres[a];
           }

        } /* next element */

    (E->exchange_node_d)(E,DTdot,lev);

    for (m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=1;i<=E->lmesh.nno;i++) {
        if(!(E->node[m][i] & (TBX | TBY | TBZ))){
	  DTdot[m][i] *= E->TMass[m][i];         /* lumped mass matrix */
	}	else
	  DTdot[m][i] = 0.0;         /* lumped mass matrix */
      }
    return;
}



/* ===================================================
   Petrov-Galerkin shape functions for a given element
   =================================================== */

static void pg_shape_fn(struct All_variables *E, int el,
                        struct Shape_function *PG,
                        struct Shape_function_dx *GNx,
                        float VV[4][9], double rtf[4][9],
                        double diffusion, int m)
{
    int i,j;
    int *ienm;

    double uc1,uc2,uc3;
    double u1,u2,u3,sint[9];
    double uxse,ueta,ufai,xse,eta,fai,adiff;

    double prod1,unorm,twodiff;

    ienm=E->ien[m][el].node;

    twodiff = 2.0*diffusion;

    uc1 =  uc2 = uc3 = 0.0;

    for(i=1;i<=ENODES3D;i++) {
      uc1 +=  E->N.ppt[GNPINDEX(i,1)]*VV[1][i];
      uc2 +=  E->N.ppt[GNPINDEX(i,1)]*VV[2][i];
      uc3 +=  E->N.ppt[GNPINDEX(i,1)]*VV[3][i];
      }

    uxse = fabs(uc1*E->eco[m][el].size[1]);
    ueta = fabs(uc2*E->eco[m][el].size[2]);
    ufai = fabs(uc3*E->eco[m][el].size[3]);

    xse = (uxse>twodiff)? (1.0-twodiff/uxse):0.0;
    eta = (ueta>twodiff)? (1.0-twodiff/ueta):0.0;
    fai = (ufai>twodiff)? (1.0-twodiff/ufai):0.0;


    unorm = uc1*uc1 + uc2*uc2 + uc3*uc3;

    adiff = (unorm>0.000001)?( (uxse*xse+ueta*eta+ufai*fai)/(2.0*unorm) ):0.0;

    for(i=1;i<=VPOINTS3D;i++)
       sint[i] = rtf[3][i]/sin(rtf[1][i]);

    for(i=1;i<=VPOINTS3D;i++) {
       u1 = u2 = u3 = 0.0;
       for(j=1;j<=ENODES3D;j++)  /* this line heavily used */ {
		u1 += VV[1][j] * E->N.vpt[GNVINDEX(j,i)];
		u2 += VV[2][j] * E->N.vpt[GNVINDEX(j,i)];
	   	u3 += VV[3][j] * E->N.vpt[GNVINDEX(j,i)];
	    }

       for(j=1;j<=ENODES3D;j++) {
            prod1 = (u1 * GNx->vpt[GNVXINDEX(0,j,i)]*rtf[3][i] +
                     u2 * GNx->vpt[GNVXINDEX(1,j,i)]*sint[i] +
                     u3 * GNx->vpt[GNVXINDEX(2,j,i)] ) ;

	    PG->vpt[GNVINDEX(j,i)] = E->N.vpt[GNVINDEX(j,i)] + adiff * prod1;
	    }
       }

   return;
}



/* ==========================================
   Residual force vector from heat-transport.
   Used to correct the Tdot term.
   =========================================  */

static void element_residual(struct All_variables *E, int el,
                             struct Shape_function *PG,
                             struct Shape_function_dx *GNx,
                             struct Shape_function_dA *dOmega,
                             float VV[4][9],
                             double **field, double **fielddot,
                             struct SOURCES *Q0,
                             double Eres[9], double rtf[4][9],
                             double diff, float **BC,
                             unsigned int **FLAGS, int m)
{
	int Is00 = (E->parallel.me_loc[1]==0)&&(E->parallel.me_loc[2]==0);
	int chatty = (E->parallel.me==0)&&(el==1);
	
	if(chatty){
		fprintf(stderr,"element_residual starts\n");
	}//debug
    int i,j,a,k,node,nodes[5],d,aid,back_front,onedfns;
    double Q, QER, QER1;
    double dT[9];
    double tx1[9],tx2[9],tx3[9],sint[9];
    double v1[9],v2[9],v3[9];
    double adv_dT,t2[4];
    double T,DT;

    double prod,sfn;
    struct Shape_function1 GM;
    struct Shape_function1_dA dGamma;
    double temp,temp1,rho,cp,heating,Q_enrich;
	double fac1, fac2;
    int nz;
	const int en_flavor=E->trace.enriched_flavor;
	const int en_flavor1=E->trace.enriched_flavor1;
	if (E->parallel.me==0&&el==1)
		fprintf(stderr,"en_flavor1 = %d\n",en_flavor1); //lhy debug

    void get_global_1d_shape_fn();

    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int nno=E->lmesh.nno;
    const int lev=E->mesh.levmax;
    const int ends=enodes[dims];
    const int vpts=vpoints[dims];
    const int onedvpts = onedvpoints[dims];
    const int diffusion = (diff != 0.0);
	double **melt_fac; //lhy for calculating latent heat
	const int Is0 = (E->parallel.me==0);

	melt_fac = (double**)malloc((vpts+1)*sizeof(double));
	for(i=0;i<=vpts;i++){
		melt_fac[i] = (double*)malloc(3*sizeof(double));
		melt_fac[i][1] = 1.0; //assign initial value
		melt_fac[i][2] = 0.0;
	}

    for(i=1;i<=vpts;i++)	{
      dT[i]=0.0;
      v1[i] = tx1[i]=  0.0;
      v2[i] = tx2[i]=  0.0;
      v3[i] = tx3[i]=  0.0;
      }

    for(i=1;i<=vpts;i++)
        sint[i] = rtf[3][i]/sin(rtf[1][i]);

    for(j=1;j<=ends;j++)       {
      node = E->ien[m][el].node[j];
      T = field[m][node];
      if(E->node[m][node] & (TBX | TBY | TBZ))
	    DT=0.0;
      else
	    DT = fielddot[m][node];

      for(i=1;i<=vpts;i++)  {
          dT[i] += DT * E->N.vpt[GNVINDEX(j,i)];
          tx1[i] += GNx->vpt[GNVXINDEX(0,j,i)] * T * rtf[3][i];
          tx2[i] += GNx->vpt[GNVXINDEX(1,j,i)] * T * sint[i];
          tx3[i] += GNx->vpt[GNVXINDEX(2,j,i)] * T;
          sfn = E->N.vpt[GNVINDEX(j,i)];
		  //v1,v2,v3 is r-th-fai velocity on gauss point
          v1[i] += VV[1][j] * sfn;
          v2[i] += VV[2][j] * sfn;
          v3[i] += VV[3][j] * sfn;
      }
    }

	if (Q0->number>0){ //lhy 20180305
    	Q=0.0;
    	for(i=0;i<Q0->number;i++){
	  	temp = Q0->Q[i] * exp(-Q0->lambda[i] * (E->monitor.elapsed_time+Q0->t_offset));
	  	Q += Q0->Q[i] * exp(-Q0->lambda[i] * (E->monitor.elapsed_time+Q0->t_offset));
		}
	}
	else{
    /* heat production */
    Q = E->control.Q0;
	}
	if(E->parallel.me==0&&el==1){
		fprintf(stderr,"elapse_time: %.4e, Q: %.4e\n",E->monitor.elapsed_time,Q); //lhy debug
	}

    /* should we add a compositional contribution? */
    if(E->control.tracer_enriched){
		if(Q0->number > 0){
		QER = 0.0;
    	for(i=0;i<Q0->number;i++)
	  		QER += Q0->QER[i] * exp(-Q0->lambda[i] * (E->monitor.elapsed_time+Q0->t_offset));
		}
		else
		QER = E->control.Q0ER;

		if(en_flavor1>0){
			if(Q0->number > 0){
			QER1 = 0.0;
    		for(i=0;i<Q0->number;i++)
	  			QER1 += Q0->QER1[i] * exp(-Q0->lambda[i] * (E->monitor.elapsed_time+Q0->t_offset));
			}
			temp1 = E->composition.comp_el[m][en_flavor1-1][el];
		}
		else{
			temp1 = 0.0;
			QER1 = 0.0;
		}

      /* XXX: change Q and Q0 to be a vector of ncomp elements */

      /* Q = Q0 for C = 0, Q = Q0ER for C = 1, and linearly in
	 between  */
	 	if(en_flavor>0){
			temp = E->composition.comp_el[m][en_flavor-1][el];
		}
      	Q *= (1.0 - temp - temp1);
      	Q += temp * QER + temp1 * QER1;
	 }
	 
	if(E->parallel.me==1&&el<=E->lmesh.elz)
		fprintf(stderr,"el = %d, Q = %.4e\n", el, Q); //lhy debug
/*	
	if(E->parallel.me==1&&el<=E->lmesh.noz)
		fprintf(stderr,"proc1, el: %d, Q: %.4e\n",el,Q); //lhy debug
*/

    nz = ((el-1) % E->lmesh.elz) + 1;
    rho = 0.5 * (E->refstate.rho[nz] + E->refstate.rho[nz+1]);
    cp = 0.5 * (E->refstate.heat_capacity[nz] + E->refstate.heat_capacity[nz+1]);
	/*
	if(Is00&&el<=E->lmesh.elz){ //lhy debug
		fprintf(stderr,"procz=%d, nz=%d, rho=%.4e\n",E->parallel.me_loc[3],nz,rho);
	}
	*/

	if(E->control.melt_latent_heating!=0)
	{
		if(Is0&&el==10){ //lhy debug
			fprintf(stderr,"substract latent heat from Q, latent_heat is %.4e\n",E->melt_latent_heat[el]);
		}
		Q=Q-E->melt_latent_heat[el];
	}

    if(E->control.disptn_number == 0)
        heating = rho * Q;
    else
        /* E->heating_latent is actually the inverse of latent heating */
        heating = (rho * Q - E->heating_adi[m][el] + E->heating_visc[m][el])
            * E->heating_latent[m][el];
	/*derive melting additional factor*/
	if(E->lunar.latent_method!=0){
		process_melt_factor(E, melt_fac,el,E->lunar.latent_method);
	}
		if(Is0&&el==10){
			fprintf(stderr,"el = %d, melt_fac[1][1] = %.4e, melt_fac[1][2] = %.4e\n",el,melt_fac[1][1],melt_fac[1][2]);
		} //debug

    /* construct residual from this information */

	
    if(diffusion){
      for(j=1;j<=ends;j++) {
	Eres[j]=0.0;
	for(i=1;i<=vpts;i++){
	  if(E->lunar.latent_method){
		  fac1 = melt_fac[i][1];
		  fac2 = melt_fac[i][2];
	  }
	  else{
		  fac1 = 1.0;
		  fac2 = 0.0;
	  }
	  Eres[j] -=
	    PG->vpt[GNVINDEX(j,i)] * dOmega->vpt[i]
              * ((dT[i] + v1[i]*tx1[i] + v2[i]*tx2[i] + v3[i]*tx3[i])*rho*cp
                 - (heating +  rho*fac2)/fac1 )
              + diff * dOmega->vpt[i] * E->heating_latent[m][el]
              * (GNx->vpt[GNVXINDEX(0,j,i)]*tx1[i]*rtf[3][i] +
                 GNx->vpt[GNVXINDEX(1,j,i)]*tx2[i]*sint[i] +
                 GNx->vpt[GNVXINDEX(2,j,i)]*tx3[i] )/fac1; //ongoing test
	 /* if(Is0&&el==10){
          temp = diff * dOmega->vpt[i] * E->heating_latent[m][el]
              * (GNx->vpt[GNVXINDEX(0,j,i)]*tx1[i]*rtf[3][i] +
                 GNx->vpt[GNVXINDEX(1,j,i)]*tx2[i]*sint[i] +
                 GNx->vpt[GNVXINDEX(2,j,i)]*tx3[i] );
		  fprintf(stderr,"i:%d, j:%d, temp:%.4e\n",i,j,temp);
		  //debug
	  }*/
	  }
      }
    }

    else { /* no diffusion term */
      for(j=1;j<=ends;j++) {
	Eres[j]=0.0;
	for(i=1;i<=vpts;i++)
	  Eres[j] -= PG->vpt[GNVINDEX(j,i)] * dOmega->vpt[i]
              * (dT[i] - heating + v1[i]*tx1[i] + v2[i]*tx2[i] + v3[i]*tx3[i]);
      }
    }

    /* See brooks etc: the diffusive term is excused upwinding for
       rectangular elements  */

    /* include BC's for fluxes at (nominally horizontal) edges (X-Y plane) */

    if(FLAGS!=NULL) {
      aid = -1;
      if (FLAGS[m][E->ien[m][el].node[1]] & FBZ) {   // only check for the 1st node
          aid = 0;
	  get_global_1d_shape_fn(E,el,&GM,&dGamma,aid,m);
          }
      else if (FLAGS[m][E->ien[m][el].node[5]] & FBZ) {   // only check for the 5th node
          aid = 1;
	  get_global_1d_shape_fn(E,el,&GM,&dGamma,aid,m);
          }
      if (aid>=0)  {
        for(a=1;a<=onedvpts;a++)  {

	  for(j=1;j<=onedvpts;j++)  {
            dT[j] = 0.0;
	    for(k=1;k<=onedvpts;k++)
              dT[j] += E->M.vpt[GMVINDEX(k,j)]*BC[3][E->ien[m][el].node[k+aid*onedvpts]];
            }
	  for(j=1;j<=onedvpts;j++)  {
	    Eres[a+aid*onedvpts] += dGamma.vpt[GMVGAMMA(aid,j)] *
		E->M.vpt[GMVINDEX(a,j)] * g_1d[j].weight[dims-1] *
		dT[j];
            }

	}
      }
    }
	
	for(i=0;i<=vpts;i++){
		free(melt_fac[i]);
	}
	free(melt_fac);
	
	if(chatty){
		fprintf(stderr,"element_residual ends\n");
	}
    return;
}


/* This function filters the temperature field. The temperature above   */
/* Tmax0(==1.0) and Tmin0(==0.0) is removed, while conserving the total */
/* energy. See Lenardic and Kaula, JGR, 1993.                           */
static void filter(struct All_variables *E)
{
    double Tsum0,Tmin,Tmax,Tsum1,TDIST,TDIST1;
    int m,i;
    double Tmax1,Tmin1;
    double *rhocp, sum_rhocp, total_sum_rhocp;
    int lev, nz;

    /* min and max temperature for filtering */
    const double Tmin0 = 0.0;
    const double Tmax0 = 1.0;

    Tsum0= Tsum1= 0.0;
    Tmin= Tmax= 0.0;
    Tmin1= Tmax1= 0.0;
    TDIST= TDIST1= 0.0;
    sum_rhocp = 0.0;

    lev=E->mesh.levmax;

    rhocp = (double *)malloc((E->lmesh.noz+1)*sizeof(double));
    for(i=1;i<=E->lmesh.noz;i++)
        rhocp[i] = E->refstate.rho[i] * E->refstate.heat_capacity[i];

    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=1;i<=E->lmesh.nno;i++)  {
            nz = ((i-1) % E->lmesh.noz) + 1;

            /* compute sum(rho*cp*T) before filtering, skipping nodes
               that's shared by another processor */
            if(!(E->NODE[lev][m][i] & SKIP))
                Tsum0 += E->T[m][i]*rhocp[nz];

            /* remove overshoot */
            if(E->T[m][i]<Tmin)  Tmin=E->T[m][i];
            if(E->T[m][i]<Tmin0) E->T[m][i]=Tmin0;
            if(E->T[m][i]>Tmax) Tmax=E->T[m][i];
            if(E->T[m][i]>Tmax0) E->T[m][i]=Tmax0;

        }

    /* find global max/min of temperature */
    MPI_Allreduce(&Tmin,&Tmin1,1,MPI_DOUBLE,MPI_MIN,E->parallel.world);
    MPI_Allreduce(&Tmax,&Tmax1,1,MPI_DOUBLE,MPI_MAX,E->parallel.world);

    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=1;i<=E->lmesh.nno;i++)  {
            nz = ((i-1) % E->lmesh.noz) + 1;

            /* remvoe undershoot */
            if(E->T[m][i]<=fabs(2*Tmin0-Tmin1))   E->T[m][i]=Tmin0;
            if(E->T[m][i]>=(2*Tmax0-Tmax1))   E->T[m][i]=Tmax0;

            /* sum(rho*cp*T) after filtering */
            if (!(E->NODE[lev][m][i] & SKIP))  {
                Tsum1 += E->T[m][i]*rhocp[nz];
                if(E->T[m][i]!=Tmin0 && E->T[m][i]!=Tmax0) {
                    sum_rhocp += rhocp[nz];
                }

            }

        }

    /* find the difference of sum(rho*cp*T) before/after the filtering */
    TDIST=Tsum0-Tsum1;
    MPI_Allreduce(&TDIST,&TDIST1,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    MPI_Allreduce(&sum_rhocp,&total_sum_rhocp,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    TDIST=TDIST1/total_sum_rhocp;

    /* keep sum(rho*cp*T) the same before/after the filtering by distributing
       the difference back to nodes */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=1;i<=E->lmesh.nno;i++)   {
            if(E->T[m][i]!=Tmin0 && E->T[m][i]!=Tmax0)
                E->T[m][i] +=TDIST;
        }

    free(rhocp);
    return;
}


static void process_visc_heating(struct All_variables *E, int m,
                                 double *heating)
{
    void strain_rate_2_inv();
    int e, i;
    double visc, temp;
    float *strain_sqr;
    const int vpts = VPOINTS3D;

    strain_sqr = (float*) malloc((E->lmesh.nel+1)*sizeof(float));
    /* note that this will be negative for Atemp < 0 for time
       reversal */
    temp = E->control.disptn_number / E->control.Atemp / vpts;

    strain_rate_2_inv(E, m, strain_sqr, 0);

    for(e=1; e<=E->lmesh.nel; e++) {
        visc = 0.0;
        for(i = 1; i <= vpts; i++)
            visc += E->EVi[m][(e-1)*vpts + i];

        heating[e] = temp * visc * strain_sqr[e];
    }

    free(strain_sqr);

    return;
}


static void process_adi_heating(struct All_variables *E, int m,
                                double *heating)
{
    int e, ez, i, j;
    double matprop, temp1, temp2;
    const int ends = ENODES3D;

    temp2 = E->control.disptn_number / ends;
    for(e=1; e<=E->lmesh.nel; e++) {
        ez = (e - 1) % E->lmesh.elz + 1;
        matprop = 0.125
            * (E->refstate.thermal_expansivity[ez] +
               E->refstate.thermal_expansivity[ez + 1])
            * (E->refstate.rho[ez] + E->refstate.rho[ez + 1])
            * (E->refstate.gravity[ez] + E->refstate.gravity[ez + 1]);

        temp1 = 0.0;
        for(i=1; i<=ends; i++) {
            j = E->ien[m][e].node[i];
            temp1 += E->sphere.cap[m].V[3][j]
                * (E->T[m][j] + E->control.surface_temp);
        }

        heating[e] = matprop * temp1 * temp2;
    }

    return;
}


static void latent_heating(struct All_variables *E, int m,
                           double *heating_latent, double *heating_adi,
                           float **B, float Ra, float clapeyron,
                           float depth, float transT, float inv_width)
{
    double temp, temp0, temp1, temp2, temp3, matprop;
    int e, ez, i, j;
    const int ends = ENODES3D;
    /* 
       note that this will be negative for time-reversal
    */
    temp0 = 2.0 * inv_width * clapeyron * E->control.disptn_number * Ra / E->control.Atemp / ends;
    temp1 = temp0 * clapeyron;

    for(e=1; e<=E->lmesh.nel; e++) {
        ez = (e - 1) % E->lmesh.elz + 1;
        matprop = 0.125
            * (E->refstate.thermal_expansivity[ez] +
               E->refstate.thermal_expansivity[ez + 1])
            * (E->refstate.rho[ez] + E->refstate.rho[ez + 1])
            * (E->refstate.gravity[ez] + E->refstate.gravity[ez + 1]);

        temp2 = 0;
        temp3 = 0;
        for(i=1; i<=ends; i++) {
            j = E->ien[m][e].node[i];
            temp = (1.0 - B[m][j]) * B[m][j]
                * (E->T[m][j] + E->control.surface_temp);
            temp2 += temp * E->sphere.cap[m].V[3][j];
            temp3 += temp;
        }

        /* correction on the adiabatic cooling term */
        heating_adi[e] += matprop * temp2 * temp0;

        /* correction on the DT/Dt term */
        heating_latent[e] += temp3 * temp1;
    }
    return;
}


static void process_latent_heating(struct All_variables *E, int m,
                                   double *heating_latent, double *heating_adi)
{
    int e;

    /* reset */
    for(e=1; e<=E->lmesh.nel; e++)
        heating_latent[e] = 1.0;

    if(E->control.Ra_410 != 0.0) {
        latent_heating(E, m, heating_latent, heating_adi,
                       E->Fas410, E->control.Ra_410,
                       E->control.clapeyron410, E->viscosity.z410,
                       E->control.transT410, E->control.inv_width410);

    }

    if(E->control.Ra_670 != 0.0) {
        latent_heating(E, m, heating_latent, heating_adi,
                       E->Fas670, E->control.Ra_670,
                       E->control.clapeyron670, E->viscosity.zlm,
                       E->control.transT670, E->control.inv_width670);
    }

    if(E->control.Ra_cmb != 0.0) {
        latent_heating(E, m, heating_latent, heating_adi,
                       E->Fascmb, E->control.Ra_cmb,
                       E->control.clapeyroncmb, E->viscosity.zcmb,
                       E->control.transTcmb, E->control.inv_widthcmb);
    }


    if(E->control.Ra_410 != 0 || E->control.Ra_670 != 0.0 ||
       E->control.Ra_cmb != 0) {
        for(e=1; e<=E->lmesh.nel; e++)
            heating_latent[e] = 1.0 / heating_latent[e];
    }

    return;
}


static double total_heating(struct All_variables *E, double **heating)
{
    int m, e;
    double sum, total;

    /* sum up within each processor */
    sum = 0;
    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(e=1; e<=E->lmesh.nel; e++)
            sum += heating[m][e] * E->eco[m][e].area;
    }

    /* sum up for all processors */
    MPI_Allreduce(&sum, &total, 1,
                  MPI_DOUBLE, MPI_SUM, E->parallel.world);

    return total;
}


static void process_heating(struct All_variables *E, int psc_pass)
{
    int m;
    double total_visc_heating, total_adi_heating;

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        if(psc_pass == 0) {
            /* visc heating does not change between psc_pass, compute only
             * at first psc_pass */
            process_visc_heating(E, m, E->heating_visc[m]);
        }
        process_adi_heating(E, m, E->heating_adi[m]);
        process_latent_heating(E, m, E->heating_latent[m], E->heating_adi[m]);
    }

    /* compute total amount of visc/adi heating over all processors
     * only at last psc_pass */
    if(psc_pass == (E->advection.temp_iterations-1)) {
        total_visc_heating = total_heating(E, E->heating_visc);
        total_adi_heating = total_heating(E, E->heating_adi);

        if(E->parallel.me == 0) {
            fprintf(E->fp, "Step: %d, Total_heating(visc, adi): %g %g\n",
                    E->monitor.solution_cycles,
                    total_visc_heating, total_adi_heating);
            fprintf(stderr, "Step: %d, Total_heating(visc, adi): %g %g\n",
                    E->monitor.solution_cycles,
                    total_visc_heating, total_adi_heating);
        }
    }

    return;
}

void isotope_input(struct All_variables *E)
{
	int m = E->parallel.me;
	int n;
	input_int("isotope_n",&(E->convection.heat_sources.number),"0",m);
	input_float("isotope_t",&(E->convection.heat_sources.t_offset),"0.0",m);
	n = E->convection.heat_sources.number;
	if (n >0) {
		E->convection.heat_sources.Q = (float*)malloc(n*sizeof(float));
		E->convection.heat_sources.QER = (float*)malloc(n*sizeof(float));
		E->convection.heat_sources.QER1 = (float*)malloc(n*sizeof(float));
		E->convection.heat_sources.lambda = (float*)malloc(n*sizeof(float));
		if ( !input_float_vector("isotope_Q",n,E->convection.heat_sources.Q,m))
		{
		fprintf(stderr,"Missing input parameter: 'isotope_Q'\n"); 
		parallel_process_termination();
		}
		if ( E->control.tracer_enriched)
		if ( !input_float_vector("Q0_enriched",n,E->convection.heat_sources.QER,m))
		{
		fprintf(stderr,"Missing input parameter: 'Q0_enriched'\n"); 
		parallel_process_termination();
		}
		input_int("enriched_flavor1",&(E->trace.enriched_flavor1),"0",m);
		if (E->trace.enriched_flavor1>0)
		if ( !input_float_vector("Q0_enriched1",n,E->convection.heat_sources.QER1,m))
		{
		fprintf(stderr,"Missing input parameter: 'Q0_enriched1'\n"); 
		parallel_process_termination();
		}
		if (! input_float_vector("isotope_lambda",n,E->convection.heat_sources.lambda,m))
		{
		fprintf(stderr,"Missing input parameter: 'isotope_lambda'\n"); 
		parallel_process_termination();
		}
	}
}
