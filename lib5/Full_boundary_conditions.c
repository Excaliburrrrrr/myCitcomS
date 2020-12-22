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
#include "element_definitions.h"
#include "global_defs.h"
#include <math.h>

#include "lith_age.h"
#ifdef USE_GGRD
#include "ggrd_handling.h"
#endif
/* ========================================== */

static void horizontal_bc(struct All_variables *,float *[],int,int,float,unsigned int,char,int,int);
void assign_internal_bc(struct All_variables * );
static void velocity_apply_periodic_bcs();
static void temperature_apply_periodic_bcs();
static void read_update_cmb_T_control(struct All_variables *);
static void read_core_liquidus(struct All_variables *);
void read_temperature_boundary_from_file(struct All_variables *);
void read_velocity_boundary_from_file(struct All_variables *);

/* ========================================== */

void full_velocity_boundary_conditions(E)
     struct All_variables *E;
{
  void velocity_imp_vert_bc();
  void velocity_apply_periodicapply_periodic_bcs();

  void apply_side_sbc();

  int j,noz,lv,k,node;

  for(lv=E->mesh.gridmax;lv>=E->mesh.gridmin;lv--)
    for (j=1;j<=E->sphere.caps_per_proc;j++)     {
      noz = E->mesh.NOZ[lv];
      if(E->mesh.topvbc != 1) {	/* free slip top */
	horizontal_bc(E,E->sphere.cap[j].VB,noz,1,0.0,VBX,0,lv,j);
	horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,VBZ,1,lv,j);
	horizontal_bc(E,E->sphere.cap[j].VB,noz,2,0.0,VBY,0,lv,j);
	horizontal_bc(E,E->sphere.cap[j].VB,noz,1,E->control.VBXtopval,SBX,1,lv,j);
	horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,SBZ,0,lv,j);
	horizontal_bc(E,E->sphere.cap[j].VB,noz,2,E->control.VBYtopval,SBY,1,lv,j);
#ifdef USE_GGRD
	/* Ggrd traction control */
	if((lv==E->mesh.gridmax) && E->control.ggrd.vtop_control)
	  ggrd_read_vtop_from_file(E, TRUE);
#endif

      }
      if(E->mesh.botvbc != 1) {	/* free slip bottom */
        horizontal_bc(E,E->sphere.cap[j].VB,1,1,0.0,VBX,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,3,0.0,VBZ,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,2,0.0,VBY,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,1,E->control.VBXbotval,SBX,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,3,0.0,SBZ,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,2,E->control.VBYbotval,SBY,1,lv,j);
        }

      if(E->mesh.topvbc == 1) {	/* velocity/no slip BC */
        horizontal_bc(E,E->sphere.cap[j].VB,noz,1,E->control.VBXtopval,VBX,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,VBZ,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,2,E->control.VBYtopval,VBY,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,1,0.0,SBX,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,SBZ,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,2,0.0,SBY,0,lv,j);

#ifdef USE_GGRD
	/* Ggrd velocity control */
	if((lv==E->mesh.gridmax) && E->control.ggrd.vtop_control)
	  ggrd_read_vtop_from_file(E,TRUE);
#endif


        if(E->control.vbcs_file){ /* this should either only be called
				     once, or the input routines need
				     to be told what to do for each
				     multigrid level and cap. it might
				     be easiest to call only once and
				     have routines deal with multigrid
				  */
	  if((lv == E->mesh.gridmin) && (j == E->sphere.caps_per_proc))
	     read_velocity_boundary_from_file(E);
	}
      }

      if(E->mesh.botvbc == 1) {	/* velocity bottom BC */
        horizontal_bc(E,E->sphere.cap[j].VB,1,1,E->control.VBXbotval,VBX,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,3,0.0,VBZ,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,2,E->control.VBYbotval,VBY,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,1,0.0,SBX,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,3,0.0,SBZ,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,2,0.0,SBY,0,lv,j);
        }
      }    /* end for j and lv */

      if(E->control.side_sbcs)
	apply_side_sbc(E);

/* if(E->control.verbose) { */
/*  for (j=1;j<=E->sphere.caps_per_proc;j++) */
/*    for (node=1;node<=E->lmesh.nno;node++) */
/*       fprintf(E->fp_out,"m=%d VB== %d %g %g %g flag %u %u %u\n",j,node,E->sphere.cap[j].VB[1][node],E->sphere.cap[j].VB[2][node],E->sphere.cap[j].VB[3][node],E->node[j][node]&VBX,E->node[j][node]&VBY,E->node[j][node]&VBZ); */
/*  fflush(E->fp_out); */
/* } */

  /* If any imposed internal velocity structure it goes here */

      
      /*
	apply stress or velocity boundary conditions, read from file
	settings are to be implemented in those routines (will only do
	anything at present, if E->mesh.toplayerbc != 0
      */
      assign_internal_bc(E);

   return; }

/* ========================================== */

void full_temperature_boundary_conditions(E)
     struct All_variables *E;
{
  void temperatures_conform_bcs();
  void temperature_imposed_vert_bcs();
  int j,lev,noz;

  lev = E->mesh.levmax;
  for (j=1;j<=E->sphere.caps_per_proc;j++)    {
    noz = E->mesh.noz;
    if(E->mesh.toptbc == 1)    {
      horizontal_bc(E,E->sphere.cap[j].TB,noz,3,E->control.TBCtopval,TBZ,1,lev,j);
      horizontal_bc(E,E->sphere.cap[j].TB,noz,3,E->control.TBCtopval,FBZ,0,lev,j);
      if(E->control.tbcs_file)
          read_temperature_boundary_from_file(E);
      }
    else   {
      horizontal_bc(E,E->sphere.cap[j].TB,noz,3,E->control.TBCtopval,TBZ,0,lev,j);
      horizontal_bc(E,E->sphere.cap[j].TB,noz,3,E->control.TBCtopval,FBZ,1,lev,j);
      }

    if(E->mesh.bottbc == 1)    {
      horizontal_bc(E,E->sphere.cap[j].TB,1,3,E->control.TBCbotval,TBZ,1,lev,j);
      horizontal_bc(E,E->sphere.cap[j].TB,1,3,E->control.TBCbotval,FBZ,0,lev,j);
      }
    else        {
      horizontal_bc(E,E->sphere.cap[j].TB,1,3,E->control.TBCbotval,TBZ,0,lev,j);
      horizontal_bc(E,E->sphere.cap[j].TB,1,3,E->control.TBCbotval,FBZ,1,lev,j);
      }

    if(E->control.lith_age_time==1)  {

   /* set the regions in which to use lithosphere files to determine temperature
   note that this is called if the lithosphere age in inputted every time step
   OR it is only maintained in the boundary regions */
      lith_age_temperature_bound_adj(E,lev);
    }


    }     /* end for j */

  temperatures_conform_bcs(E);
  E->temperatures_conform_bcs = temperatures_conform_bcs;

   return; }


/*  =========================================================  */

static void horizontal_bc(struct All_variables *E,float *BC[],int ROW,int dirn,float value,
			  unsigned int mask,char onoff,int level,int m)
{
  int i,j,node,rowl;

    /* safety feature */
  if(dirn > E->mesh.nsd)
     return;

  if (ROW==1)
      rowl = 1;
  else
      rowl = E->lmesh.NOZ[level];

  if ( ( (ROW==1) && (E->parallel.me_loc[3]==0) ) ||
       ( (ROW==E->mesh.NOZ[level]) && (E->parallel.me_loc[3]==E->parallel.nprocz-1) ) ) {

    /* turn bc marker to zero */
    if (onoff == 0)          {
      for(j=1;j<=E->lmesh.NOY[level];j++)
    	for(i=1;i<=E->lmesh.NOX[level];i++)     {
    	  node = rowl+(i-1)*E->lmesh.NOZ[level]+(j-1)*E->lmesh.NOX[level]*E->lmesh.NOZ[level];
    	  E->NODE[level][m][node] = E->NODE[level][m][node] & (~ mask);
    	  }        /* end for loop i & j */
      }

    /* turn bc marker to one */
    else        {
      for(j=1;j<=E->lmesh.NOY[level];j++)
        for(i=1;i<=E->lmesh.NOX[level];i++)       {
    	  node = rowl+(i-1)*E->lmesh.NOZ[level]+(j-1)*E->lmesh.NOX[level]*E->lmesh.NOZ[level];
    	  E->NODE[level][m][node] = E->NODE[level][m][node] | (mask);
    	  if(level==E->mesh.levmax)   /* NB */
    	    BC[dirn][node] = value;
    	  }     /* end for loop i & j */
      }

    }             /* end for if ROW */

  return;
}


static void velocity_apply_periodic_bcs(E)
    struct All_variables *E;
{
  fprintf(E->fp,"Periodic boundary conditions\n");

  return;
  }

static void temperature_apply_periodic_bcs(E)
    struct All_variables *E;
{
 fprintf(E->fp,"Periodic temperature boundary conditions\n");

  return;
  }



/*lhy update_cmb_T
 * read in controls when initially set the problem
 * when restarted, this should also be run
 */
void initial_update_cmb_T(struct All_variables *E){
	read_update_cmb_T_control(E);
	if(E->control.inner_core_latent){
		read_core_liquidus(E); //read core liquidus from file
	}
}
/*lhy update_cmb_T
 * read in controling parameters*/
static void read_update_cmb_T_control(struct All_variables *E){
	int m = E->parallel.me;
	/*method for updating cmb tempreture
	 * 0 : don't update
	 * 1 : due to heat flux on cmb and latent heat of solidifying inner core
	 */
	input_int("coreT_method",&(E->control.coreT),"0",m); 
	input_int("inner_core_latent",&(E->control.inner_core_latent),"0",m); //1: include inner core latent heat
	input_double("cpc",&(E->data.Cpc),"1000.0",m); //core heat capacity
	input_double("Lc",&(E->data.Lc),"300000.0",m); //latent heat
	input_double("inner_core_radius",&(E->core.rin),"0.0",m); //initial inner core undimentional radius, change when restarted
	input_double("Tcmb",&(E->core.Tcmb),"0.0",m); //initial cmb undimentional temperature, change when restarted
	input_int("cnr",&(E->core.cnr),"50",m); //number of radial point on the core liquidus profile
	input_string("core_liquidus_file",E->control.core_liquidus_file,"core_solidus",m); //file name for the core liquidus profile
	if(E->control.inner_core_latent){
		E->core.crr = (double*)malloc((E->core.cnr+1)*sizeof(double));
		E->core.liq = (double*)malloc((E->core.cnr+1)*sizeof(double));
		E->core.liq_diff = (double*)malloc((E->core.cnr+1)*sizeof(double));
	}
}
/*lhy update_cmb_T
 * read in core liquidus temperature gradient*/
static void read_core_liquidus(struct All_variables *E){
	int cnr,i;
	double v1,v2,v3;
	char input_file[255], input_s[1000];
	FILE *fp;
	cnr = E->core.cnr;
	sprintf(input_file,"%s/%s",E->control.sol_liq_dir,E->control.core_liquidus_file);
	fp=fopen(input_file,"r");
	if (fp == NULL) {
		fprintf(stderr,"(read_sol_liq_from_file.c #1) Cannot open %s\n",input_file);
		parallel_process_termination();
	}
	for (i=1;i<=cnr;i++){
		fgets(input_s,1000,fp);
		/*v1: undimentional radius. v2: undimentional temperature. v3: undimentional temperature gradient to radius*/
		if(sscanf(input_s,"%lf %lf %lf",&(v1),&(v2),&(v3)) != 3){
			fprintf(stderr,"Error while reading file '%s'\n", input_file);
			exit(8);
		}
		E->core.crr[i] = v1;
		E->core.liq[i] = v2;
		E->core.liq_diff[i] = v3;
	}
	fclose(fp);
}
/*lhy update_cmb_T
 * get core liquidus temperature gradient
 * linearly derive value between two existing values
 * */
double get_diff_liq(struct All_variables *E,double rin){
	int cnr, i;
	double diff_liq;
	double *crr, *liqArr;
	
	cnr=E->core.cnr;
	rin=E->core.rin;
	crr=E->core.crr;
	liqArr=E->core.liq_diff;
	
	if(E->parallel.me == 0)
		fprintf(stderr,"rin: %.4e\n",rin); //lhy debug
	for(i=1;i<cnr;i++){
		if((rin > crr[i]-1e-5)&&(rin < crr[i+1])){
			diff_liq = (rin-crr[i+1])/(crr[i]-crr[i+1])*liqArr[i] +\
					(crr[i]-rin)/(crr[i]-crr[i+1])*liqArr[i+1];
			break;
		}
	}
	return diff_liq;
}
/*lhy update_cmb_T
 * update cmb temperature by average heat flux on cmb
 * can also include: 1. effect of temperature profile in outer core
 * 				     2. latent heat release by solidification of inner core
 * this function is called in main function
 */
void update_cmb_T(struct All_variables *E){
	double rhom, rhoc, Cm, Cc, rc, refT, deltat, deltar, gradT, deltaT;
	double L, Rc, rin, diffTc, diff_liq, coeffQ, coeffL;
	int IsEutex,j,lev;
  	void temperatures_conform_bcs();	
	double get_diff_liq();
	if (E->parallel.me == 0)
		fprintf(stderr,"update_cmb_T\n"); 

	lev = E->mesh.levmax;
	rhom = E->data.density;
	rhoc = E->data.density_below;
	Cm = E->data.Cp;
	Cc = E->data.Cpc;
	rc = E->sphere.ri; 
	L = E->data.Lc;
	refT = E->data.ref_temperature;
	Rc = rc*E->data.radius;
	gradT = E->sphere.cap[1].heat_flux;

	deltat = E->advection.timestep;

	IsEutex = 0;
	if (E->parallel.me == 0)
		fprintf(stderr,"gradT = %.4e, deltat = %.4e\n", gradT, deltat); //lhy debug
	/*see nanzhang etal 2013 for reference
	 * undimentionalize all components in eq 12
	 */
	coeffQ = 3*rhom*Cm/(rhoc*Cc)/rc*gradT*deltat;

	if (E->control.coreT == 1)
		diffTc = 1;
	else if (E->control.coreT == 2){
	}
	
	/*inner core latent heat*/
	if ((E->control.inner_core_latent)&&(E->core.Tcmb<E->core.liq[1])){
		rin = E->core.rin;
		diff_liq=get_diff_liq(E,rin);
		if (E->parallel.me == 0)
			fprintf(stderr,"diffliq: %.4e\n",diff_liq); //lhy debug
		/*Eutectic point is point in phase transform diagram where T approximately remain contant and S begin to merge in solid*/
		if (fabs(diff_liq) < 1e-5)
			IsEutex = 1;
		else
			coeffL = 3*L/(refT*Cc)*pow(rin,2)/pow(rc,3)/diff_liq;
	}
	else
		coeffL = 0;
	
	if (IsEutex == 1)
		deltaT = 0.0;
	else
		deltaT = (-1.0)*coeffQ/(diffTc-coeffL);
	
	if ((E->control.inner_core_latent)&&(E->core.Tcmb<E->core.liq[1])){
		if (IsEutex == 1)
			deltar = rhom*Cm*refT/rhoc/L*pow(rc/rin,2)*gradT*deltat;
	 	else
			deltar = deltaT/diff_liq;
	}
	else
		deltar = 0;
	
	E->core.Tcmb += deltaT;
	E->core.rin += deltar;
	E->control.TBCbotval = E->core.Tcmb; //may seems redundant, but core.Tmb serves to save initial cmb temperature. 

	if (E->parallel.me == 0)
		fprintf(stderr,"deltaT = %.4e, deltar = %.4e\n",deltaT,deltar); //lhy debug
	/*change and conform boundary condition with updated cmb temperature*/
	for(j=1; j<=E->sphere.caps_per_proc; j++){
		if(E->mesh.bottbc == 1){
	  		if (E->parallel.me == 0) //lhy debug
		  		fprintf(stderr,"TBCbotval = %.4e, rin = %.4e\n",E->control.TBCbotval,E->core.rin);
      	horizontal_bc(E,E->sphere.cap[j].TB,1,3,E->control.TBCbotval,TBZ,1,lev,j);
      	horizontal_bc(E,E->sphere.cap[j].TB,1,3,E->control.TBCbotval,FBZ,0,lev,j);		
		}
	}
  	temperatures_conform_bcs(E);	
}


/* version */
/* $Id$ */

/* End of file  */
