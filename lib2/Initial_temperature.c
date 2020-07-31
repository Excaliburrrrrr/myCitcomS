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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA *
 *</LicenseText>
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <mpi.h>
#include <math.h>
#include <sys/types.h>

#include "global_defs.h"
#include "lith_age.h"
#include "parsing.h"
#include "element_definitions.h"
#include "global_defs.h"
#include "parsing.h"
#include "output.h"

void parallel_process_termination();
void sphere_expansion();
void temperatures_conform_bcs(struct All_variables *);
double modified_plgndr_a(int, int, double);
void rtp2xyzd(double,double,double,double *);
void sphere_expansion();
void correct_for_mc(struct All_variables *,float **,double *);

#include "initial_temperature.h"
static void debug_tic(struct All_variables *);
static void read_tic_from_file(struct All_variables *);
static void read_tic_from_file_im(struct All_variables *);
static void construct_tic_from_input(struct All_variables *);
static void add_perturbations_at_all_layers(struct All_variables *);

#ifdef USE_GZDIR
void restart_tic_from_gzdir_file(struct All_variables *);
#endif
#ifdef USE_GGRD
#include "ggrd_handling.h"
#endif
static double ibc_layer(struct All_variables *, double);
static void all_minerals(struct All_variables *);
static void all_chemicals(struct All_variables *);
static double none_dimensional_parameters(struct All_variables *);
static void write_specified_model_parameters(struct All_variables *, char*);
static void lunar_read_parametes(struct All_variables *);
static double T_dependent_viscosity(struct All_variables *, double, double, double, double, double);
static void write_double_array_into_file(FILE *, char *, double *, int, int);
static double visc_from_chemical(struct All_variables *, double *, double *);
static double iso_stress(double, double, double);
static void read_temp_from_radial(struct All_variables *);
static double ibc_layer(struct All_variables *, double);
static void all_minerals(struct All_variables *);
static void all_chemicals(struct All_variables *);
static double none_dimensional_parameters(struct All_variables *);
static void write_specified_model_parameters(struct All_variables *, char*);
static void lunar_read_parametes(struct All_variables *);
static double T_dependent_viscosity(struct All_variables *, double, double, double, double, double);
static void write_double_array_into_file(FILE *, char *, double *, int, int);
static double T_dependent_viscosity(struct All_variables *, double, double, double, double, double);
void compute_horiz_avg(struct All_variables *);
static reset_IBC(struct All_variables *);
void count_tracers_of_flavors(struct All_variables *);
void fill_composition(struct All_variables *);
void cart_to_sphere(struct All_variables *, double, double, double, double *, double *, double *);
void element_coordinate(struct All_variables *, int, double*, double*);
void read_sol_liq_from_file(struct All_variables *);
void read_sol_liq_from_file1(struct All_variables *);

void tic_input(struct All_variables *E)
{
  void melting_initialize();
  int m = E->parallel.me;
  int noz = E->lmesh.noz;
  int n;
#ifdef USE_GGRD
  int tmp;
#endif

  input_int("tic_method", &(E->convection.tic_method), "0,0,2", m);

  input_boolean("tic_max_T", &(E->convection.tic_max_T), "off", m);

  input_int("sol_liq", &(E->convection.sol_liq), "0,0,2", m); //lhy solidus_liquidus

  input_float("tic_max_T_value", &(E->convection.tic_max_T_value), "1.0", m);

#ifdef USE_GGRD			/* for backward capability */
  input_int("ggrd_tinit", &tmp, "0", m);
  if(tmp){
    E->convection.tic_method = 4; /*  */
    E->control.ggrd.use_temp = 1;
  }
#endif
  /* When tic_method is 0 (default), the temperature is a linear profile +
     perturbation at some layers.

     When tic_method is -1, the temperature is read in from the
     [datafile_old].velo.[rank].[solution_cycles_init] files.

     When tic_method is 1, the temperature is isothermal (== bottom b.c.) +
     uniformly cold plate (thickness specified by 'half_space_age').

     When tic_method is 2, (tic_method==1) + a hot blob. A user can specify
     the location and radius of the blob, and also the amplitude of temperature
     change in the blob relative to the ambient mantle temperautre
     (E->control.mantle_temp).
        - blob_center: A comma-separated list of three float numbers.
        - blob_radius: A dmensionless length, typically a fraction
                       of the Earth's radius.
        - blob_dT    : Dimensionless temperature.

     When tic_method is 3, the temperature is a linear profile + perturbation
     for whole mantle.

     tic_method is 4: read in initial temperature distribution from a set of netcdf grd
                      files. this required the GGRD extension to be compiled in

  */

    /* This part put a temperature anomaly at depth where the global
       node number is equal to load_depth. The horizontal pattern of
       the anomaly is given by spherical harmonic ll & mm. */

    input_int("num_perturbations", &n, "0,0,PERTURB_MAX_LAYERS", m);

    if (n > 0) {
      E->convection.number_of_perturbations = n;

      if (! input_float_vector("perturbmag", n, E->convection.perturb_mag, m) ) {
	fprintf(stderr,"Missing input parameter: 'perturbmag'\n");
	parallel_process_termination();
      }
      if (! input_int_vector("perturbm", n, E->convection.perturb_mm, m) ) {
	fprintf(stderr,"Missing input parameter: 'perturbm'\n");
	parallel_process_termination();
      }
      if (! input_int_vector("perturbl", n, E->convection.perturb_ll, m) ) {
	fprintf(stderr,"Missing input parameter: 'perturbl'\n");
	parallel_process_termination();
      }
      if (! input_int_vector("perturblayer", n, E->convection.load_depth, m) ) {
	fprintf(stderr,"Missing input parameter: 'perturblayer'\n");
	parallel_process_termination();
      }
    }
    else {
      E->convection.number_of_perturbations = 1;
      E->convection.perturb_mag[0] = 1;
      E->convection.perturb_mm[0] = 2;
      E->convection.perturb_ll[0] = 2;
      E->convection.load_depth[0] = (noz+1)/2;
    }

    input_float("half_space_age", &(E->convection.half_space_age), "40.0,1e-3,nomax", m);
    input_float("mantle_temp",&(E->control.mantle_temp),"1.0",m);


    switch(E->convection.tic_method){
    case 2:			/* blob */
      if( ! input_float_vector("blob_center", 3, E->convection.blob_center, m)) {
	assert( E->sphere.caps == 12 || E->sphere.caps == 1 );
	if(E->sphere.caps == 12) { /* Full version: just quit here */
	  fprintf(stderr,"Missing input parameter: 'blob_center'.\n");
	  parallel_process_termination();
	}
	else if(E->sphere.caps == 1) { /* Regional version: put the blob at the center */
	  fprintf(stderr,"Missing input parameter: 'blob_center'. The blob will be placed at the center of the domain.\n");
	  E->convection.blob_center[0] = 0.5*(E->control.theta_min+E->control.theta_max);
	  E->convection.blob_center[1] = 0.5*(E->control.fi_min+E->control.fi_max);
	  E->convection.blob_center[2] = 0.5*(E->sphere.ri+E->sphere.ro);
	}
      }
      input_float("blob_radius", &(E->convection.blob_radius), "0.063,0.0,1.0", m);
      input_float("blob_dT", &(E->convection.blob_dT), "0.18,nomin,nomax", m);
      input_boolean("blob_bc_persist",&(E->convection.blob_bc_persist),"off",m);
      break;
    case 4:
      /*
	case 4: initial temp from grd files
      */
#ifdef USE_GGRD
      /*
	 read in some more parameters

      */
      /* scale the anomalies with PREM densities */
      input_boolean("ggrd_tinit_scale_with_prem",
		    &(E->control.ggrd.temp.scale_with_prem),"off",E->parallel.me);
      /* limit T to 0...1 */
      input_boolean("ggrd_tinit_limit_trange",
		    &(E->control.ggrd.temp.limit_trange),"on",E->parallel.me);
      /* scaling factor for the grids */
      input_double("ggrd_tinit_scale",
		   &(E->control.ggrd.temp.scale),"1.0",E->parallel.me); /* scale */
      /* temperature offset factor */
      input_double("ggrd_tinit_offset",
		   &(E->control.ggrd.temp.offset),"0.0",E->parallel.me); /* offset */
      /*
	 do we want a different scaling for the lower mantle?
      */
      input_float("ggrd_lower_depth_km",&(E->control.ggrd_lower_depth_km),"7000",
		  E->parallel.me); /* depth, in km, below which
				      different scaling applies */
      input_float("ggrd_lower_scale",&(E->control.ggrd_lower_scale),"1.0",E->parallel.me);
      input_float("ggrd_lower_offset",&(E->control.ggrd_lower_offset),"1.0",E->parallel.me);

      /* grid name, without the .i.grd suffix */
      input_string("ggrd_tinit_gfile",
		   E->control.ggrd.temp.gfile,"",E->parallel.me); /* grids */
      input_string("ggrd_tinit_dfile",
		   E->control.ggrd.temp.dfile,"",E->parallel.me); /* depth.dat layers of grids*/
      /* override temperature boundary condition? */
      input_boolean("ggrd_tinit_override_tbc",
		    &(E->control.ggrd.temp.override_tbc),"off",E->parallel.me);
      input_string("ggrd_tinit_prem_file",
		   E->control.ggrd.temp.prem.model_filename,"hc/prem/prem.dat",
		   E->parallel.me); /* PREM model filename */

      /* non-linear scaling, downweighing negative anomalies? */
      input_boolean("ggrd_tinit_nl_scale",&(E->control.ggrd_tinit_nl_scale),"off",E->parallel.me);

#else
      fprintf(stderr,"tic_method 4 only works for USE_GGRD compiled code\n");
      parallel_process_termination();
#endif
      break;
    } /* no default needed */
    return;
}



void convection_initial_temperature(struct All_variables *E)
{
  int i,m;
  void report();
  void add_T_impact();
  void melting_initialize();

  report(E,"Initialize temperature field");

  if (E->convection.tic_method == -1) {
      /* read temperature from file */
#ifdef USE_GZDIR
      if(strcmp(E->output.format, "ascii-gz") == 0)
          restart_tic_from_gzdir_file(E);
      else
#endif
		  read_temp_from_radial(E);
 //         read_tic_from_file(E);
  		add_perturbations_at_layers(E);  //zwb 0515
	}
  else if (E->convection.tic_method == -2) {
      /* read temperature from file */

#ifdef USE_GZDIR
      if(strcmp(E->output.format, "ascii-gz") == 0)
          restart_tic_from_gzdir_file(E);
      else
#endif

          read_tic_from_file(E);
      add_perturbations_at_all_layers(E);
  }
  else if (E->convection.tic_method == -3) {
      /* read temperature from file */
		  read_temp_from_radial(E);
          //read_tic_from_file(E);
		  add_T_impact(E);
  }
  else if (E->convection.tic_method == -4) {
      /* read temperature from file */
#ifdef USE_GZDIR
      if(strcmp(E->output.format, "ascii-gz") == 0)
          restart_tic_from_gzdir_file(E);
      else
#endif
		  read_temp_from_radial(E);
          //read_tic_from_file(E);
      read_tic_from_file_im(E);
	  if(E->control.T_IBC_reset>0.0)
		  reset_IBC(E);
  }
  else if (E->convection.tic_method == -5) {
		  add_T_impact(E);

    /*for (m = 1; m <= E->sphere.caps_per_proc; m++)
	for (i = 1;i <= E->lmesh.nno; i++)
	{
	      E->T[m][i] = 0.0*E->T[m][i]; //changed
	}*/
          //read_tic_from_file_im(E);
  }
  else if (E->control.lith_age)
      lith_age_construct_tic(E);
  else
      construct_tic_from_input(E);

  /* Note: it is the callee's responsibility to conform tbc. */
  /* like a call to temperatures_conform_bcs(E); */

  if (E->convection.tic_max_T == 1)
	  for(m=1;m<=E->sphere.caps_per_proc;m++)
	  for (i=1;i<=E->lmesh.nno;i++)
		  if (E->T[m][i]>E->convection.tic_max_T_value)
	      E->T[m][i] = E->convection.tic_max_T_value;

  if (E->control.verbose)
    debug_tic(E);

  if (E->convection.sol_liq){
      melting_initialize(E);
	  read_sol_liq_from_file(E);
  }



  return;
}


static void debug_tic(struct All_variables *E)
{
  int m, j;

  fprintf(E->fp_out,"output_temperature\n");
  for(m=1;m<=E->sphere.caps_per_proc;m++)        {
    fprintf(E->fp_out,"for cap %d\n",E->sphere.capid[m]);
    for (j=1;j<=E->lmesh.nno;j++)
      fprintf(E->fp_out,"X = %.6e Z = %.6e Y = %.6e T[%06d] = %.6e \n",E->sx[m][1][j],E->sx[m][2][j],E->sx[m][3][j],j,E->T[m][j]);
  }
  fflush(E->fp_out);

  return;
}


static void read_tic_from_file(struct All_variables *E)
{
  int ii, ll, mm;
  float tt;
  int i, m;
  char output_file[255], input_s[1000];
  FILE *fp;

  float v1, v2, v3, g;

  ii = E->monitor.solution_cycles_init;
  sprintf(output_file,"%s.velo.%d.%d",E->control.old_P_file,E->parallel.me,ii);
  fp=fopen(output_file,"r");
  if (fp == NULL) {
    fprintf(stderr,"(Initial_temperature.c #1) Cannot open %s\n",output_file);
    parallel_process_termination();
  }



 if (E->parallel.me==0)
    fprintf(E->fp,"Reading %s for initial temperature\n",output_file);

  fgets(input_s,1000,fp);
  sscanf(input_s,"%d %d %f",&ll,&mm,&tt);

  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    fgets(input_s,1000,fp);
    sscanf(input_s,"%d %d",&ll,&mm);
    for(i=1;i<=E->lmesh.nno;i++)  {
      fgets(input_s,1000,fp);
      if(sscanf(input_s,"%g %g %g %f",&(v1),&(v2),&(v3),&(g)) != 4) {
        fprintf(stderr,"Error while reading file '%s'\n", output_file);
        exit(8);
      }
      /* Truncate the temperature to be within (0,1). */
      /* This might not be desirable in some situations. */
      E->T[m][i] = max(0.0,min(g,1.0));
    }
  }
  fclose (fp);

  temperatures_conform_bcs(E);

  return;
}

static void read_tic_from_file_im(struct All_variables *E)
{
  const int Is0 = (E->parallel.me==0);
  int ii, ll, mm;
  float tt;
  int i, m, j;
  int ncomp = 4;
  char output_file[255], input_s[1000];
  FILE *fp;

  double v1, v2, v3, T, g;

  const float Tref = E->data.ref_temperature;

  ii = E->monitor.solution_cycles_init;
  sprintf(output_file,"%s.coord.%d",E->control.old_P_file,E->parallel.me);
  fp=fopen(output_file,"r");
  if (fp == NULL) {
    fprintf(stderr,"(Initial_temperature.c #1) Cannot open %s\n",output_file);
    parallel_process_termination();
  }
  if (E->parallel.me==0)
    fprintf(E->fp,"Reading %s for initial temperature\n",output_file);

  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    for(i=1;i<=E->lmesh.nno;i++)  {
      fgets(input_s,1000,fp);
      if(sscanf(input_s,"%lf %lf %lf %lf",&(v1),&(v2),&(v3),&(T)) != 4) {
        fprintf(stderr,"Error while reading file '%s'\n", output_file);
        exit(8);
      }
	  T = T/Tref;
      /* Truncate the temperature to be within (0,1). */
      /* This might not be desirable in some situations. */
      //if (E->sx[m][3][i]<(1-E->viscosity.zlith))  //nanzhang2019
        E->T[m][i] += T;
      //E->T[m][i] = max(0.0,min(T,1.0));

    }
  }
  if(E->parallel.me==0){
	  for(i=1;i<=10;i++){
	  	fprintf(stderr,"T[%d] = %.4e\n",i,E->T[1][i]);
	  }//debug
  }
  fclose (fp);

  temperatures_conform_bcs(E);

  return;
}
/*lhy solidus_liquidus
 * allocate variables for melting calculation
 * two arrays for saving solidus and liquidus
 * one array for saving melting percent
 */
void melting_initialize(struct All_variables *E)
{
  int m,i;
  for(m=1;m<=E->sphere.caps_per_proc;m++) {
	  E->sol[m]=(double*)malloc((E->mesh.noz+1)*sizeof(double));
	  E->liq[m]=(double*)malloc((E->mesh.noz+1)*sizeof(double));
	  E->sol_l[m]=(double*)malloc((E->lmesh.noz+1)*sizeof(double));
	  E->liq_l[m]=(double*)malloc((E->lmesh.noz+1)*sizeof(double));
	  E->gsol[m]=(double*)malloc((E->mesh.noz+1)*sizeof(double));
  }

  for(m=1;m<=E->sphere.caps_per_proc;m++){
	  E->melting[m] = (double*)malloc((E->lmesh.nno+1)*sizeof(double));
	  E->melt_el[m] = (double*)malloc((E->lmesh.nel+1)*sizeof(double));
  }
  for(i = 0;i <= E->lmesh.nno;i++)	E->melting[1][i] = 0.0;
  for(i = 0;i <= E->lmesh.nel;i++)	E->melt_el[1][i] = 0.0;
  E->comp_melting_volume = (double*)malloc(E->composition.ncomp*sizeof(double));


}

/*lhy solidus_liqudus
 * read solidus and liqudus from file
 * assing value to two array of size E->mesh.noz + 1 (not using index 0)
 */
void read_sol_liq_from_file(struct All_variables *E)
{
  int ii, ll, mm;
  float tt;
  int i, m, noz, nz;
  char output_file[255], input_s[1000];
  FILE *fp;

  float v1, v2, v3, v4, g;

  ii = E->monitor.solution_cycles_init;
  noz = E->lmesh.noz;
  sprintf(output_file,"%s/%s",E->control.sol_liq_dir,E->control.sol_liq_file);
  fp=fopen(output_file,"r");
  if (fp == NULL) {
    fprintf(stderr,"(read_sol_liq_from_file.c #1) Cannot open %s\n",output_file);
    parallel_process_termination();
  }

  for(m=1;m<=E->sphere.caps_per_proc;m++) {
	if (E->parallel.me == 0)
		fprintf(stderr, "Solidus Liquidus Field: %s\n", output_file); //lhy debug
    for(i=1;i<=E->mesh.noz;i++)  {
      fgets(input_s,1000,fp);
	  if(E->lunar.latent_method == 2){
      	if(sscanf(input_s,"%g %g %g %g",&(v1),&(v2),&(v3),&(v4)) != 4) {
        	fprintf(stderr,"Error while reading file '%s'\n", output_file);
        	exit(8);
      	}
      	E->sol[m][i] = v2;
	  	E->liq[m][i] = v3;
      	E->gsol[m][i] = v4;
	  }
	  else{
      	if(sscanf(input_s,"%g %g %g",&(v1),&(v2),&(v3)) != 3) {
        	fprintf(stderr,"Error while reading file '%s'\n", output_file);
        	exit(8);
      	}
      	E->sol[m][i] = v2;
	  	E->liq[m][i] = v3;
	  }
  	}
    for(i=1;i<=noz;i++){
	nz = i+E->parallel.me_loc[3]*(noz-1);
	E->sol_l[m][i] = E->sol[m][nz];
	E->liq_l[m][i] = E->liq[m][nz];
    }
  }
	fclose(fp);
}

/*lhy solidus_liqudus
 * read solidus and liqudus from file
 * interpolate with respect to rr
 * assing value to two array of size E->mesh.noz + 1 (not using index 0)
 */
void read_sol_liq_from_file1(struct All_variables *E)
{
  const float tiny = 1e-7;
  const int size_of_input=1000;
  const int Is0 = (E->parallel.me == 0);
  int ii, ll, mm, num_input_nodes, input_element;
  float tt;
  float *rr;
  float *input_rad, *input_sol, *input_liq, *input_gsol;
  int i, j, m, noz, nz, keep_going, kk;
  char output_file[255], input_s[1000];
  FILE *fp;

  float v1, v2, v3, v4, g;
  float rad_bottom;
  float rad_top;
  float rad;
  float eta,delrad;
  float shape1,shape2;

  ii = E->monitor.solution_cycles_init;
  noz = E->lmesh.noz;
  sprintf(output_file,"%s/%s",E->control.sol_liq_dir,E->control.sol_liq_file);

  input_rad = (float*)malloc(size_of_input*sizeof(float));
  input_sol = (float*)malloc(size_of_input*sizeof(float));
  input_liq = (float*)malloc(size_of_input*sizeof(float));
  input_gsol = (float*)malloc(size_of_input*sizeof(float));
  fp=fopen(output_file,"r");
  if (fp == NULL) {
    fprintf(stderr,"(read_sol_liq_from_file.c #1) Cannot open %s\n",output_file);
    parallel_process_termination();
  }
  if(Is0)
      fprintf(stderr, "Read Solidus and Liquidus: %s\n",
              output_file);
  keep_going = 1;
  kk=0;
    while(keep_going)
    {
       if (fgets(input_s,200,fp)==NULL)
       {
          keep_going=0;
       }
       else
       {
          kk++;
          if (kk>(size_of_input-1))
          {
            fprintf(stderr,"ERROR(read_temp_from_radial) - file too big. Increase size_of_input\n");
            fflush(stderr);
            exit(10);
          }
          if(sscanf(input_s,"%f %f %f %f",&(v1),&(v2),&(v3),&(v4)) != 4) {
        	fprintf(stderr,"Error while reading file '%s'\n", output_file);
        	exit(8);
      	  }
          input_rad[kk] = v1;
          input_sol[kk] = v2;
          input_liq[kk] = v3;
          input_gsol[kk] = v4;
/* some control checks on input file */
          if ((kk>1)&&(input_rad[kk]<=input_rad[kk-1]))
          {
              fprintf(stderr,"ERROR(read_temp_from_radial)-rad does not increase? -check input file\n");
              fprintf(stderr,"rad[kk-1]: %f rad[kk]: %f kk: %d\n",input_rad[kk-1],input_rad[kk],kk);
              exit(10);
          }
       }
    }

    num_input_nodes=kk;
    if (num_input_nodes<2)
    {
       fprintf(stderr,"ERROR(read_temp_from_radial) - need at least 2 input points!\n");
       fflush(stderr);
       exit(10);
    }

/* interpolate citcom nz nodes */
    for (j=1;j<=E->sphere.caps_per_proc;j++)
    {

    for (kk=1;kk<=noz;kk++)
    {

       rad=E->sx[j][3][kk];
       if (rad>E->sphere.ro) rad = E->sphere.ro; //nanzhang
/* find which input element */

       input_element=0;
       for (mm=1;mm<=num_input_nodes-1;mm++)
       {
          rad_bottom=input_rad[mm];
          rad_top=input_rad[mm+1];
          if ( (rad>rad_bottom-tiny) && (rad<rad_top+tiny))
          {
             input_element=mm;
             goto foundit;
          }
       }

foundit:;

/* find local coordinate,eta, in input element. */
/* here, local coordinate extends from 0 to 1. */
       delrad=rad_top-rad_bottom;
       eta=(rad-rad_bottom)/delrad;

       if ((eta<-1e-6)||(eta>(1.0+1e-4)))
       {
          fprintf(stderr,"ERROR(read_temp_radial) z from m %d %d %f %f %f %f %f\n",kk,mm,eta,rad_bottom,rad_top,rad,delrad);
          fflush(stderr);
          exit(10);
       }

/* find shape functions at local coordinate, eta */
       shape1=input_sol[input_element]*(1.0-eta);
       shape2=input_sol[input_element+1]*(eta);
       E->sol[j][kk]=shape1+shape2;
       shape1=input_liq[input_element]*(1.0-eta);
       shape2=input_liq[input_element+1]*(eta);
       E->liq[j][kk]=shape1+shape2;
       shape1=input_gsol[input_element]*(1.0-eta);
       shape2=input_gsol[input_element+1]*(eta);
       E->gsol[j][kk]=shape1+shape2;
    }

    for(i=1;i<=noz;i++){
    	nz = i+E->parallel.me_loc[3]*(noz-1);
    	E->sol_l[j][i] = E->sol[j][nz];
    	E->liq_l[j][i] = E->liq[j][nz];
    }

    }
	fclose(fp);
    free(input_rad);
    free(input_sol);
    free(input_liq);
    free(input_gsol);
    return;

}

/*lhy solidus_liquidus
 * get melting percent from node temperature
 * using linear method fai = (T-Ts)/(Tl-Ts)
 * where fai is the melting percent in a solid matrix
 */
void get_melting_status(struct All_variables *E){
	int m, i, j, noz, lnz, nz;
	int procz, nprocz;
	double deltaT;
	double **comp_melt;
	procz = E->parallel.me_loc[3];
	nprocz = E->parallel.nprocz;
	noz = E->lmesh.noz;
	/*allocate space*/
	comp_melt=(double**)malloc((E->sphere.caps_per_proc+1)*sizeof(double*));
  	for(m=0;m<=E->sphere.caps_per_proc;m++) {
		comp_melt[m]=(double*)malloc((E->lmesh.nno+1)*sizeof(double));
	}


	if(E->parallel.me == 0)
		fprintf(stderr,"get_melting_status\n");
	/*calculate melt fraction on each node*/
  	for(m=1;m<=E->sphere.caps_per_proc;m++) {
    	for(i=1;i<=E->lmesh.nno;i++)  {
			lnz = i%noz;
			nz = ((lnz == 0)?noz:lnz) + procz * (noz - 1);
			deltaT = E->T[m][i] - E->sol[m][nz];
			E->melting[m][i] = ((deltaT < 0.0)? 0.0 : deltaT/(E->liq[m][nz] - E->sol[m][nz]));
			//if(E->parallel.me == 1)
				//fprintf(stderr,"%d %d %.4e %.4e %.4e\n",i,nz,E->sol[m][nz],E->liq[m][nz],E->melting[m][i]); //lhy debug
		}
  	}
	/*get total volume*/
	E->melting_volume = return_bulk_value_d(E,E->melting,0);
	/*map melting to each compositon*/
	for(j=0;j<E->composition.ncomp;j++){
  		for(m=1;m<=E->sphere.caps_per_proc;m++) {
    		for(i=1;i<=E->lmesh.nno;i++)  {
				comp_melt[m][i] = E->melting[m][i]*E->composition.comp_node[m][j][i];
			}
		}
		E->comp_melting_volume[j] = return_bulk_value_d(E,comp_melt,0);
	}
	/*free space*/
  	for(m=0;m<=E->sphere.caps_per_proc;m++) {
		free(comp_melt[m]);
	}
	free(comp_melt);
	return;
}

static void linear_temperature_profile(struct All_variables *E)
{
    int m, i, j, k, node;
    int nox, noy, noz;
    double r1;

    nox = E->lmesh.nox;
    noy = E->lmesh.noy;
    noz = E->lmesh.noz;

    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=1; i<=noy; i++)
            for(j=1; j<=nox;j ++)
                for(k=1; k<=noz; k++) {
                    node = k + (j-1)*noz + (i-1)*nox*noz;
                    r1 = E->sx[m][3][node];
                    E->T[m][node] = E->control.TBCbotval - (E->control.TBCtopval + E->control.TBCbotval)*(r1 - E->sphere.ri)/(E->sphere.ro - E->sphere.ri);
                }

    return;
}


static void conductive_temperature_profile(struct All_variables *E)
{
    int m, i, j, k, node;
    int nox, noy, noz;
    double r1;

    nox = E->lmesh.nox;
    noy = E->lmesh.noy;
    noz = E->lmesh.noz;

    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=1; i<=noy; i++)
            for(j=1; j<=nox;j ++)
                for(k=1; k<=noz; k++) {
                    node = k + (j-1)*noz + (i-1)*nox*noz;
                    r1 = E->sx[m][3][node];
                    E->T[m][node] = (E->control.TBCtopval*E->sphere.ro
                                     - E->control.TBCbotval*E->sphere.ri)
                        / (E->sphere.ro - E->sphere.ri)
                        + (E->control.TBCbotval - E->control.TBCtopval)
                        * E->sphere.ro * E->sphere.ri / r1
                        / (E->sphere.ro - E->sphere.ri);
                }

    return;
}


static void constant_temperature_profile(struct All_variables *E, double mantle_temp)
{
    int m, i;

    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=1; i<=E->lmesh.nno; i++)
            E->T[m][i] = mantle_temp;

    return;
}


static void add_top_tbl(struct All_variables *E, double age_in_myrs, double mantle_temp)
{
    int m, i, j, k, node;
    int nox, noy, noz;
    double r1, dT, tmp;

    nox = E->lmesh.nox;
    noy = E->lmesh.noy;
    noz = E->lmesh.noz;

    dT = (mantle_temp - E->control.TBCtopval);
    tmp = 0.5 / sqrt(age_in_myrs / E->data.scalet);

    fprintf(stderr, "%e %e\n", dT, tmp);
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=1; i<=noy; i++)
            for(j=1; j<=nox;j ++)
                for(k=1; k<=noz; k++) {
                    node = k + (j-1)*noz + (i-1)*nox*noz;
                    r1 = E->sx[m][3][node];
                    E->T[m][node] -= dT * erfc(tmp * (E->sphere.ro - r1));
                }

    return;
}


static void add_bottom_tbl(struct All_variables *E, double age_in_myrs, double mantle_temp)
{
    int m, i, j, k, node;
    int nox, noy, noz;
    double r1, dT, tmp;

    nox = E->lmesh.nox;
    noy = E->lmesh.noy;
    noz = E->lmesh.noz;

    dT = (E->control.TBCbotval - mantle_temp);
    tmp = 0.5 / sqrt(age_in_myrs / E->data.scalet);

    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=1; i<=noy; i++)
            for(j=1; j<=nox;j ++)
                for(k=1; k<=noz; k++) {
                    node = k + (j-1)*noz + (i-1)*nox*noz;
                    r1 = E->sx[m][3][node];
                    E->T[m][node] += dT * erfc(tmp * (r1 - E->sphere.ri));
                }

    return;
}


static void add_perturbations_at_layers(struct All_variables *E)
{
    /* This function put a temperature anomaly at depth where the global
       node number is equal to load_depth. The horizontal pattern of
       the anomaly is given by wavenumber (in regional model) or
       by spherical harmonic (in global model). */

    int m, i, j, k, node;
    int p, ll, mm, kk;
    int nox, noy, noz, gnoz;
    double t1, f1, tlen, flen, con;

    nox = E->lmesh.nox;
    noy = E->lmesh.noy;
    noz = E->lmesh.noz;
    gnoz = E->mesh.noz;

    for (p=0; p<E->convection.number_of_perturbations; p++) {
        ll = E->convection.perturb_ll[p];
        mm = E->convection.perturb_mm[p];
        kk = E->convection.load_depth[p];
        con = E->convection.perturb_mag[p];

        if ( (kk < 1) || (kk > gnoz) ) continue; /* layer kk is outside domain */

        k = kk - E->lmesh.nzs + 1; /* convert global nz to local nz */
        if ( (k < 1) || (k > noz) ) continue; /* layer k is not inside this proc. */
        if (E->parallel.me_loc[1] == 0 && E->parallel.me_loc[2] == 0
            && E->sphere.capid[1] == 1 )
            fprintf(stderr,"Initial temperature perturbation:  layer=%d  mag=%g  l=%d  m=%d\n", kk, con, ll, mm);

        if(E->sphere.caps == 1) {
            /* regional mode, add sinosoidal perturbation */

            tlen = M_PI / (E->control.theta_max - E->control.theta_min);
            flen = M_PI / (E->control.fi_max - E->control.fi_min);

            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(i=1; i<=noy; i++)
                    for(j=1; j<=nox;j ++) {
                        node = k + (j-1)*noz + (i-1)*nox*noz;
                        t1 = (E->sx[m][1][node] - E->control.theta_min) * tlen;
                        f1 = (E->sx[m][2][node] - E->control.fi_min) * flen;

                        E->T[m][node] += con * cos(ll*t1) * cos(mm*f1);
                    }
        }
        else {
            /* global mode, add spherical harmonics perturbation */

            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(i=1; i<=noy; i++)
                    for(j=1; j<=nox;j ++) {
                        node = k + (j-1)*noz + (i-1)*nox*noz;
                        t1 = E->sx[m][1][node];
                        f1 = E->sx[m][2][node];

                        E->T[m][node] += con * modified_plgndr_a(ll,mm,t1) * cos(mm*f1);
                    }
        } /* end if */
    } /* end for p */

    return;
}


static void add_perturbations_at_all_layers(struct All_variables *E)
{
    /* This function put a temperature anomaly for whole mantle with
       a sinosoidal amplitude in radial dependence. The horizontal pattern
       of the anomaly is given by wavenumber (in regional model) or
       by spherical harmonic (in global model). */

    int m, i, j, k, node;
    int p, ll, mm;
    int nox, noy, noz, gnoz;
    double r1, t1, f1, tlen, flen, rlen, con;

    nox = E->lmesh.nox;
    noy = E->lmesh.noy;
    noz = E->lmesh.noz;
    gnoz = E->mesh.noz;

    rlen = M_PI / (E->sphere.ro - E->sphere.ri);

    for (p=0; p<E->convection.number_of_perturbations; p++) {
        ll = E->convection.perturb_ll[p];
        mm = E->convection.perturb_mm[p];
        con = E->convection.perturb_mag[p];

        if (E->parallel.me_loc[1] == 0 && E->parallel.me_loc[2] == 0
            && E->sphere.capid[1] == 1 )
            fprintf(stderr,"Initial temperature perturbation:  mag=%g  l=%d  m=%d\n", con, ll, mm);

        if(E->sphere.caps == 1) {
            /* regional mode, add sinosoidal perturbation */

            tlen = M_PI / (E->control.theta_max - E->control.theta_min);
            flen = M_PI / (E->control.fi_max - E->control.fi_min);

            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(i=1; i<=noy; i++)
                    for(j=1; j<=nox;j ++)
                        for(k=1; k<=noz; k++) {
                            node = k + (j-1)*noz + (i-1)*nox*noz;
                            t1 = (E->sx[m][1][node] - E->control.theta_min) * tlen;
                            f1 = (E->sx[m][2][node] - E->control.fi_min) * flen;
                            r1 = E->sx[m][3][node];

                            E->T[m][node] += con * cos(ll*t1) * cos(mm*f1)
                                * sin((r1-E->sphere.ri) * rlen);
                        }
        }
        else {
            /* global mode, add spherical harmonics perturbation */

            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(i=1; i<=noy; i++)
                    for(j=1; j<=nox;j ++)
                        for(k=1; k<=noz; k++) {
                            node = k + (j-1)*noz + (i-1)*nox*noz;
                            t1 = E->sx[m][1][node];
                            f1 = E->sx[m][2][node];
                            r1 = E->sx[m][3][node];

                            E->T[m][node] += con * modified_plgndr_a(ll,mm,t1)
                                * (cos(mm*f1) + sin(mm*f1))
                                * sin((r1-E->sphere.ri) * rlen);
                        }
        } /* end if */
    } /* end for p */

    return;
}


static void add_spherical_anomaly(struct All_variables *E)
{
    int i, j ,k , m, node;
    int nox, noy, noz;

    double theta_center, fi_center, r_center,x_center[4],dx[4];
    double radius, amp, r1,rout,rin;
    const double e_4 = 1e-4;
    double distance;

    noy = E->lmesh.noy;
    nox = E->lmesh.nox;
    noz = E->lmesh.noz;

    rout = E->sphere.ro;
    rin = E->sphere.ri;


    theta_center = E->convection.blob_center[0];
    fi_center    = E->convection.blob_center[1];
    r_center     = E->convection.blob_center[2];
    radius       = E->convection.blob_radius;
    amp          = E->convection.blob_dT;

    if(E->parallel.me == 0)
      fprintf(stderr,"center=(%e %e %e) radius=%e dT=%e\n",
	      theta_center, fi_center, r_center, radius, amp);

    rtp2xyzd(r_center, theta_center, fi_center, (x_center+1));

    /* compute temperature field according to nodal coordinate */
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=1; i<=noy; i++)
            for(j=1; j<=nox;j ++)
                for(k=1; k<=noz; k++) {
                    node = k + (j-1)*noz + (i-1)*nox*noz;
		    dx[1] = E->x[m][1][node] - x_center[1];
		    dx[2] = E->x[m][2][node] - x_center[2];
		    dx[3] = E->x[m][3][node] - x_center[3];
                    distance = sqrt(dx[1]*dx[1] + dx[2]*dx[2] + dx[3]*dx[3]);

                    if (distance < radius){
		      E->T[m][node] += amp * exp(-1.0*distance/radius);

		      if(E->convection.blob_bc_persist){
			r1 = E->sx[m][3][node];
			if((fabs(r1 - rout) < e_4) || (fabs(r1 - rin) < e_4)){
			  /* at bottom or top of box, assign as TBC */
			  E->sphere.cap[m].TB[1][node]=E->T[m][node];
			  E->sphere.cap[m].TB[2][node]=E->T[m][node];
			  E->sphere.cap[m].TB[3][node]=E->T[m][node];
			}
		      }
		    }
                }
    return;
}


static void construct_tic_from_input(struct All_variables *E)
{
    double mantle_temperature;

    switch (E->convection.tic_method){
    case 0:
        /* a linear temperature profile + perturbations at some layers */
        linear_temperature_profile(E);
        add_perturbations_at_layers(E);
        break;

    case 1:
        /* T=1 for whole mantle +  cold lithosphere TBL */
        mantle_temperature = 1;
        constant_temperature_profile(E, mantle_temperature);
        add_top_tbl(E, E->convection.half_space_age, mantle_temperature);
        break;

    case 2:
        /* T='mantle_temp' for whole mantle + cold lithosphere TBL
           + a spherical anomaly at lower center */
        mantle_temperature = E->control.mantle_temp;
        constant_temperature_profile(E, mantle_temperature);
        add_top_tbl(E, E->convection.half_space_age, mantle_temperature);
        add_spherical_anomaly(E);
        break;

    case 3:
        /* a conductive temperature profile + perturbations at all layers */
        conductive_temperature_profile(E);
        add_perturbations_at_all_layers(E);
        break;

    case 4:
        /* read initial temperature from grd files */
#ifdef USE_GGRD
        ggrd_temp_init_general(E,1);
#else
        fprintf(stderr,"tic_method 4 only works for USE_GGRD compiled code\n");
        parallel_process_termination();
#endif
        break;

    case 10:
        /* T='mantle_temp' for whole mantle + cold lithosphere TBL
           + perturbations at some layers */

        mantle_temperature = E->control.mantle_temp;
        constant_temperature_profile(E, mantle_temperature);
        add_top_tbl(E, E->convection.half_space_age, mantle_temperature);
        add_perturbations_at_all_layers(E);
        break;

    case 11:
        /* T='mantle_temp' for whole mantle + hot CMB TBL
           + perturbations at some layers */

        mantle_temperature = E->control.mantle_temp;
        constant_temperature_profile(E, mantle_temperature);
        add_bottom_tbl(E, E->convection.half_space_age, mantle_temperature);
        add_perturbations_at_all_layers(E);
        break;

    case 12:
        /* T='mantle_temp' for whole mantle + cold lithosphere TBL
           + hot CMB TBL + perturbations at some layers */

        mantle_temperature = E->control.mantle_temp;
        constant_temperature_profile(E, mantle_temperature);
        add_top_tbl(E, E->convection.half_space_age, mantle_temperature);
        add_bottom_tbl(E, E->convection.half_space_age, mantle_temperature);
        add_perturbations_at_all_layers(E);
        break;

    case 90:
        /* for benchmarking purpose */
        /* a constant temperature (0) + single perturbation at mid-layer
           as a delta function in r */

        if((E->parallel.nprocz % 2) == 0) {
            if(E->parallel.me==0)
                fprintf(stderr, "ERROR: tic_method=%d -- nprocz is even, cannot put perturbation on processor boundary!\n",
                        E->convection.tic_method);

            parallel_process_termination();
        }

        constant_temperature_profile(E, 0);

        {
            /* adjust the amplitude of perturbation, so that
             * its integral in r is 1 */
            int mid, k;

            E->convection.number_of_perturbations = 1;

            mid = (E->mesh.noz+1) / 2;
            E->convection.load_depth[0] = mid;

            k = mid - E->lmesh.nzs + 1; /* convert to local nz */
            E->convection.perturb_mag[0] = 0;
            if ( (k > 1) && (k < E->lmesh.noz) ) {
                /* layer k is inside this proc. */
                E->convection.perturb_mag[0] = 2 / (E->sx[1][3][k+1] - E->sx[1][3][k-1]);
            }

        }
        add_perturbations_at_layers(E);
        break;

    case 100:
        /* user-defined initial temperature goes here */
        fprintf(stderr,"Need user definition for initial temperture: 'tic_method=%d'\n",
                E->convection.tic_method);
        parallel_process_termination();
        break;

    default:
        /* unknown option */
        fprintf(stderr,"Invalid value: 'tic_method=%d'\n", E->convection.tic_method);
        parallel_process_termination();
        break;
    }

    temperatures_conform_bcs(E);

    /* debugging the code of expanding spherical harmonics */
    /* debug_sphere_expansion(E);*/
    return;
}

void add_T_impact(struct All_variables *E)
{
    int m, i;
    const double r2d = 180.0 / M_PI, d2r = M_PI / 180.0;
    double colat, lon, dist;
    double Db   = E->impact.Db*1e-3;		// diameter of impact crater [km]
    double Dc   = E->impact.Dc*1e-3;		// simple-to-complex diameter [km]
    double vi   = E->impact.vi*1e-3;		// impact velocity [km/s]
    double g    = E->data.grav_acc;		// gravitational acceleration of[m/s^2]
    double uc   = vi / 2;	// particle velocity [km/s]
    double rho0 = E->data.density;		// mantle density [kg/m^3]
    double C    = E->impact.C*1e-3;		// mantle P-wave speed [m/s]
    double S    = E->impact.S;		// constant [1]
    double Cp   = E->data.Cp;		// specific heat [J/(kg*K)]
    double R    = E->data.radius*1e-3;		// planet radius [km]
    double Dtr;			// diameter of transient basin [m]
    double Dp;			// projectile size [m]
    double rc;			// radius of the isobaric core [m]
    double im_colat, im_lon;	// impact location [degree]
    double d, d1, d2, Dt, Pdel, beta, f, DT, P0, Ps, tau, radi;

    double dist_sphere();

    im_colat = 0;
    im_lon   = 0;

    Dtr = pow(Db / 1.02 * pow(Dc, 0.086), 1 / 1.086);
    Dp  = 0.69 * pow(Dtr * 1e3, 1.28) * pow(vi * 1e3, -0.56) * pow(g, 0.28) * 1e-3;

    rc  = 0.225 * Dp * pow(vi, 0.211);
    d   = 0.7 * 0.5 * Dp;
    tau = 0.5 * Dp / vi;

    if (E->parallel.me == 0)
	fprintf(stderr, "Dtr, Dp, rc, d, tau = %f %f %f %f %f\n", Dtr, Dp, rc, d, tau);
    for (m = 1; m <= E->sphere.caps_per_proc; m++)
    {
	for (i = 1;i <= E->lmesh.nno; i++)
	{
	    //E->T[m][i] = 0.0;

            colat = E->sx[m][1][i] * r2d;
	    lon   = E->sx[m][2][i] * r2d;
	    radi  = E->sx[m][3][i] * R;
	    dist  = dist_sphere(colat, lon, im_colat, im_lon) * d2r;	// radian

	    d1 = sqrt(pow(R - d, 2) + pow(radi, 2) - 2 * (R - d) * radi * cos(dist));
	    d2 = sqrt(pow(R + d, 2) + pow(radi, 2) - 2 * (R + d) * radi * cos(dist));

	    P0 = rho0 * g * fabs(R - radi);
	    Ps = rho0 * (C + S * uc) * uc * 1e3;
	    if (d1 > rc)
		Ps *= pow(rc / d1, -1.84 + 2.61 * log10(vi));

	    Dt = (d2 - d1) / C;
	    if (Dt < tau)
		Ps *= (Dt / tau);

	    Pdel = (Ps - P0) * 1e3;	// [Pa]
	    beta = pow(C * 1e3, 2) * rho0 / (2 * S);
	    f    = -Pdel / beta / (1 - sqrt(2 * Pdel / beta + 1));
	    DT   = Pdel * (1 - 1 / f) / (2 * rho0 * S * Cp)
		   - pow(C * 1e3 / S, 2) * (f - log(f) - 1) / Cp;

	    if (dist * r2d > 58 && dist * r2d < 62)
		    fprintf(stderr, "Ps, DT, radi = %f %f %f\n", Ps, DT, radi);

	    E->T[m][i] += DT / E->data.ref_temperature;
	}
    }
}

static void all_minerals(struct All_variables *E){
	//fprintf(stderr,"all_minerals begin\n"); //debug
	int num, i;
	E->mineral.num = 6; /*1 for harzburgite mantle , 2 for ilmenite, 3 for clinopyroxene,4 for anorthsite crust*/
	num = E->mineral.num;
/*initialize*/
	E->mineral.Density = (double*)malloc((num+1)*sizeof(double));
	for(i = 0; i <= E->mineral.num; i++)
		E->mineral.Density[i] = 0.0;
	E->mineral.buoyancy = (double*)malloc((num+1)*sizeof(double));
	E->mineral.ActE = (double*)malloc((num+1)*sizeof(double));
	E->mineral.actE = (double*)malloc((num+1)*sizeof(double));
	E->mineral.RefT = (double*)malloc((num+1)*sizeof(double));
	E->mineral.refT = (double*)malloc((num+1)*sizeof(double));
	E->mineral.offsetT = (double*)malloc((num+1)*sizeof(double));
	E->mineral.Eta0 = (double*)malloc((num+1)*sizeof(double));
	E->mineral.eta0 = (double*)malloc((num+1)*sizeof(double));
/*assign value*/
	/*mineral density kg/m^3*/
	if(E->lunar.model_type==1){
		 E->mineral.Density[1] = 3240.0; //for mg-suite model
		 E->mineral.Density[4] = 2800.0; //crust
		 E->mineral.Density[5] = 3210.0; //crust
	}
	else if(E->lunar.model_type==2){
		E->mineral.Density[1] = 3400.0;
		E->mineral.Density[4] = 2800.0; //crust
		E->mineral.Density[5] = 3400.0; //crust
		E->mineral.Density[6] = 3300.0; //last liquid of LMO, from Dygert et al 2017
	}

	else{
		E->mineral.Density[1] = 3400.0; //mantle reference
	}
	E->mineral.Density[2] = 4790.0; //ilmenite
	E->mineral.Density[3] = 3560.0; //cpx



	/*mineral activition energy J/mol*/
	E->mineral.ActE[1] = 500e3;
	E->mineral.ActE[2] = 281e3;
	E->mineral.ActE[3] = 500e3;
	E->mineral.preactE = 0.3;
	/*ref temperature*/
	E->mineral.RefT[1] = 1300 + 273.15;
	E->mineral.RefT[2] = 1300 + 273.15;
	E->mineral.RefT[3] = 1300 + 273.15;
	/*ref viscostiy*/
	E->mineral.Eta0[1] = E->data.ref_viscosity;
	E->mineral.Eta0[2] = 3e16;
	E->mineral.Eta0[3] = E->data.ref_viscosity;
//	fprintf(stderr,"all_minerals end\n"); //debug
}

/*get fraction of mineral for each chemical layer*/
static void all_chemicals(struct All_variables *E){
	//fprintf(stderr,"all_chemicals begin\n"); //debug
	int i, j;
	double thick;
	const int nflavors = E->trace.nflavors;
	E->chemical.ic_flavor = 2;
	if(E->lunar.model_type>0)
		E->trace.ic_flavor = E->chemical.ic_flavor;
	else
		E->trace.ic_flavor = 1;
	E->chemical.frac = (double**)malloc(nflavors*sizeof(double*));
	E->chemical.interface = (double**)malloc(nflavors*sizeof(double*));
	E->chemical.buoyancy = (double*)malloc(nflavors*sizeof(double));
	for(i=0;i<nflavors;i++){
		E->chemical.frac[i] = (double*)malloc((E->mineral.num + 1)*sizeof(double));
		E->chemical.interface[i] = (double*)malloc(3*sizeof(double));
	}
	for(i=0;i<nflavors;i++){
		for(j=0;j<=E->mineral.num;j++){
			E->chemical.frac[i][j] = 0.0;
		}
	}
	/*0 for harzburgite mantle*/
	i = (E->lunar.model_type==0)?1:2;
	thick =  1.0-E->viscosity.zlith-E->trace.z_interface[i-1];
	E->chemical.frac[0][1] = 1.0;
	if(E->lunar.model_type>0){
		E->chemical.frac[1][5] = 1.0;
		E->chemical.interface[1][1] = E->trace.z_interface[0];
		E->chemical.interface[1][2] = E->trace.z_interface[1]-0.3*thick;
	}

	i = (E->lunar.model_type==0)?1:2;
	/*1 for ilmenite layer below lithosphere*/
	E->chemical.frac[i][2] = ibc_layer(E,E->trace.z_interface[i-1]);
	E->chemical.frac[i][3] = E->chemical.frac[i][2]*(1-0.11)/0.11;  //0.11 is the volume fraction of ilmenite in the 20km ibc
	if(E->lunar.model_type==0){
		E->chemical.frac[i][1] = 1.0-E->chemical.frac[i][2]-E->chemical.frac[i][3];
	}
	else{
		//if take the mg-suite model, give this value to middle mantle fraction
		E->chemical.frac[i][5] = 1.0-E->chemical.frac[i][2]-E->chemical.frac[i][3];
	}
	E->chemical.interface[i][1] = E->trace.z_interface[i-1]-0.3*thick;
	E->chemical.interface[i][2] = 1.0-E->viscosity.zlith+0.1*thick;
	/*crust*/
	if(E->lunar.model_type>0){ //only in mg-suite model
		E->chemical.frac[3][4] = 1.0;
		E->chemical.interface[3][1] = 1.0-E->viscosity.zlith+0.1*thick;
		E->chemical.interface[3][2] = 1.0;
	}
	// ibc liquid
	if(E->lunar.model_type == 2){
		E->chemical.frac[nflavors - 1][6] = 1.0;
	}
	// derive chemical buoyancy
	for(i=0;i<E->trace.nflavors;i++){
		/*derive chemical buoyancy*/
		E->chemical.buoyancy[i] = 0.0;
		for(j=1;j<=E->mineral.num;j++){
			E->chemical.buoyancy[i] += E->chemical.frac[i][j]*E->mineral.buoyancy[j];
		}
	}
}

/*initialize model specified parameters*/
void specified_model_initialize(struct All_variables *E)
{
	//fprintf(stderr,"2\n"); //debug
	int Is0 = (E->parallel.me==0);
	int up;
	/*
	if (Is0){
		fprintf(stderr,"specified_model_initialize begin\n");
	}
	*/
	char filename[1000];
	int i;
	/*read in parameters from file*/
	lunar_read_parametes(E);
	/*give value for dimensional values*/
	/*if(Is0){
		fprintf(stderr,"begin all_minerals\n");
	}
	*/
	all_minerals(E);
	/*if(Is0){
		fprintf(stderr,"begin none_dimensional_parameters\n");
	}*/
	/*get nondimensional values*/
	none_dimensional_parameters(E);
	/*derive chemical values*/
	/*if(Is0){
		fprintf(stderr,"begin all_chemicals\n");
	}*/
	all_chemicals(E);
	/*ralate with global defined parameters*/
	if (E->lunar.bfit){ //buoyancy_from_ic_thickness
		if(E->control.t_ic_sol<0)
			up = E->trace.nflavors-1;
		else{
			up = E->trace.nflavors-2;
			E->composition.buoyancy_ratio[up] = 0.0;
		}
		for(i=0;i<up;i++){
			E->composition.buoyancy_ratio[i] = E->chemical.buoyancy[i+1];
		}
	}
	/*write specified variables to output*/
	if(E->parallel.me == 0){
		sprintf(filename,"%s/specified_parameters",E->control.data_dir);
		write_specified_model_parameters(E,filename);
	}
	if(Is0){
		fprintf(stderr,"specified_model_initialize end\n");//debug
	}
}

static void lunar_read_parametes(struct All_variables *E){
	int input_int();
	int input_double();
	int m = E->parallel.me;
	input_int("model_type",&(E->lunar.model_type),"0",m);
	input_int("latent_method",&(E->lunar.latent_method),"0",m);
	input_double("latent_heat",&(E->lunar.latent_heat),"6e5",m);
	input_int("buoyancy_from_ic_thickness",&(E->lunar.bfit),"0",m);
	input_int("smooth_lower_half",&(E->lunar.smooth_lower_half),"0",m);
	input_double("lower_interface",&(E->lunar.lower_interface),"0.5977",m);
    input_int("smooth_upper_layer",&(E->lunar.smooth_upper_layer),"0",m);       //zwb 20200715
	input_double("upper_interface",&(E->lunar.upper_interface),"0.9529",m);         //zwb 20200715
	input_int("smooth_downwelling",&(E->lunar.smooth_downwelling),"0",m);       //zwb 20200730
    input_float_vector("smooth_cdepv",E->trace.nflavors, (E->lunar.smooth_cdepv),m);        //zwb 20200730
	if(m==0){
		if(E->lunar.model_type == 1){
			fprintf(stderr,"start a mg-suite model\n");
		}
		else if(E->lunar.model_type == 2){
			fprintf(stderr,"start a syn-crys model\n");
		}
		else{
			fprintf(stderr,"start a ibc-chemical or impact model\n");
		}
	}
	/*
	if(m == 0){
		fprintf(stderr,"latent_method = %d\n",E->lunar.latent_method);
		fprintf(stderr,"latent_heat = %.4e\n",E->lunar.latent_heat);
	}*/ //debug

}

void process_melt_factor(struct All_variables *E, double **melt_fac, int e, int method){
	//fprintf(stderr,"process_melt_factor begin\n"); //debug
	int noz, lnz ,nz, dims, ends, vpts, i, j, node, node1, lnz1, nz1;
	double fac, gTs;
	double T[9], Tdiff[9], Ts[9], Tl[9], ur[9];
	double tiny = 1e-10;
	int quiry = -1; //debug
	//assign values
	dims = E->mesh.nsd;
	ends = enodes[dims];
	vpts = vpoints[dims];
	noz = E->lmesh.noz;
	//fprintf(stderr,"L = %.4e, cp = %.4e, refT = %.4e\n",E->lunar.latent_heat,E->data.Cp,E->data.ref_temperature);
	fac = E->lunar.latent_heat/(E->data.Cp*E->data.ref_temperature);
	//derive gTs(gradient of solidus*/
	node = E->ien[1][e].node[1];
	node1 = E->ien[1][e].node[5];
	lnz = node%noz;
	nz = ((lnz == 0)?noz:lnz) + E->parallel.me_loc[3] * (noz - 1);
	lnz1 =node1%noz;
	nz1 = ((lnz1 == 0)?noz:lnz1) + E->parallel.me_loc[3] * (noz - 1);
	gTs = (E->gsol[1][nz] + E->gsol[1][nz1])/2.0;
	//fprintf(stderr,"fac = %.4e, gTs = %.4e\n",fac,gTs); //debug
	/*get value at gauss points*/
	for(i=1;i<=vpts;i++){
		T[i] = 0.0;
		Tdiff[i] = 0.0;
		Ts[i] = 0.0;
		Tl[i] = 0.0;
		ur[i] = 0.0;
	}
	for(j=1;j<=ends;j++){
		node = E->ien[1][e].node[j];
		lnz = node%noz;
		nz = ((lnz == 0)?noz:lnz) + E->parallel.me_loc[3] * (noz - 1);
		for(i=1;i<=vpts;i++){
			T[i] += E->N.vpt[GNVINDEX(j,i)]*E->T[1][node];
			Tdiff[i] += E->N.vpt[GNVINDEX(j,i)]*E->Tdiff[1][node];
			ur[i] += E->N.vpt[GNVINDEX(j,i)]*E->sphere.cap[1].V[3][node];
			Ts[i] += E->N.vpt[GNVINDEX(j,i)]*E->sol[1][nz];
			Tl[i] += E->N.vpt[GNVINDEX(j,i)]*E->liq[1][nz];
		}
	}
	//fprintf(stderr,"T = %.4e, ur = %.4e, Ts = %.4e, Tl = %.4e\n",T[quiry],ur[quiry],Ts[quiry],Tl[quiry]);
	for(i=1;i<=vpts;i++){
		if(T[i]>Ts[i]+tiny){
			if(method == 1){
				melt_fac[i][1] = 1.0;
				melt_fac[i][2] = -1.0*fac/(Tl[i]-Ts[i])*Tdiff[i];
			}
			else if(method == 2){
				melt_fac[i][1] = 1.0 + fac/(Tl[i]-Ts[i]);
				melt_fac[i][2] = fac*gTs*ur[i]/(Tl[i]-Ts[i]);
			}
		}
		else{
			melt_fac[i][1] = 1.0;
			melt_fac[i][2] = 0.0;
		}
	}
	//fprintf(stderr,"meltfac1 = %.4e, meltfac2 = %.4e\n",melt_fac[quiry][1],melt_fac[quiry][2]);
	//fprintf(stderr,"process_melt_factor end\n"); //debug
	return;
}
/*lhy 20180719
 * get nondimensional density at every node and calculate nodimentional total mass*/
void rho_at_nodes(struct All_variables *E){
	double Tref,T0;
	double coeff1, coeff2, fac1, fac2;
	int m,i,j,nz;
	int Is0 = (E->parallel.me==0);
	int query = 12;

	if(Is0){//debug
		fprintf(stderr,"rho_at_nodes starts\n");
	}

	for(m=1;m<=E->sphere.caps_per_proc;m++){
		E->rho_node[m] = (double*)malloc((E->lmesh.nno+1)*sizeof(double));
	}
	/*assign coefficients*/
	Tref = E->data.ref_temperature;
	T0 = 273.15;
	coeff1 = E->data.therm_exp*Tref;
	coeff2 = (T0-E->data.Tsurf)/Tref;
	if(Is0){//debug
		fprintf(stderr,"coeff1 = %.4e, coeff2 = %.4e\n",coeff1,coeff2);
	}
	/*calculate density*/
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=1;i<=E->lmesh.nno;i++) {
		nz = ((i-1) % E->lmesh.noz) + 1;
		/*fac1 represents density variation by chemical*/
		fac1 = 1.0;
      	for(j=0;j<E->composition.ncomp;j++) {
			fac1 += coeff1*E->composition.buoyancy_ratio[j]*E->composition.comp_node[m][j][i];
		}
		/*fac2 represents density variation by temperature*/
		fac2 = (1-coeff1*(E->T[m][i]-coeff2));
		E->rho_node[m][i] = E->refstate.rho[nz]*fac1*fac2;
		if(Is0&&(query==i)){//debug
			fprintf(stderr,"node = %d, fac1 = %.4e, fac2 = %.4e, rho = %.4e\n",
					i,fac1,fac2, E->rho_node[m][i]);
		}
	  }

	E->total_mass = return_bulk_value_d(E,E->rho_node,0);

	if(Is0){//debug
		fprintf(stderr,"total_mass = %.4e\n",E->total_mass);
		fprintf(stderr,"rho_at_nodes ends\n");
	}
	return;
}
/*lhy 20180719
 * calculate mass center*/
void mass_center(struct All_variables *E){
	int m,i;
	int Is0 = (E->parallel.me==0);
	int query = 12;
  	double **rho_x, **rho_y, **rho_z;

	if(Is0){//debug
		fprintf(stderr,"mass_center starts\n");
	}
/*initialize*/
	rho_x = (double**)malloc((E->sphere.caps_per_proc+1)*sizeof(double*));
	rho_y = (double**)malloc((E->sphere.caps_per_proc+1)*sizeof(double*));
	rho_z = (double**)malloc((E->sphere.caps_per_proc+1)*sizeof(double*));
	for(m=0;m<=E->sphere.caps_per_proc;m++){
		rho_x[m] = (double*)malloc((E->lmesh.nno+1)*sizeof(double));
		rho_y[m] = (double*)malloc((E->lmesh.nno+1)*sizeof(double));
		rho_z[m] = (double*)malloc((E->lmesh.nno+1)*sizeof(double));
	}

	for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=1;i<=E->lmesh.nno;i++) {
		  rho_x[m][i] = E->rho_node[m][i]*E->x[m][1][i];
		  rho_y[m][i] = E->rho_node[m][i]*E->x[m][2][i];
		  rho_z[m][i] = E->rho_node[m][i]*E->x[m][3][i];
		  /*if(Is0&&(query==i)){//debug
			  fprintf(stderr,"node = %d, rho_x = %.4e, rho_y = %.4e, rho_z = %.4e\n",i,rho_x[m][i],rho_y[m][i],rho_z[m][i]);
		  }*/
	  }

	E->mc[1] = return_bulk_value_d(E,rho_x,0)/E->total_mass;
	E->mc[2] = return_bulk_value_d(E,rho_y,0)/E->total_mass;
	E->mc[3] = return_bulk_value_d(E,rho_z,0)/E->total_mass;

	for(m=0;m<=E->sphere.caps_per_proc;m++){
		free(rho_x[m]);
		free(rho_y[m]);
		free(rho_z[m]);
	}
	free(rho_x);
	free(rho_y);
	free(rho_z);

	/*if(Is0){//debug
		fprintf(stderr,"mass_center ends\n");
	}*/
	return;
}

/*lhy 20180720
 * correct topography with mass center*/
void correct_for_mc(struct All_variables *E,float **topo,double *mc){
	int i,s;
	double d,x,y,z;
	double scaling;
	int Is0 = (E->parallel.me==1);
/*
	if(Is0){//debug
		fprintf(stderr,"correct_for_mc starts\n");
		fprintf(stderr,"mc[1] = %.4e, mc[2] = %.4e, mc[3] = %.4e\n",mc[1],mc[2],mc[3]);
	}
*/
    /*get scaling value for mc*/
	scaling = E->control.Atemp/(E->data.therm_exp*E->data.ref_temperature);
	for(i=1;i<=E->lmesh.nsf;i++) {
        s = i*E->lmesh.noz;
		x = E->x[1][1][s];
		y = E->x[1][2][s];
		z = E->x[1][3][s];
		d = sqrt(pow(x-mc[1],2.0)+pow(y-mc[2],2.0)+pow(z-mc[3],2.0))-E->sphere.ro;
		E->slice.tpgmc[1][i] = topo[1][i]+d*scaling;
	}
	/*if(Is0){//debug
		fprintf(stderr,"correct_for_mc ends\n");
	}*/
}
/*lhy 20170720
 * calculate sph expansion for topography*/
void topo_sphere_expansion(double *topo_sph_sum, double **topo_sph,struct All_variables *E,float **topo){
	int i,mm,ll,pp;
	int ll_max = E->output.llmax;
	int Is0 = (E->parallel.me==1);

	/*if(Is0){//debug
		fprintf(stderr,"topo_sphere_expansion starts\n");
	}*/
	/* do sphere expansion for topography*/
	sphere_expansion(E, topo, E->sphere.harm_tpgt[0], E->sphere.harm_tpgt[1]);
	/*add topography to get multitude for individual degrees*/
	for(i=1;i<=E->lmesh.nsf;i++)   {
		topo_sph_sum[i] = 0.0;
		for(ll = 0; ll<=ll_max ; ll++){
			topo_sph[ll][i] = 0.0;
			for(mm = 0; mm<=ll; mm++){
				pp = E->sphere.hindex[ll][mm];
				topo_sph[ll][i] += E->sphere.harm_tpgt[0][pp]*E->sphere.tablesplm[1][i][pp]*E->sphere.tablescosf[1][i][mm]+
				E->sphere.harm_tpgt[1][pp]*E->sphere.tablesplm[1][i][pp]*E->sphere.tablessinf[1][i][mm];
			}
			topo_sph_sum[i] += topo_sph[ll][i];
		}
	}

	/*if(Is0){//debug
		fprintf(stderr,"topo_sphere_expansion ends\n");
	}*/
	return;
}

/*get system viscosity using model defined method for mixed material viscosity*/
void visc_from_TC_by_mineral(struct All_variables *E,float **EEta){
	int is0 = (E->parallel.me==1);
	if(is0){
		fprintf(stderr,"visc_from_TC_by_mineral begins\n");
	}//debug*/
	int m,i,jj,kk,ends,vpts,nel,p,nflavors,mnum,mtr;
	double temp;
	double *TT, *cc_loc, *etam;
	double **CC;
	int query=56; //debug
	/*initialize*/
	vpts = vpoints[E->mesh.nsd];
	ends = enodes[E->mesh.nsd];
	nel = E->lmesh.nel;
	nflavors = E->trace.nflavors;
	mnum = E->mineral.num;
	CC = (double**)malloc(nflavors*sizeof(double*));
	for(p=0;p<E->trace.nflavors;p++){
		CC[p] = (double*)malloc((ends+1)*sizeof(double));
	}
	cc_loc = (double*)malloc(nflavors*sizeof(double));
	TT = (double*)malloc((ends+1)*sizeof(double));
	etam = (double*)malloc((mnum+1)*sizeof(double));
	/*system viscosity*/
	for(m=1;m<=E->sphere.caps_per_proc;m++){
		for(i=1;i<=nel;i++){
			int isO = (i==query); //debug
			/*get temperature and composition on nodes*/
			if(isO&&is0){
				fprintf(stderr, "element %d\n",nel);
			} //debug
            for(kk=1;kk<=ends;kk++) {
				if(isO&&is0){
					fprintf(stderr,"node %d\n",kk);
				}//debug
		  		TT[kk] = E->T[m][E->ien[m][i].node[kk]];
				CC[0][kk] = 1.0;
				for(p=0; p<E->composition.ncomp; p++) {
                	CC[p+1][kk] = E->composition.comp_node[m][p][E->ien[m][i].node[kk]];
                	if(CC[p+1][kk] < 0)CC[p+1][kk]=0.0;
                	if(CC[p+1][kk] > 1)CC[p+1][kk]=1.0;
					if(E->lunar.smooth_lower_half){
						if(E->sx[1][3][E->ien[m][i].node[kk]]<E->lunar.lower_interface){
							CC[p+1][kk] = 0.0;
						}
					}
					CC[0][kk] -= CC[p+1][kk];
				}
				if(isO&&is0){
					fprintf(stderr,"temperature = %.4e\n",TT[kk]);
					for(p=0; p<nflavors; p++){
						fprintf(stderr,"comp%d = %.4e\n",p,CC[p][kk]);
					}
				}//debug
            }
            for(jj=1;jj<=vpts;jj++) {
				if(isO&&is0){
					fprintf(stderr,"vpt %d\n",jj);
				}//debug
				/*get temperature and composition on gaussin points*/
            	temp=0.0;
            	for(p=0; p<nflavors; p++) {
                	cc_loc[p] = 0.0;
				}
                for(kk=1;kk<=ends;kk++)   {
                	temp += TT[kk] * E->N.vpt[GNVINDEX(kk,jj)];
            		for(p=0; p<nflavors; p++) {
                    	cc_loc[p] += CC[p][kk] * E->N.vpt[GNVINDEX(kk, jj)];
                	}
                }
				if(isO&&is0){
					fprintf(stderr,"temperature = %.4e\n",temp);
					for(p=0; p<nflavors; p++){
						fprintf(stderr,"comp%d = %.4e\n",p,cc_loc[p]);
					}
				}//debug
				/*get viscosity for all minerals*/
				for(mtr=1;mtr<=E->mineral.num;mtr++){
					etam[mtr] = T_dependent_viscosity(E,temp,E->mineral.eta0[mtr],E->mineral.actE[mtr],E->mineral.refT[mtr],E->mineral.offsetT[mtr]);
					if(is0&&isO){
						fprintf(stderr,"material = %d, viscosity = %.4e\n",mtr,etam[mtr]);
					}//debug
				}
				/*apply model defined rule for mixed material viscosity*/
				 EEta[m][(i-1)*vpts+jj] = visc_from_chemical(E, cc_loc, etam);
				if(is0&&isO){
					fprintf(stderr,"gauss point viscosity = %.4e\n",EEta[m][(i-1)*vpts+jj]);
				}//debug
			}
		}
	}
	/*free space*/
	free(etam);
	free(TT);
	free(cc_loc);
	for(p=0;p<E->trace.nflavors;p++){
		free(CC[p]);
	}
	free(CC);
	if(is0){
		fprintf(stderr,"visc_from_TC_by_mineral ends\n");
	}//debug*/
}

/*derive chemical from chemical fraction*/
static double visc_from_chemical(struct All_variables *E, double *chemical, double *cvisc){
	double *mineral;
	double visc,x;
	int i,j;
	/*initialize*/
	mineral = (double*)malloc((E->mineral.num+1)*sizeof(double));
	/*derive mineral fraction*/
	for(i=0;i<E->trace.nflavors;i++){
		for(j=1;j<=E->mineral.num;j++){
			mineral[j] += chemical[i]*E->chemical.frac[i][j];
		}
	}
	/*derive viscosity using ilmenite viscosity and isostress model*/
	x = mineral[2]; /*fraction of ilmentite*/
	//fprintf(stderr,"ilmenite fraction: %.4e\n",x); //debug
	if(x>0.04){
		visc = iso_stress(cvisc[1],cvisc[2], 1.0-x);
	}
	else{
		visc = cvisc[1];
	}
	return visc;
}

/*get nondimensional parameters for minerals*/
static double none_dimensional_parameters(struct All_variables *E){
//	fprintf(stderr,"none_dimensional_parameters begin\n"); //debug
	int num,i;
	double grav,radius,refrho,rayleigh,kappa,alpha,refvisc,preE,Rgas,reftemp,surftemp;
	num = E->mineral.num;

	/*assign dimensional variables*/
	grav = E->data.grav_acc;
	radius = E->data.radius;
	refrho = E->data.density;
	rayleigh = E->control.Atemp;
	alpha = E->data.therm_exp;
	kappa = E->data.therm_diff;
	refvisc = E->data.ref_viscosity;
	preE = E->mineral.preactE;
	reftemp = E->data.ref_temperature;
	surftemp = E->data.Tsurf;
	Rgas = 8.3144;
/*	fprintf(stderr,"grav = %.4e, radius = %.4e, refrho = %.4e, rayleigh = %.4e,\
	alpha = %.4e, kappa = %.4e, refvisc = %.4e, preE = %.4e, reftemp = %.4e, surftemp = %.4e\n",\
	grav, radius, refrho, rayleigh, alpha, kappa, refvisc, preE, reftemp, surftemp); */
	/*derive nondimensional variables*/
	for(i=1;i<=num;i++){
		E->mineral.buoyancy[i] = (E->mineral.Density[i]-refrho)/(alpha*refrho*reftemp);
		E->mineral.eta0[i] = E->mineral.Eta0[i]/refvisc;
		E->mineral.actE[i] = preE*E->mineral.ActE[i]/(reftemp*Rgas);
		E->mineral.refT[i] = (E->mineral.RefT[i]-surftemp)/reftemp;
		E->mineral.offsetT[i] = surftemp/reftemp;
	}
//	fprintf(stderr,"none_dimensional_parameters end\n"); //debug
}

/*get ilmenite fraction in the ibc_layer layer thickness is derived via E->trace_zinterface*/
static double ibc_layer(struct All_variables *E, double rad){
/*rad is lower interface of ibc layer*/
	double h0,rad0,rad1,rad2,frac0,frac;
	h0 = 20e3/E->data.radius; //thickness of ibc in the LMO crystallization
	frac0 = 0.11; //volume fraction 0.11 corresponding to weight fraction 0.15
	rad0 = 1-E->viscosity.zlith-h0;
	rad2 = 1-E->viscosity.zlith;
	frac = frac0*(pow(rad2,3.0)-pow(rad0,3.0))/(pow(rad2,3.0)-pow(rad,3.0));
	return frac;
}

/*write model specified parameters into a output file*/
static void write_specified_model_parameters(struct All_variables *E, char* filename){
	FILE *fp;
	double *frac,*chemical,*etam,*eta,*interface;
	int i,j,mnum,cnum,index;
	double T;
	fprintf(stderr,"writing to file %s\n",filename);
	/*initialize*/
	mnum = E->mineral.num;
	if(E->control.t_ic_sol>0.0){
		cnum = E->trace.nflavors-1;
	}
	else{
		cnum = E->trace.nflavors;
	}
	T = (E->mineral.RefT[1] - E->data.Tsurf)/E->data.ref_temperature;
	frac = (double*)malloc(mnum*cnum*sizeof(double));
	interface = (double*)malloc(2*(cnum-1)*sizeof(double));
	chemical = (double*)malloc(cnum*sizeof(double));
	eta = (double*)malloc(cnum*sizeof(double));
	etam = (double*)malloc((mnum+1)*sizeof(double));
	/*output*/
	fp = fopen(filename, "w");
	if(fp==NULL){
		fprintf(E->fp,"(write_specified_model_parameters.c #1) cannot open f%s for writing\n",filename);
		parallel_process_termination();
	}
	/*print mineral parameters*/
	fprintf(fp,"#mineral\n");
	write_double_array_into_file(fp,"mineral_buoyancy",E->mineral.buoyancy,mnum,1);
	write_double_array_into_file(fp,"mineral_eta0",E->mineral.eta0,mnum,1);
	write_double_array_into_file(fp,"mineral_actE",E->mineral.actE,mnum,1);
	write_double_array_into_file(fp,"mineral_refT",E->mineral.refT,mnum,1);
	write_double_array_into_file(fp,"mineral_offsetT",E->mineral.offsetT,mnum,1);
	/*printf chemical parameters and interfaces*/
	for(i=0;i<cnum;i++){
		for(j=0;j<mnum;j++){
			index = i*mnum+j;
			frac[index] = E->chemical.frac[i][j+1];
		}
		if(i>=1){
			for(j=0;j<2;j++){
				index = 2*(i-1)+j;
				interface[index] = E->chemical.interface[i][j+1];
			}
		}
	}
	fprintf(fp,"#chemical\n");
	write_double_array_into_file(fp,"chemical_fraction",frac,cnum*mnum,0);
	write_double_array_into_file(fp,"chemical_interfaces",interface,(cnum-1)*2,0);
	if(E->control.t_ic_sol>0.0){
		write_double_array_into_file(fp,"chemical_buoyancy",E->chemical.buoyancy,cnum+1,0);
	}
	else{
		write_double_array_into_file(fp,"chemical_buoyancy",E->chemical.buoyancy,cnum,0);
	}
	/*derive chemical viscosity*/
	for(i=1;i<=E->mineral.num;i++){
		etam[i] = T_dependent_viscosity(E, T, E->mineral.eta0[i], E->mineral.actE[i], E->mineral.refT[i], E->mineral.offsetT[i]);
	}
	for(i=0;i<cnum;i++){
		for(j=0;j<cnum;j++){
			chemical[j] = 0.0;
		}
		chemical[i] = 1.0;
		eta[i] = visc_from_chemical(E, chemical, etam);
	}
	write_double_array_into_file(fp,"chemical_eta0",eta,cnum,0);
	/*free space*/
	free(interface);
	free(eta);
	free(etam);
	free(chemical);
	free(frac);
	fclose(fp);
}

/*write a double array into a opened file*/
static void write_double_array_into_file(FILE *fp, char *name, double *array, int num, int type){
	//fprintf(stderr,"array name: %s\n",name); //debug
	int start,end,i;
	/*determine type, 0 for array starts from index 0, 1 for array starts from index 1*/
	if (type == 0){
		start = 0;
		end = num-1;
	}
	else{
		start =1;
		end = num;
	}
	/*output array into file fp*/
	fprintf(fp,"%s=",name);
	for(i=start;i<=end;i++){
		if(i>start){
			fprintf(fp,",");
		}
		fprintf(fp,"%.4e",array[i]);
	}
	fprintf(fp,"\n");
}
static double T_dependent_viscosity(struct All_variables *E, double T, double eta0, double actE, double refT, double offsetT){
	double eta;
	eta = eta0 * exp(actE/(T+offsetT)-actE/(refT+offsetT));
	return eta;
}
/*return viscosity of mixed material with isostress model*/
static double iso_stress(double visc1,double visc2, double x){
	double visc,strain;
	strain = x/visc1 + (1.0-x)/visc2;
	visc = 1.0/strain;
	return visc;
}

static void read_temp_from_radial(struct All_variables *E)
{
const int Is0 = (E->parallel.me == 0);
const double tiny = 1e-7;
char input_s[200];

int keep_going;
int kk;
int num_input_nodes;
int input_element;
int j,mm;
int nx,ny,nz;
int node;

double rad_bottom;
double rad_top;
double rad;
double eta,delrad;
double shape1,shape2,temperature;

double *input_rad;
double *input_t;
double *citcom_t;

int size_of_input=1000;
int noz=E->lmesh.noz;
int nox=E->lmesh.nox;
int noy=E->lmesh.noy;

FILE *fp_read;

    input_rad=(double *)malloc(size_of_input*sizeof(double));
    input_t=(double *)malloc(size_of_input*sizeof(double));
    citcom_t=(double *)malloc((noz + 1)*sizeof(double));

/* read input data */
    if ((fp_read=fopen(E->control.background_profile_file,"r"))==NULL)
    {
       fprintf(stderr,"ERROR(read temp from radial)-no file\n");
       fflush(stderr);
       exit(10);
    }

    if (E->parallel.me==0) fprintf(stderr,"Initial Radial Temperature Field: %s\n",E->control.background_profile_file);

    keep_going=1;
    kk=0;

    while(keep_going)
    {
       if (fgets(input_s,200,fp_read)==NULL)
       {
          keep_going=0;
       }
       else
       {
          kk++;
          if (kk>(size_of_input-1))
          {
            fprintf(stderr,"ERROR(read_temp_from_radial) - file too big. Increase size_of_input\n");
            fflush(stderr);
            exit(10);
          }
      //sscanf(input_s,"%lf %lf",&input_t[kk],&input_rad[kk]);
          sscanf(input_s,"%lf %lf",&input_rad[kk],&input_t[kk]);//nanzhang
/* some control checks on input file */
          if ((kk>1)&&(input_rad[kk]<=input_rad[kk-1]))
          {
              fprintf(stderr,"ERROR(read_temp_from_radial)-rad does not increase? -check input file\n");
              fprintf(stderr,"rad[kk-1]: %f rad[kk]: %f kk: %d\n",input_rad[kk-1],input_rad[kk],kk);
              exit(10);
          }
          /*if ((input_rad[kk]<E->sphere.ri)||(input_rad[kk]>E->sphere.ro))
          {
              fprintf(stderr,"ERROR(read_temp_from_radial)-wierd rad? -check input file\n");
              fprintf(stderr,"input rad: %f\n",input_rad[kk]);
              fflush(stderr);
              exit(10);
          } */
       }
    }

    num_input_nodes=kk;

        if (num_input_nodes<2)
    {
       fprintf(stderr,"ERROR(read_temp_from_radial) - need at least 2 input points!\n");
       fflush(stderr);
       exit(10);
    }

    fclose(fp_read);
/* interpolate citcom nz nodes */
    for (j=1;j<=E->sphere.caps_per_proc;j++)
    {

    for (kk=1;kk<=noz;kk++)
    {
       rad=E->sx[j][3][kk];
       if (rad>E->sphere.ro) rad = E->sphere.ro; //nanzhang
/* find which input element */

       input_element=0;
       for (mm=1;mm<=num_input_nodes-1;mm++)
       {
          rad_bottom=input_rad[mm];
          rad_top=input_rad[mm+1];
          if ( (rad>rad_bottom-tiny) && (rad<rad_top+tiny))
          {
             input_element=mm;
             goto foundit;
          }
       }

foundit:;

/* find local coordinate,eta, in input element. */
/* here, local coordinate extends from 0 to 1. */
       delrad=rad_top-rad_bottom;
       eta=(rad-rad_bottom)/delrad;

       if ((eta<-1e-6)||(eta>(1.0+1e-4)))
       {
          fprintf(stderr,"ERROR(read_temp_radial) z from m %d %d %f %f %f %f %f\n",kk,mm,eta,rad_bottom,rad_top,rad,delrad);
          fflush(stderr);
          exit(10);
       }

/* find shape functions at local coordinate, eta */
       shape1=input_t[input_element]*(1.0-eta);
       shape2=input_t[input_element+1]*(eta);
       temperature=shape1+shape2;
       citcom_t[kk]=temperature;
    }

/* now fill citcom nodes */
    for (ny=1;ny<=noy;ny++)
    {
    for (nx=1;nx<=nox;nx++)
    {
    for (nz=1;nz<=noz;nz++)
    {
       node=nz+(nx-1)*noz+(ny-1)*(nox*noz);
       E->T[j][node]=citcom_t[nz];
    }
    }
    }
    } /* end j */
    free(citcom_t);
    free(input_rad);
    free(input_t);
    return;

}

static reset_IBC(struct All_variables *E){
	const int Is1 = (E->parallel.me==1);
	if(Is1)
		fprintf(stderr,"reset_IBC begin\n"); 	// debug
	const int ends = enodes[E->mesh.nsd];
	const int vpts = vpoints[E->mesh.nsd];
	int j, ic, kk, e, n, k, elz;
	float T, Tavg, weight;
	float T_IBC_reset = E->control.T_IBC_reset;
	compute_horiz_avg(E);
	ic = E->trace.ic_flavor;
	const int number_of_tracers = E->trace.ntracers[1];
	for (kk=1;kk<=number_of_tracers;kk++) {
		if(E->trace.extraq[1][0][kk] == ic){
			e = E->trace.ielement[1][kk];
			T = 0.0;
			for(k=1;k<=ends;k++){
				weight = 0.0;
				n = E->ien[1][e].node[k];
				for(j=1;j<=vpts;j++)
					weight += E->N.vpt[GNVINDEX(k,j)] * E->gDA[1][e].vpt[j];
				T += E->T[1][n]*weight;
			}
			T = T/E->eco[1][e].area;
			elz = (e-1)%E->lmesh.elz + 1;
			Tavg = 0.0;
			Tavg += E->Have.T[elz]*E->gDA[1][e].vpt[1];
			Tavg += E->Have.T[elz+1]*E->gDA[1][e].vpt[5];
			Tavg = Tavg/(E->gDA[1][e].vpt[1]+E->gDA[1][e].vpt[5]);
			if(T-Tavg > T_IBC_reset){
         		E->trace.extraq[1][0][kk] = ic-1;
			}
		}
	}
    if (E->trace.nflavors > 0)
        count_tracers_of_flavors(E);
	fill_composition(E);
    map_composition_to_nodes(E);
}

// surf area of each points
void output_surf_area(E)
     struct All_variables *E;
{
  char output_file[255];
  int i, j, k, nproc;
  int el, elz, elx, ely;
  const int nsd = E->mesh.nsd;
  const int Is0 = (E->parallel.me == 0);
  double x_e[4], rtf_e[4];
  FILE *fp;
  if(Is0)
  	fprintf(stderr, "output_surf_area\n");
  if (E->parallel.me % E->parallel.nprocz == E->parallel.nprocz - 1)
  {
	  sprintf(output_file, "%s.surf_area.%d", E->control.data_file,
     		  E->parallel.me);
  }
  else
  {
	  return;
  }
  fp = fopen(output_file, "w");
  elz = E->lmesh.elz;
  elx = E->lmesh.elx;
  ely = E->lmesh.ely;
  i = 0;
  for (k = 1; k <= ely; k++)
  for (j = 1; j <= elx; j++)
  {
      el = elz + (j-1)*elz + (k-1)*elx*elz;
	  i++;
	  element_coordinate(E, el, x_e, rtf_e);   // coordinate of element
	  fprintf(fp, "%d %.4e %.4e %.4e %.4e\n", i, rtf_e[1], rtf_e[2],
	  		  rtf_e[3], E->eco[1][el].area);
  }
  fclose(fp);
  return;
}

// another method
void output_surf_area1(E)
     struct All_variables *E;
{
  	const int Is0 = (E->parallel.me == 0);
  	char output_file[255];
	FILE *fp;
	int es, nint;
	double area;
  	if(Is0)
  		fprintf(stderr, "output_surf_area1\n");
  	if (E->parallel.me % E->parallel.nprocz == E->parallel.nprocz - 1)
  	{
	  	sprintf(output_file, "%s.surf_area.%d", E->control.data_file,
     		  	E->parallel.me);
  	}
  	else
  	{
	  	return;
  	}
  	fp = fopen(output_file, "w");
	for (es=1;es<=E->lmesh.snel;es++)   {
		area = 0.0;
		for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)
			area += E->surf_det[1][nint][es];
	  	fprintf(fp, "%d %.4e\n", es, area);
	}
	fclose(fp);
	return;
}


// FE functions by lhy

// Element coordinate

void element_coordinate(E, el, x_e, rtf_e)
	struct All_variables *E;
	int el;
	double *x_e, *rtf_e;
{
	const int nsd = E->mesh.nsd;
	const int e_nodes = enodes[nsd];
	int i, j, lnode;
	for(i = 0; i <= nsd; i++) // initial 0.0
		x_e[i] = 0;
	for(i = 1; i <= nsd; i++) // cartesion of element
	for(j = 1; j <= e_nodes; j++)
	{
		lnode = E->ien[1][el].node[j];
		x_e[i] += E->x[1][i][lnode];
	}
	for(i = 0; i <= nsd; i++)
		x_e[i] /= e_nodes;
	cart_to_sphere(E, x_e[1], x_e[2], x_e[3],
				   &rtf_e[1], &rtf_e[2], &rtf_e[3]); // spherical of element
	return;
}
