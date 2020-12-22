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
/*

  Tracer_setup.c

      A program which initiates the distribution of tracers
      and advects those tracers in a time evolving velocity field.
      Called and used from the CitCOM finite element code.
      Written 2/96 M. Gurnis for Citcom in cartesian geometry
      Modified by Lijie in 1998 and by Vlad and Eh in 2005 for the
      regional version of CitcomS. In 2003, Allen McNamara wrote the
      tracer module for the global version of CitcomS. In 2007, Eh Tan
      merged the two versions of tracer codes together.
*/

/*#define MELT_DEBUG*/
#include <math.h>
#include <string.h>
#include "global_defs.h"
#include "parsing.h"
#include "parallel_related.h"
#include "composition_related.h"
#include "element_definitions.h"

#ifdef USE_GGRD
#include "ggrd_handling.h"
#endif

#ifdef USE_GZDIR
int open_file_zipped(char *, FILE **,struct All_variables *);
void gzip_file(char *);
#endif
void parallel_process_termination();
int icheck_that_processor_shell(struct All_variables *E,
                                       int j, int nprocessor, double rad);
void expand_later_array(struct All_variables *E, int j);
void expand_tracer_arrays(struct All_variables *E, int j);
void tracer_post_processing(struct All_variables *E);
void allocate_tracer_arrays(struct All_variables *E,
                            int j, int number_of_tracers);
void count_tracers_of_flavors(struct All_variables *E);

int full_icheck_cap(struct All_variables *E, int icap,
                    double x, double y, double z, double rad);
int regional_icheck_cap(struct All_variables *E, int icap,
                        double x, double y, double z, double rad);

static float depth_pressure(double depth_citcom);
static double get_temperature(struct All_variables *E,int j, int nelem, double theta, double phi, double rad);
static float degree_of_dry_melting(float T, float P, float Mcpx, float *Ts);
static float katz_wet_melting(float T, float P, float Mcpx, float water);
static float kelley_wetmelting(float T, float P, float water, float Mcpx, int fertility);
static float linear_melting(struct All_variables *E, float T, float rad);
static void find_tracers(struct All_variables *E);
static void predict_tracers(struct All_variables *E);
static void correct_tracers(struct All_variables *E);
static void make_tracer_array(struct All_variables *E);
static void generate_random_tracers(struct All_variables *E,
                                    int tracers_cap, int j);
static void read_tracer_file(struct All_variables *E);
static void read_old_tracer_file(struct All_variables *E);
static void check_sum(struct All_variables *E);
static int isum_tracers(struct All_variables *E);
static void init_tracer_flavors(struct All_variables *E);
static void reduce_tracer_arrays(struct All_variables *E);
static void put_away_later(struct All_variables *E, int j, int it);
static void eject_tracer(struct All_variables *E, int j, int it);
static void correct_ic_sol(struct All_variables *);
static float linear_melting_postp(struct All_variables *, float, float);
int read_double_vector(FILE *, int , double *);
void cart_to_sphere(struct All_variables *,
                    double , double , double ,
                    double *, double *, double *);
void sphere_to_cart(struct All_variables *,
                    double , double , double ,
                    double *, double *, double *);
int icheck_processor_shell(struct All_variables *,
                           int , double );
void update_melt(struct All_variables *);
static double melt_frac_flavor_method(struct All_variables *,int,double,double);
static float U_dot_grad_F(struct All_variables *E,int e,float *F);
static double get_element_area(struct All_variables *E, int m, int iel);


void tracer_input(struct All_variables *E)
{
    void full_tracer_input();
    void myerror();
    void report();
    char message[100];
    int m=E->parallel.me;
	const int Is0 = (E->parallel.me==0);
    int i;
    input_boolean("tracer",&(E->control.tracer),"off",m);
    input_boolean("tracer_enriched",
		  &(E->control.tracer_enriched),"off",m);
	input_int("enriched_flavor",&(E->trace.enriched_flavor),"0",m);
    if(E->control.tracer_enriched){
      if(!E->control.tracer)	/* check here so that we can get away
				   with only one if statement in
				   Advection_diffusion */
	myerror(E,"need to switch on tracers for tracer_enriched");
      
      /*input_float("Q0_enriched",&(E->control.Q0ER),"0.0",m);
      snprintf(message,100,"using compositionally enriched heating: C = 0: %g C = 1: %g (only one composition!)",
	       E->control.Q0,E->control.Q0ER);
      report(E,message);*/
      //
      // this check doesn't work at this point in the code, and we didn't want to put it into every call to
      // Advection diffusion
      //
      //if(E->composition.ncomp != 1)
      //myerror(E,"enriched tracers cannot deal with more than one composition yet");

    }
    if(E->control.tracer) {

        /* tracer_ic_method=0 (random generated array) */
        /* tracer_ic_method=1 (all proc read the same file) */
        /* tracer_ic_method=2 (each proc reads its restart file) */
        input_int("tracer_ic_method",&(E->trace.ic_method),"0,0,nomax",m);

        if (E->trace.ic_method==0){
            input_int("tracers_per_element",&(E->trace.itperel),"10,0,nomax",m);
	}
        else if (E->trace.ic_method==1)
            input_string("tracer_file",E->trace.tracer_file,"tracer.dat",m);
        else if (E->trace.ic_method==2) {
            /* Use 'datadir_old', 'datafile_old', and 'solution_cycles_init' */
            /* to form the filename */
        }
        else {
            fprintf(stderr,"Sorry, tracer_ic_method only 0, 1 and 2 available\n");
            parallel_process_termination();
        }


        /* How many flavors of tracers */
        /* If tracer_flavors > 0, each element will report the number of
         * tracers of each flavor inside it. This information can be used
         * later for many purposes. One of it is to compute composition,
         * either using absolute method or ratio method. */
        input_int("tracer_flavors",&(E->trace.nflavors),"0,0,nomax",m);
		if(E->control.t_ic_sol>0)
			E->trace.nflavors += 1;
		if(Is0) //lhy debug
			fprintf(stderr,"tracer_nflavors = %d\n",E->trace.nflavors);

	/* 0: default from layers 
	   1: from netcdf grds
	   
	   
	   99: from grds, overriding checkpoints during restart
	   (1 and 99 require ggrd)
	*/

        input_int("ic_method_for_flavors",
		  &(E->trace.ic_method_for_flavors),"0,0,nomax",m);


        if (E->trace.nflavors > 1) {
	      /* default method */
          //  case 0:			
	      /* flavors initialized from layers */
			E->trace.z_interface = (double*) malloc((E->trace.nflavors-1)*sizeof(double));
            for(i=0; i<E->trace.nflavors-1; i++){
            	E->trace.z_interface[i] = 0.7;
			}

            input_int("flavor_method",&(E->trace.flavor_method),"0",m);
			input_double_vector("z_interface", E->trace.nflavors-1,E->trace.z_interface, m);
		}
                //break;
		/* 
		   two grd init method, second will override restart
		*/
/*#ifdef USE_GGRD
            case 1:
	    case 99:*/		/* will override restart */
	      /* from grid in top n materials, this will override
		 the checkpoint input */
//	      input_string("ictracer_grd_file",E->trace.ggrd_file,"",m); /* file from which to read */
//	      input_int("ictracer_grd_layers",&(E->trace.ggrd_layers),"2",m);
//									      >0 : which top layers to use, layer <= ictracer_grd_layers
//									      <0 : only use one layer layer == -ictracer_grd_layers
	      //break;
	      
//#endif
 /*           default:
                fprintf(stderr,"ic_method_for_flavors %i undefined (1 and 99 only for ggrd mode)\n",E->trace.ic_method_for_flavors);
                parallel_process_termination();
                break;
*/
        }

        /* Warning level */
        input_boolean("itracer_warnings",&(E->trace.itracer_warnings),"on",m);


        if(E->parallel.nprocxy == 12)
            full_tracer_input(E);


        composition_input(E);
	


    return;
}


void tracer_initial_settings(struct All_variables *E)
{
   void full_keep_within_bounds();
   void full_tracer_setup();
   void full_get_velocity();
   int full_iget_element();
   void regional_keep_within_bounds();
   void regional_tracer_setup();
   void regional_get_velocity();
   int regional_iget_element();

   E->trace.advection_time = 0;
   E->trace.find_tracers_time = 0;
   E->trace.lost_souls_time = 0;

   if(E->parallel.nprocxy == 1) {
       E->problem_tracer_setup = regional_tracer_setup;

       E->trace.keep_within_bounds = regional_keep_within_bounds;
       E->trace.get_velocity = regional_get_velocity;
       E->trace.iget_element = regional_iget_element;
   }
   else {
       E->problem_tracer_setup = full_tracer_setup;

       E->trace.keep_within_bounds = full_keep_within_bounds;
       E->trace.get_velocity = full_get_velocity;
       E->trace.iget_element = full_iget_element;
   }
}



/*****************************************************************************/
/* This function is the primary tracing routine called from Citcom.c         */
/* In this code, unlike the original 3D cartesian code, force is filled      */
/* during Stokes solution. No need to call thermal_buoyancy() after tracing. */


void tracer_advection(struct All_variables *E)
{
    double CPU_time0();
    double begin_time = CPU_time0();
	void calculate_melting_flux();

	if(E->control.t_ic_sol>0)
		correct_ic_sol(E);

    /* advect tracers */
    predict_tracers(E);
    correct_tracers(E);

	//mml melt fraction
	if(E->control.melting_model!=0){
		calculate_melting_flux(E);
	}

    /* check that the number of tracers is conserved */
    check_sum(E);

    /* count # of tracers of each flavor */
    if (E->trace.nflavors > 0)
        count_tracers_of_flavors(E);

    /* update the composition field */
    if (E->composition.on) {
        fill_composition(E);
    }

    if (E->convection.sol_liq){
		update_melt(E);
    }
    E->trace.advection_time += CPU_time0() - begin_time;

    tracer_post_processing(E);

    return;
}

void update_melt(struct All_variables *E){
	/*
	if (E->parallel.me==0){ //debug
		fprintf(stderr,"update_melt begin\n");
	}
	*/
    int j,kk,flavor,e,ele_tracers;
	const int nel = E->lmesh.nel;
    double x,y,z,r,frac,F_new,T;
	int numtracers;
    for (j=1; j<=E->sphere.caps_per_proc; j++) {
		numtracers = E->trace.ntracers[j];
		for (e=1; e<=nel; e++){
			E->melt_el[j][e] = 0.0; 
		}
        /* Fill arrays */
        for (kk=1; kk<=numtracers; kk++) {
            e = E->trace.ielement[j][kk];
            flavor = E->trace.extraq[j][0][kk];
		    T = 0.0; //undifined
            x=E->trace.basicq[j][3][kk];
            y=E->trace.basicq[j][4][kk];
            z=E->trace.basicq[j][5][kk];
			r=E->trace.basicq[j][2][kk];
	   		F_new = melt_frac_flavor_method(E,flavor,T,r); 
	    	if(F_new>E->trace.melt[j][0][kk]){
				E->trace.melt[j][1][kk] = F_new-E->trace.melt[j][0][kk];
				E->trace.melt[j][0][kk] = F_new;
	    	}
	    	else{
				E->trace.melt[j][1][kk] = 0.0;
	    	}
			E->melt_el[j][e] += E->trace.melt[j][1][kk]; 
        }
		/*average to total tracers in an element*/
		for (e=1; e<=nel; e++){
			ele_tracers = 0;
        	for (flavor=0; flavor<E->trace.nflavors; flavor++){
            	ele_tracers += E->trace.ntracer_flavor[j][flavor][e];
			}
			E->melt_el[j][e] = E->melt_el[j][e]/ele_tracers; 
		}
    }	
	/*
	if (E->parallel.me==0){ //debug
		fprintf(stderr,"update_melt end\n");
	}
	*/
}

static double melt_frac_flavor_method(struct All_variables *E,int flavor, double T, double r){
    double T_sol,T_liq,r0,r1;
    double tiny=1e-8,frac;
	int i;
    /*get solidus and liquidus temperature at r*/
    for(i=1;i<=E->lmesh.noz-1;i++)  {
	r0 = E->sx[1][3][i];
 	r1 = E->sx[1][3][i+1];
	if(r<r1+tiny && r>r0-tiny)
	{
            T_sol=(r-r1)/(r0-r1)*E->sol_l[1][i]+(r-r0)/(r1-r0)*E->sol_l[1][i+1];
            T_liq=(r-r1)/(r0-r1)*E->liq_l[1][i]+(r-r0)/(r1-r0)*E->liq_l[1][i+1];
	}
    }
    frac = T>T_sol-tiny?(T-T_sol)/(T_liq-T_sol):0.0;
    return frac;
}

/********* TRACER POST PROCESSING ****************************************/

void tracer_post_processing(struct All_variables *E)
{
    int i;

    /* reset statistical counters */

    E->trace.istat_isend=0;
    E->trace.istat_elements_checked=0;
    E->trace.istat1=0;

    /* write timing information every 20 steps */
    if ((E->monitor.solution_cycles % 20) == 0) {
        fprintf(E->trace.fpt, "STEP %d\n", E->monitor.solution_cycles);

        fprintf(E->trace.fpt, "Advecting tracers takes %f seconds.\n",
                E->trace.advection_time - E->trace.find_tracers_time);
        fprintf(E->trace.fpt, "Finding element takes %f seconds.\n",
                E->trace.find_tracers_time - E->trace.lost_souls_time);
        fprintf(E->trace.fpt, "Exchanging lost tracers takes %f seconds.\n",
                E->trace.lost_souls_time);
    }

    if(E->control.verbose){
      fprintf(E->trace.fpt,"Number of times for all element search  %d\n",E->trace.istat1);

      fprintf(E->trace.fpt,"Number of tracers sent to other processors: %d\n",E->trace.istat_isend);

      fprintf(E->trace.fpt,"Number of times element columns are checked: %d \n",E->trace.istat_elements_checked);

      /* compositional and error fraction data files */
      //TODO: move
      if (E->composition.on) {
        fprintf(E->trace.fpt,"Empty elements filled with old compositional "
                "values: %d (%f percent)\n", E->trace.istat_iempty,
                (100.0*E->trace.istat_iempty)/E->lmesh.nel);
        E->trace.istat_iempty=0;


        get_bulk_composition(E);

        if (E->parallel.me==0) {

            fprintf(E->fp,"composition: %e",E->monitor.elapsed_time);
            for (i=0; i<E->composition.ncomp; i++)
                fprintf(E->fp," %e", E->composition.bulk_composition[i]);
            fprintf(E->fp,"\n");

            fprintf(E->fp,"composition_error_fraction: %e",E->monitor.elapsed_time);
            for (i=0; i<E->composition.ncomp; i++)
                fprintf(E->fp," %e", E->composition.error_fraction[i]);
            fprintf(E->fp,"\n");

        }
      }
      fflush(E->trace.fpt);
    }

    return;
}


/*********** PREDICT TRACERS **********************************************/
/*                                                                        */
/* This function predicts tracers performing an euler step                */
/*                                                                        */
/*                                                                        */
/* Note positions used in tracer array                                    */
/* [positions 0-5 are always fixed with current coordinates               */
/*  Positions 6-8 contain original Cartesian coordinates.                 */
/*  Positions 9-11 contain original Cartesian velocities.                 */
/*                                                                        */


static void predict_tracers(struct All_variables *E)
{

    int numtracers;
    int j;
    int kk;
    int nelem;

    double dt;
    double theta0,phi0,rad0;
    double x0,y0,z0;
    double theta_pred,phi_pred,rad_pred;
    double x_pred,y_pred,z_pred;
    double velocity_vector[4];

    void cart_to_sphere();


    dt=E->advection.timestep;


    for (j=1;j<=E->sphere.caps_per_proc;j++) {

        numtracers=E->trace.ntracers[j];

        for (kk=1;kk<=numtracers;kk++) {

            theta0=E->trace.basicq[j][0][kk];
            phi0=E->trace.basicq[j][1][kk];
            rad0=E->trace.basicq[j][2][kk];
            x0=E->trace.basicq[j][3][kk];
            y0=E->trace.basicq[j][4][kk];
            z0=E->trace.basicq[j][5][kk];

            nelem=E->trace.ielement[j][kk];
            (E->trace.get_velocity)(E,j,nelem,theta0,phi0,rad0,velocity_vector);

            x_pred=x0+velocity_vector[1]*dt;
            y_pred=y0+velocity_vector[2]*dt;
            z_pred=z0+velocity_vector[3]*dt;


            /* keep in box */

            cart_to_sphere(E,x_pred,y_pred,z_pred,&theta_pred,&phi_pred,&rad_pred);
            (E->trace.keep_within_bounds)(E,&x_pred,&y_pred,&z_pred,&theta_pred,&phi_pred,&rad_pred);

            /* Current Coordinates are always kept in positions 0-5. */

            E->trace.basicq[j][0][kk]=theta_pred;
            E->trace.basicq[j][1][kk]=phi_pred;
            E->trace.basicq[j][2][kk]=rad_pred;
            E->trace.basicq[j][3][kk]=x_pred;
            E->trace.basicq[j][4][kk]=y_pred;
            E->trace.basicq[j][5][kk]=z_pred;

            /* Fill in original coords (positions 6-8) */

            E->trace.basicq[j][6][kk]=x0;
            E->trace.basicq[j][7][kk]=y0;
            E->trace.basicq[j][8][kk]=z0;

            /* Fill in original velocities (positions 9-11) */

            E->trace.basicq[j][9][kk]=velocity_vector[1];  /* Vx */
            E->trace.basicq[j][10][kk]=velocity_vector[2];  /* Vy */
            E->trace.basicq[j][11][kk]=velocity_vector[3];  /* Vz */


        } /* end kk, predicting tracers */
    } /* end caps */

    /* find new tracer elements and caps */

    find_tracers(E);

    return;

}


/*********** CORRECT TRACERS **********************************************/
/*                                                                        */
/* This function corrects tracers using both initial and                  */
/* predicted velocities                                                   */
/*                                                                        */
/*                                                                        */
/* Note positions used in tracer array                                    */
/* [positions 0-5 are always fixed with current coordinates               */
/*  Positions 6-8 contain original Cartesian coordinates.                 */
/*  Positions 9-11 contain original Cartesian velocities.                 */
/*                                                                        */


static void correct_tracers(struct All_variables *E)
{

    int j;
    int kk;
    int nelem;


    double dt;
    double x0,y0,z0;
    double theta_pred,phi_pred,rad_pred;
    double x_pred,y_pred,z_pred;
    double theta_cor,phi_cor,rad_cor;
    double x_cor,y_cor,z_cor;
    double velocity_vector[4];
    double Vx0,Vy0,Vz0;
    double Vx_pred,Vy_pred,Vz_pred;

    void cart_to_sphere();


    dt=E->advection.timestep;


    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        for (kk=1;kk<=E->trace.ntracers[j];kk++) {

            theta_pred=E->trace.basicq[j][0][kk];
            phi_pred=E->trace.basicq[j][1][kk];
            rad_pred=E->trace.basicq[j][2][kk];
            x_pred=E->trace.basicq[j][3][kk];
            y_pred=E->trace.basicq[j][4][kk];
            z_pred=E->trace.basicq[j][5][kk];

            x0=E->trace.basicq[j][6][kk];
            y0=E->trace.basicq[j][7][kk];
            z0=E->trace.basicq[j][8][kk];

            Vx0=E->trace.basicq[j][9][kk];
            Vy0=E->trace.basicq[j][10][kk];
            Vz0=E->trace.basicq[j][11][kk];

            nelem=E->trace.ielement[j][kk];

            (E->trace.get_velocity)(E,j,nelem,theta_pred,phi_pred,rad_pred,velocity_vector);

            Vx_pred=velocity_vector[1];
            Vy_pred=velocity_vector[2];
            Vz_pred=velocity_vector[3];

            x_cor=x0 + dt * 0.5*(Vx0+Vx_pred);
            y_cor=y0 + dt * 0.5*(Vy0+Vy_pred);
            z_cor=z0 + dt * 0.5*(Vz0+Vz_pred);

            cart_to_sphere(E,x_cor,y_cor,z_cor,&theta_cor,&phi_cor,&rad_cor);
            (E->trace.keep_within_bounds)(E,&x_cor,&y_cor,&z_cor,&theta_cor,&phi_cor,&rad_cor);

            /* Fill in Current Positions (other positions are no longer important) */

            E->trace.basicq[j][0][kk]=theta_cor;
            E->trace.basicq[j][1][kk]=phi_cor;
            E->trace.basicq[j][2][kk]=rad_cor;
            E->trace.basicq[j][3][kk]=x_cor;
            E->trace.basicq[j][4][kk]=y_cor;
            E->trace.basicq[j][5][kk]=z_cor;

        } /* end kk, correcting tracers */
    } /* end caps */

    /* find new tracer elements and caps */

    find_tracers(E);

    return;
}


/************ FIND TRACERS *************************************/
/*                                                             */
/* This function finds tracer elements and moves tracers to    */
/* other processor domains if necessary.                       */
/* Array ielement is filled with elemental values.                */

static void find_tracers(struct All_variables *E)
{

    int iel;
    int kk;
    int j;
    int it;
    int iprevious_element;
    int num_tracers;

    double theta,phi,rad;
    double x,y,z;
    double time_stat1;
    double time_stat2;

    void put_away_later();
    void eject_tracer();
    void reduce_tracer_arrays();
    void sphere_to_cart();
    void full_lost_souls();
    void regional_lost_souls();

    double CPU_time0();
    double begin_time = CPU_time0();


    for (j=1;j<=E->sphere.caps_per_proc;j++) {


        /* initialize arrays and statistical counters */

        E->trace.ilater[j]=E->trace.ilatersize[j]=0;

        E->trace.istat1=0;
        for (kk=0;kk<=4;kk++) {
            E->trace.istat_ichoice[j][kk]=0;
        }

        //TODO: use while-loop instead of for-loop
        /* important to index by it, not kk */

        it=0;
        num_tracers=E->trace.ntracers[j];

        for (kk=1;kk<=num_tracers;kk++) {

            it++;

            theta=E->trace.basicq[j][0][it];
            phi=E->trace.basicq[j][1][it];
            rad=E->trace.basicq[j][2][it];
            x=E->trace.basicq[j][3][it];
            y=E->trace.basicq[j][4][it];
            z=E->trace.basicq[j][5][it];

            iprevious_element=E->trace.ielement[j][it];

            iel=(E->trace.iget_element)(E,j,iprevious_element,x,y,z,theta,phi,rad);
            /* debug *
            fprintf(E->trace.fpt,"BB. kk %d %d %d %d %f %f %f %f %f %f\n",kk,j,iprevious_element,iel,x,y,z,theta,phi,rad);
            fflush(E->trace.fpt);
            */

            E->trace.ielement[j][it]=iel;

            if (iel == -99) {
                /* tracer is inside other processors */
                put_away_later(E,j,it);
                eject_tracer(E,j,it);
                it--;
            } else if (iel == -1) {
                /* tracer is inside this processor,
                 * but cannot find its element.
                 * Throw away the tracer. */

                if (E->trace.itracer_warnings) exit(10);


                eject_tracer(E,j,it);
                it--;
            }

        } /* end tracers */

    } /* end j */


    /* Now take care of tracers that exited cap */

    /* REMOVE */
    /*
      parallel_process_termination();
    */

    if (E->parallel.nprocxy == 12)
        full_lost_souls(E);
    else
        regional_lost_souls(E);

    /* Free later arrays */

    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        if (E->trace.ilatersize[j]>0) {
            for (kk=0;kk<=((E->trace.number_of_tracer_quantities)-1);kk++) {
                free(E->trace.rlater[j][kk]);
            }
        }
    } /* end j */


    /* Adjust Array Sizes */

    reduce_tracer_arrays(E);

    E->trace.find_tracers_time += CPU_time0() - begin_time;

    return;
}


/***********************************************************************/
/* This function computes the number of tracers in each element.       */
/* Each tracer can be of different "flavors", which is the 0th index   */
/* of extraq. How to interprete "flavor" is left for the application.  */

void count_tracers_of_flavors(struct All_variables *E)
{

    int j, flavor, e, kk;
    int numtracers;

    for (j=1; j<=E->sphere.caps_per_proc; j++) {

        /* first zero arrays */
        for (flavor=0; flavor<E->trace.nflavors; flavor++)
            for (e=1; e<=E->lmesh.nel; e++)
                E->trace.ntracer_flavor[j][flavor][e] = 0;

        numtracers=E->trace.ntracers[j];

        /* Fill arrays */
        for (kk=1; kk<=numtracers; kk++) {
            e = E->trace.ielement[j][kk];
            flavor = E->trace.extraq[j][0][kk];
            E->trace.ntracer_flavor[j][flavor][e]++;
        }
    }

    /* debug */
    /**
    for (j=1; j<=E->sphere.caps_per_proc; j++) {
        for (e=1; e<=E->lmesh.nel; e++) {
            fprintf(E->trace.fpt, "element=%d ntracer_flaver =", e);
            for (flavor=0; flavor<E->trace.nflavors; flavor++) {
                fprintf(E->trace.fpt, " %d",
                        E->trace.ntracer_flavor[j][flavor][e]);
            }
            fprintf(E->trace.fpt, "\n");
        }
    }
    fflush(E->trace.fpt);
    */

    return;
}



void initialize_tracers(struct All_variables *E)
{

    if (E->trace.ic_method==0)
        make_tracer_array(E);
    else if (E->trace.ic_method==1)
        read_tracer_file(E);
    else if (E->trace.ic_method==2)
        read_old_tracer_file(E);
    else {
        fprintf(E->trace.fpt,"Not ready for other inputs yet\n");
        fflush(E->trace.fpt);
        parallel_process_termination();
    }
	
    /* total number of tracers  */

    E->trace.ilast_tracer_count = isum_tracers(E);
    fprintf(E->trace.fpt, "Sum of Tracers: %d\n", E->trace.ilast_tracer_count);
    if(E->parallel.me==0)
        fprintf(stderr, "Sum of Tracers: %d\n", E->trace.ilast_tracer_count);


    /* find elements */

    //find_tracers(E);


    /* count # of tracers of each flavor */

    if (E->trace.nflavors > 0)
        count_tracers_of_flavors(E);

    return;
}


/************** MAKE TRACER ARRAY ********************************/
/* Here, each processor will generate tracers somewhere          */
/* in the sphere - check if its in this cap  - then check radial */

static void make_tracer_array(struct All_variables *E)
{

    int tracers_cap;
    int j;
    double processor_fraction;

    void generate_random_tracers();
    void init_tracer_flavors();

    if (E->parallel.me==0) fprintf(stderr,"Making Tracer Array\n");

    for (j=1;j<=E->sphere.caps_per_proc;j++) {

        processor_fraction=E->lmesh.volume/E->mesh.volume;
        tracers_cap=E->mesh.nel*E->trace.itperel*processor_fraction;
        /*
          fprintf(stderr,"AA: proc frac: %f (%d) %d %d %f %f\n",processor_fraction,tracers_cap,E->lmesh.nel,E->parallel.nprocz, E->sx[j][3][E->lmesh.noz],E->sx[j][3][1]);
        */

        fprintf(E->trace.fpt,"\nGenerating %d Tracers\n",tracers_cap);

        generate_random_tracers(E, tracers_cap, j);



    }/* end j */
	find_tracers(E);


    /* Initialize tracer flavors */
    if (E->trace.nflavors) init_tracer_flavors(E);

    return;
}



static void generate_random_tracers(struct All_variables *E,
                                    int tracers_cap, int j)
{
    void cart_to_sphere();
    int kk;
    int ival;
    int number_of_tries=0;
    int max_tries;

    double x,y,z;
    double theta,phi,rad;
    double xmin,xmax,ymin,ymax,zmin,zmax;
    double random1,random2,random3;


    allocate_tracer_arrays(E,j,tracers_cap);

    /* Finding the min/max of the cartesian coordinates. */
    /* One must loop over E->X to find the min/max, since the 8 corner */
    /* nodes may not be the min/max. */
    xmin = ymin = zmin = E->sphere.ro;
    xmax = ymax = zmax = -E->sphere.ro;
    for (kk=1; kk<=E->lmesh.nno; kk++) {
        x = E->x[j][1][kk];
        y = E->x[j][2][kk];
        z = E->x[j][3][kk];

        xmin = ((xmin < x) ? xmin : x);
        xmax = ((xmax > x) ? xmax : x);
        ymin = ((ymin < y) ? ymin : y);
        ymax = ((ymax > y) ? ymax : y);
        zmin = ((zmin < z) ? zmin : z);
        zmax = ((zmax > z) ? zmax : z);
    }

    /* Tracers are placed randomly in cap */
    /* (intentionally using rand() instead of srand() )*/
    while (E->trace.ntracers[j]<tracers_cap) {

        number_of_tries++;
        max_tries=100*tracers_cap;

        if (number_of_tries>max_tries) {
            fprintf(E->trace.fpt,"Error(make_tracer_array)-too many tries?\n");
            fprintf(E->trace.fpt,"%d %d %d\n",max_tries,number_of_tries,RAND_MAX);
            fflush(E->trace.fpt);
            exit(10);
        }

#if 1
        random1=drand48();
        random2=drand48();
        random3=drand48();
#else  /* never called */
        random1=(1.0*rand())/(1.0*RAND_MAX);
        random2=(1.0*rand())/(1.0*RAND_MAX);
        random3=(1.0*rand())/(1.0*RAND_MAX);
#endif

        x=xmin+random1*(xmax-xmin);
        y=ymin+random2*(ymax-ymin);
        z=zmin+random3*(zmax-zmin);

        /* first check if within shell */

        cart_to_sphere(E,x,y,z,&theta,&phi,&rad);

        if (rad>=E->sx[j][3][E->lmesh.noz]) continue;
        if (rad<E->sx[j][3][1]) continue;


        /* check if in current cap */
        if (E->parallel.nprocxy==1)
            ival=regional_icheck_cap(E,0,theta,phi,rad,rad);
        else
            ival=full_icheck_cap(E,0,x,y,z,rad);

        if (ival!=1) continue;

        /* Made it, so record tracer information */

        (E->trace.keep_within_bounds)(E,&x,&y,&z,&theta,&phi,&rad);

        E->trace.ntracers[j]++;
        kk=E->trace.ntracers[j];

        E->trace.basicq[j][0][kk]=theta;
        E->trace.basicq[j][1][kk]=phi;
        E->trace.basicq[j][2][kk]=rad;
        E->trace.basicq[j][3][kk]=x;
        E->trace.basicq[j][4][kk]=y;
        E->trace.basicq[j][5][kk]=z;

    } /* end while */

    return;
}


/******** READ TRACER ARRAY *********************************************/
/*                                                                      */
/* This function reads tracers from input file.                         */
/* All processors read the same input file, then sort out which ones    */
/* belong.                                                              */

static void read_tracer_file(struct All_variables *E)
{

    char input_s[1000];

    int number_of_tracers, ncolumns;
    int kk;
    int icheck;
    int iestimate;
    int icushion;
    int i, j;


    int icheck_processor_shell();
    void sphere_to_cart();
    void cart_to_sphere();
    void expand_tracer_arrays();

    double x,y,z;
    double theta,phi,rad;
    double buffer[100];

    FILE *fptracer;

    fptracer=fopen(E->trace.tracer_file,"r");

    fgets(input_s,200,fptracer);
    if(sscanf(input_s,"%d %d",&number_of_tracers,&ncolumns) != 2) {
        fprintf(stderr,"Error while reading file '%s'\n", E->trace.tracer_file);
        exit(8);
    }
    fprintf(E->trace.fpt,"%d Tracers, %d columns in file \n",
            number_of_tracers, ncolumns);

    /* some error control */
    if (E->trace.number_of_extra_quantities+3 != ncolumns) {
        fprintf(E->trace.fpt,"ERROR(read tracer file)-wrong # of columns\n");
        fflush(E->trace.fpt);
        exit(10);
    }


    /* initially size tracer arrays to number of tracers divided by processors */

    icushion=100;

    /* for absolute tracer method */
    E->trace.number_of_tracers = number_of_tracers;

    iestimate=number_of_tracers/E->parallel.nproc + icushion;

    for (j=1;j<=E->sphere.caps_per_proc;j++) {

        allocate_tracer_arrays(E,j,iestimate);

        for (kk=1;kk<=number_of_tracers;kk++) {
            int len, ncol;
            ncol = 3 + E->trace.number_of_extra_quantities;

            len = read_double_vector(fptracer, ncol, buffer);
            if (len != ncol) {
                fprintf(E->trace.fpt,"ERROR(read tracer file) - wrong input file format: %s\n", E->trace.tracer_file);
                fflush(E->trace.fpt);
                exit(10);
            }

            theta = buffer[0];
            phi = buffer[1];
            rad = buffer[2];

            sphere_to_cart(E,theta,phi,rad,&x,&y,&z);


            /* make sure theta, phi is in range, and radius is within bounds */

            (E->trace.keep_within_bounds)(E,&x,&y,&z,&theta,&phi,&rad);

            /* check whether tracer is within processor domain */

            icheck=1;
            if (E->parallel.nprocz>1) icheck=icheck_processor_shell(E,j,rad);
            if (icheck!=1) continue;

            if (E->parallel.nprocxy==1)
                icheck=regional_icheck_cap(E,0,theta,phi,rad,rad);
            else
                icheck=full_icheck_cap(E,0,x,y,z,rad);

            if (icheck==0) continue;

            /* if still here, tracer is in processor domain */


            E->trace.ntracers[j]++;

            if (E->trace.ntracers[j]>=(E->trace.max_ntracers[j]-5)) expand_tracer_arrays(E,j);

            E->trace.basicq[j][0][E->trace.ntracers[j]]=theta;
            E->trace.basicq[j][1][E->trace.ntracers[j]]=phi;
            E->trace.basicq[j][2][E->trace.ntracers[j]]=rad;
            E->trace.basicq[j][3][E->trace.ntracers[j]]=x;
            E->trace.basicq[j][4][E->trace.ntracers[j]]=y;
            E->trace.basicq[j][5][E->trace.ntracers[j]]=z;

            for (i=0; i<E->trace.number_of_extra_quantities; i++)
                E->trace.extraq[j][i][E->trace.ntracers[j]]=buffer[i+3];

        } /* end kk, number of tracers */

        fprintf(E->trace.fpt,"Number of tracers in this cap is: %d\n",
                E->trace.ntracers[j]);

        /** debug **
        for (kk=1; kk<=E->trace.ntracers[j]; kk++) {
            fprintf(E->trace.fpt, "tracer#=%d sph_coord=(%g,%g,%g)", kk,
                    E->trace.basicq[j][0][kk],
                    E->trace.basicq[j][1][kk],
                    E->trace.basicq[j][2][kk]);
            fprintf(E->trace.fpt, "   extraq=");
            for (i=0; i<E->trace.number_of_extra_quantities; i++)
                fprintf(E->trace.fpt, " %g", E->trace.extraq[j][i][kk]);
            fprintf(E->trace.fpt, "\n");
        }
        fflush(E->trace.fpt);
        */

    } /* end j */

    fclose(fptracer);

    icheck=isum_tracers(E);

    if (icheck!=number_of_tracers) {
        fprintf(E->trace.fpt,"ERROR(read_tracer_file) - tracers != number in file\n");
        fprintf(E->trace.fpt,"Tracers in system: %d\n", icheck);
        fprintf(E->trace.fpt,"Tracers in file: %d\n", number_of_tracers);
        fflush(E->trace.fpt);
        exit(10);
    }

    return;
}


/************** READ OLD TRACER FILE *************************************/
/*                                                                       */
/* This function read tracers written from previous calculation          */
/* and the tracers are read as seperate files for each processor domain. */

static void read_old_tracer_file(struct All_variables *E)
{

    char output_file[200];
    char input_s[1000];

    int i,j,kk,rezip;
    int idum1,ncolumns;
    int numtracers;

    double rdum1;
    double theta,phi,rad;
    double x,y,z;
    double buffer[100];

    void sphere_to_cart();

    FILE *fp1;

    if (E->trace.number_of_extra_quantities>99) {
        fprintf(E->trace.fpt,"ERROR(read_old_tracer_file)-increase size of extra[]\n");
        fflush(E->trace.fpt);
        parallel_process_termination();
    }



    /* deal with different output formats */
#ifdef USE_GZDIR
    if(strcmp(E->output.format, "ascii-gz") == 0){
      sprintf(output_file,"%s/%d/tracer.%d.%d",
	      E->control.data_dir_old,E->monitor.solution_cycles_init,E->parallel.me,E->monitor.solution_cycles_init);
      rezip = open_file_zipped(output_file,&fp1,E);
    }else{
      sprintf(output_file,"%s.tracer.%d.%d",E->control.old_P_file,E->parallel.me,E->monitor.solution_cycles_init);
      if ( (fp1=fopen(output_file,"r"))==NULL) {
        fprintf(E->trace.fpt,"ERROR(read_old_tracer_file)-gziped file not found %s\n",output_file);
        fflush(E->trace.fpt);
        exit(10);
      }
    }
#else
    sprintf(output_file,"%s.tracer.%d.%d",E->control.old_P_file,E->parallel.me,E->monitor.solution_cycles_init);
    if ( (fp1=fopen(output_file,"r"))==NULL) {
        fprintf(E->trace.fpt,"ERROR(read_old_tracer_file)-file not found %s\n",output_file);
        fflush(E->trace.fpt);
        exit(10);
    }
#endif

    fprintf(stderr,"Read old tracers from %s\n",output_file);


    for(j=1;j<=E->sphere.caps_per_proc;j++) {
        fgets(input_s,200,fp1);
        if(sscanf(input_s,"%d %d %d %lf",
                  &idum1, &numtracers, &ncolumns, &rdum1) != 4) {
            fprintf(stderr,"Error while reading file '%s'\n", output_file);
            exit(8);
        }


        /* some error control */
        if (E->trace.number_of_extra_quantities+3 != ncolumns) {
            fprintf(E->trace.fpt,"ERROR(read_old_tracer_file)-wrong # of columns\n");
            fflush(E->trace.fpt);
            exit(10);
        }

        /* allocate memory for tracer arrays */

        allocate_tracer_arrays(E,j,numtracers);
        E->trace.ntracers[j]=numtracers;

        for (kk=1;kk<=numtracers;kk++) {
            int len, ncol;
            ncol = 3 + E->trace.number_of_extra_quantities;

            len = read_double_vector(fp1, ncol, buffer);
            if (len != ncol) {
                fprintf(E->trace.fpt,"ERROR(read_old_tracer_file) - wrong input file format: %s\n", output_file);
                fflush(E->trace.fpt);
                exit(10);
            }

            theta = buffer[0];
            phi = buffer[1];
            rad = buffer[2];

            sphere_to_cart(E,theta,phi,rad,&x,&y,&z);

            /* it is possible that if on phi=0 boundary, significant digits can push phi over 2pi */

            (E->trace.keep_within_bounds)(E,&x,&y,&z,&theta,&phi,&rad);

            E->trace.basicq[j][0][kk]=theta;
            E->trace.basicq[j][1][kk]=phi;
            E->trace.basicq[j][2][kk]=rad;
            E->trace.basicq[j][3][kk]=x;
            E->trace.basicq[j][4][kk]=y;
            E->trace.basicq[j][5][kk]=z;

            for (i=0; i<E->trace.number_of_extra_quantities; i++)
                E->trace.extraq[j][i][kk]=buffer[i+3];

        }

        /** debug **
        for (kk=1; kk<=E->trace.ntracers[j]; kk++) {
            fprintf(E->trace.fpt, "tracer#=%d sph_coord=(%g,%g,%g)", kk,
                    E->trace.basicq[j][0][kk],
                    E->trace.basicq[j][1][kk],
                    E->trace.basicq[j][2][kk]);
            fprintf(E->trace.fpt, "   extraq=");
            for (i=0; i<E->trace.number_of_extra_quantities; i++)
                fprintf(E->trace.fpt, " %g", E->trace.extraq[j][i][kk]);
            fprintf(E->trace.fpt, "\n");
        }
        fflush(E->trace.fpt);
        */

        fprintf(E->trace.fpt,"Read %d tracers from file %s\n",numtracers,output_file);
        fflush(E->trace.fpt);

    }
    fclose(fp1);
#ifdef USE_GZDIR
    if(strcmp(E->output.format, "ascii-gz") == 0)
      if(rezip)			/* rezip */
	gzip_file(output_file);
#endif

    return;
}





/*********** CHECK SUM **************************************************/
/*                                                                      */
/* This functions checks to make sure number of tracers is preserved    */

static void check_sum(struct All_variables *E)
{

    int number, iold_number;

    number = isum_tracers(E);

    iold_number = E->trace.ilast_tracer_count;

    if (number != iold_number) {
        fprintf(E->trace.fpt,"ERROR(check_sum)-break in conservation %d %d\n",
                number,iold_number);
        fflush(E->trace.fpt);
        if (E->trace.itracer_warnings)
            parallel_process_termination();
    }

    E->trace.ilast_tracer_count = number;

    return;
}


/************* ISUM TRACERS **********************************************/
/*                                                                       */
/* This function uses MPI to sum all tracers and returns number of them. */

static int isum_tracers(struct All_variables *E)
{
    int imycount;
    int iallcount;
    int j;

    iallcount = 0;

    imycount = 0;
    for (j=1; j<=E->sphere.caps_per_proc; j++)
        imycount = imycount + E->trace.ntracers[j];
	
	//fprintf(stderr,"proc = %d, imycount = %d\n",E->parallel.me,imycount);//debug
    MPI_Allreduce(&imycount,&iallcount,1,MPI_INT,MPI_SUM,E->parallel.world);

    return iallcount;
}



/********** CART TO SPHERE ***********************/
void cart_to_sphere(struct All_variables *E,
                    double x, double y, double z,
                    double *theta, double *phi, double *rad)
{

    double temp;
    double myatan();

    temp=x*x+y*y;

    *rad=sqrt(temp+z*z);
    *theta=atan2(sqrt(temp),z);
    *phi=myatan(y,x);
	
    return;
}

/********** SPHERE TO CART ***********************/
void sphere_to_cart(struct All_variables *E,
                    double theta, double phi, double rad,
                    double *x, double *y, double *z)
{

    double sint,cost,sinf,cosf;
    double temp;

    sint=sin(theta);
    cost=cos(theta);
    sinf=sin(phi);
    cosf=cos(phi);

    temp=rad*sint;

    *x=temp*cosf;
    *y=temp*sinf;
    *z=rad*cost;

    return;
}



static void init_tracer_flavors(struct All_variables *E)
{
	void read_chemicals_from_files();
	const int ncomp = E->composition.ncomp;
	const int nel = E->lmesh.nel;
	const int Is0 = (E->parallel.me==0);
	const int Is1 = (E->parallel.me==1);
    int i,j,e,kk, number_of_tracers;
	int flavor;
    double temp, temp1;
    double rad;
	double tiny = 1e-15;
	int *nf_comp;
	float **comp_e;
	const int ic_flavor = E->trace.ic_flavor;
	const double ic_inter =  E->chemical.interface[ic_flavor][1];
	const double ic_thick = E->chemical.interface[ic_flavor][2] 
								- E->chemical.interface[ic_flavor][1];
	const double t_scale = pow(E->data.radius,2.0)/E->data.therm_diff/(365*24*3600*1e6);
	if(Is0){
		fprintf(stderr,"ic_flavor = %d, ic_inter = %.4e, ic_thickness = %.4e, t_scale = %.4e\n", ic_flavor,ic_inter,ic_thick,t_scale);
	}

    switch(E->trace.ic_method_for_flavors){
    case 0:
      /* ic_method_for_flavors == 0 (layered structure) */
      /* any tracer above z_interface[i] is of flavor i */
      /* any tracer below z_interface is of flavor (nflavors-1) */
      for (j=1;j<=E->sphere.caps_per_proc;j++) {

	number_of_tracers = E->trace.ntracers[j];
	for (kk=1;kk<=number_of_tracers;kk++) {
	  rad = E->trace.basicq[j][2][kk];

          //flavor = E->trace.nflavors - 1;
          flavor = 0; //modified by lhy

          //for (i=0; i<E->trace.nflavors-1; i++) 
          for (i=1; i<E->trace.nflavors; i++) { //modified by lhy
			  switch(E->trace.flavor_method){ //lhy 170918 for new flavor method
              case 0:
			  if (rad > E->trace.z_interface[i-1] && rad < (1.0 - E->viscosity.zlith))
                  flavor = i; //lhy for intermediate layer
			  break;
			  case 1:
			  if (rad > E->trace.z_interface[i-1]) 
                  flavor = i; //lhy for intermediate layer
              break;
			  case 2:
			  if (rad < E->trace.z_interface[i-1]) 
                  flavor = i; //lhy for intermediate layer
			  break;
			  case 3:
			  if ((rad > E->chemical.interface[i][1] - tiny)&&(rad < E->chemical.interface[i][2] + tiny)){
				  flavor = i;
				  if((E->control.t_ic_sol>0.0)&&(i==ic_flavor)){
					  flavor = E->trace.nflavors-1;
          			  E->trace.extraq[j][6][kk] = (rad-ic_inter)/ic_thick*E->control.t_ic_sol/t_scale;
					  /*if(Is1){ //lhy debug
						  fprintf(stderr,"flavor = %d, rad = %.4e, t_sol = %.4e\n",flavor,rad,E->trace.extraq[j][6][kk]);
					  }*/
				  }
			  }
		  	  }
          E->trace.extraq[j][0][kk] = flavor;
		 }
    	}
	  }
      break;

    case 2:			/* from grd in top n layers */
		nf_comp = (int*)malloc((ncomp+1)*sizeof(int));
		nf_comp[0] = 1; //mantle
		nf_comp[1] = 3; //ilmenite
		nf_comp[2] = 0; //crust
		comp_e = (float**)malloc((ncomp+1)*sizeof(float*));
		for(i=0;i<=ncomp;i++){
			comp_e[i] = (float*)malloc((nel+1)*sizeof(float));
		}
		//read composition from file
		read_chemicals_from_files(E,comp_e,nf_comp);
        number_of_tracers=E->trace.ntracers[1];
        /*assign flavor to each tracer */
		if(E->parallel.me==0){
			fprintf(stderr,"assign flavor to tracers\n");//debug
		}
        for (kk=1; kk<=number_of_tracers; kk++) {
            e = E->trace.ielement[1][kk];
			flavor = 0;
			temp = comp_e[0][e];
			temp1 = comp_e[0][e];
			for(i=1;i<=ncomp;i++){
				if(comp_e[i][e]>temp){
					temp = comp_e[i][e];
					flavor = i;
				}
				temp1 += comp_e[i][e];
			}
			if(1-temp1>temp){
				flavor = 0;
			}
          	E->trace.extraq[1][0][kk] = flavor;
        }
		for(i=0;i<=ncomp;i++){
			free(comp_e[i]);
		}
		free(comp_e);
		if(E->parallel.me==0){
			fprintf(stderr,"init_tracer_flavors done\n");//debug
		}

		break;
    case 99:			/* (will override restart) */
#ifndef USE_GGRD
      fprintf(stderr,"ic_method_for_flavors %i requires the ggrd routines from hc, -DUSE_GGRD\n",
	      E->trace.ic_method_for_flavors);
      parallel_process_termination();
#else
      ggrd_init_tracer_flavors(E);
#endif
      break;


    default:

      fprintf(stderr,"ic_method_for_flavors %i undefined\n",E->trace.ic_method_for_flavors);
      parallel_process_termination();
      break;
	}	
    return;
}


/******************* get_neighboring_caps ************************************/
/*                                                                           */
/* Communicate with neighboring processors to get their cap boundaries,      */
/* which is later used by (E->trace.icheck_cap)()                            */
/*                                                                           */

void get_neighboring_caps(struct All_variables *E)
{
    void sphere_to_cart();

    const int ncorners = 4; /* # of top corner nodes */
    int i, j, n, d, kk, lev, idb;
    int num_ngb, neighbor_proc, tag;
    MPI_Status status[200];
    MPI_Request request[200];

    int node[ncorners];
    double xx[ncorners*2], rr[12][ncorners*2];
    int nox,noy,noz;
    double x,y,z;
    double theta,phi,rad;

    nox=E->lmesh.nox;
    noy=E->lmesh.noy;
    noz=E->lmesh.noz;

    node[0]=nox*noz*(noy-1)+noz;
    node[1]=noz;
    node[2]=noz*nox;
    node[3]=noz*nox*noy;

    lev = E->mesh.levmax;
    tag = 45;

    for (j=1; j<=E->sphere.caps_per_proc; j++) {

        /* loop over top corners to get their coordinates */
        n = 0;
        for (i=0; i<ncorners; i++) {
            for (d=0; d<2; d++) {
                xx[n] = E->sx[j][d+1][node[i]];
                n++;
            }
        }

        idb = 0;
        num_ngb = E->parallel.TNUM_PASS[lev][j];
        for (kk=1; kk<=num_ngb; kk++) {
            neighbor_proc = E->parallel.PROCESSOR[lev][j].pass[kk];

            MPI_Isend(xx, n, MPI_DOUBLE, neighbor_proc,
                      tag, E->parallel.world, &request[idb]);
            idb++;

            MPI_Irecv(rr[kk], n, MPI_DOUBLE, neighbor_proc,
                      tag, E->parallel.world, &request[idb]);
            idb++;
        }

        /* Storing the current cap information */
        for (i=0; i<n; i++)
            rr[0][i] = xx[i];

        /* Wait for non-blocking calls to complete */

        MPI_Waitall(idb, request, status);

        /* Storing the received cap information
         * XXX: this part assumes:
         *      1) E->sphere.caps_per_proc==1
         *      2) E->mesh.nsd==3
         */
        for (kk=0; kk<=num_ngb; kk++) {
            n = 0;
            for (i=1; i<=ncorners; i++) {
                theta = rr[kk][n++];
                phi = rr[kk][n++];
                rad = E->sphere.ro;

                sphere_to_cart(E, theta, phi, rad, &x, &y, &z);

                E->trace.xcap[kk][i] = x;
                E->trace.ycap[kk][i] = y;
                E->trace.zcap[kk][i] = z;
                E->trace.theta_cap[kk][i] = theta;
                E->trace.phi_cap[kk][i] = phi;
                E->trace.rad_cap[kk][i] = rad;
                E->trace.cos_theta[kk][i] = cos(theta);
                E->trace.sin_theta[kk][i] = sin(theta);
                E->trace.cos_phi[kk][i] = cos(phi);
                E->trace.sin_phi[kk][i] = sin(phi);
            }
        } /* end kk, number of neighbors */

        /* debugging output *
        for (kk=0; kk<=num_ngb; kk++) {
            if (kk==0)
                neighbor_proc = E->parallel.me;
            else
                neighbor_proc = E->parallel.PROCESSOR[lev][1].pass[kk];

            for (i=1; i<=ncorners; i++) {
                fprintf(E->trace.fpt, "pass=%d rank=%d corner=%d "
                        "sx=(%g, %g, %g)\n",
                        kk, neighbor_proc, i,
                        E->trace.theta_cap[kk][i],
                        E->trace.phi_cap[kk][i],
                        E->trace.rad_cap[kk][i]);
            }
        }
        fflush(E->trace.fpt);
        */
    }

    return;
}


/**************** INITIALIZE TRACER ARRAYS ************************************/
/*                                                                            */
/* This function allocates memories to tracer arrays.                         */

void allocate_tracer_arrays(struct All_variables *E,
                            int j, int number_of_tracers)
{

    int kk;
    int ii;

    /* max_ntracers is physical size of tracer array */
    /* (initially make it 25% larger than required */

    E->trace.max_ntracers[j]=number_of_tracers+number_of_tracers/4;
    E->trace.ntracers[j]=0;

    /* make tracer arrays */

    if ((E->trace.ielement[j]=(int *) malloc(E->trace.max_ntracers[j]*sizeof(int)))==NULL) {
        fprintf(E->trace.fpt,"ERROR(make tracer array)-no memory 1a\n");
        fflush(E->trace.fpt);
        exit(10);
    }
    for (kk=1;kk<E->trace.max_ntracers[j];kk++)
        E->trace.ielement[j][kk]=-99;


    for (kk=0;kk<E->trace.number_of_basic_quantities;kk++) {
        if ((E->trace.basicq[j][kk]=(double *)malloc(E->trace.max_ntracers[j]*sizeof(double)))==NULL) {
            fprintf(E->trace.fpt,"ERROR(initialize tracer arrays)-no memory 1b.%d\n",kk);
            fflush(E->trace.fpt);
            exit(10);
        }
    }

    for (kk=0;kk<E->trace.number_of_extra_quantities;kk++) {
        if ((E->trace.extraq[j][kk]=(double *)malloc(E->trace.max_ntracers[j]*sizeof(double)))==NULL) {
            fprintf(E->trace.fpt,"ERROR(initialize tracer arrays)-no memory 1c.%d\n",kk);
            fflush(E->trace.fpt);
            exit(10);
        }
    }
    for (kk=0;kk<=1;kk++) {
	//lhy tracer melt
        if ((E->trace.melt[j][kk]=(double *)malloc(E->trace.max_ntracers[j]*sizeof(double)))==NULL) {
            fprintf(E->trace.fpt,"ERROR(initialize tracer arrays)-no memory 1c.%d\n",kk);
            fflush(E->trace.fpt);
            exit(10);
        }
	else{
	    for(ii=0;ii<E->trace.max_ntracers[j];ii++){
		E->trace.melt[j][kk][ii] = 0.0;
	    }
	}
	}

    if (E->trace.nflavors > 0) {
        E->trace.ntracer_flavor[j]=(int **)malloc(E->trace.nflavors*sizeof(int*));
        for (kk=0;kk<E->trace.nflavors;kk++) {
            if ((E->trace.ntracer_flavor[j][kk]=(int *)malloc((E->lmesh.nel+1)*sizeof(int)))==NULL) {
                fprintf(E->trace.fpt,"ERROR(initialize tracer arrays)-no memory 1c.%d\n",kk);
                fflush(E->trace.fpt);
                exit(10);
            }
        }
    }


    fprintf(E->trace.fpt,"Physical size of tracer arrays (max_ntracers): %d\n",
            E->trace.max_ntracers[j]);
    fflush(E->trace.fpt);

    return;
}



/****** EXPAND TRACER ARRAYS *****************************************/

void expand_tracer_arrays(struct All_variables *E, int j)
{

    int inewsize;
    int kk;
    int icushion;

    /* expand basicq and ielement by 20% */

    icushion=100;

    inewsize=E->trace.max_ntracers[j]+E->trace.max_ntracers[j]/5+icushion;

    if ((E->trace.ielement[j]=(int *)realloc(E->trace.ielement[j],inewsize*sizeof(int)))==NULL) {
        fprintf(E->trace.fpt,"ERROR(expand tracer arrays )-no memory (ielement)\n");
        fflush(E->trace.fpt);
        exit(10);
    }

    for (kk=0;kk<=((E->trace.number_of_basic_quantities)-1);kk++) {
        if ((E->trace.basicq[j][kk]=(double *)realloc(E->trace.basicq[j][kk],inewsize*sizeof(double)))==NULL) {
            fprintf(E->trace.fpt,"ERROR(expand tracer arrays )-no memory (%d)\n",kk);
            fflush(E->trace.fpt);
            exit(10);
        }
    }

    for (kk=0;kk<=((E->trace.number_of_extra_quantities)-1);kk++) {
        if ((E->trace.extraq[j][kk]=(double *)realloc(E->trace.extraq[j][kk],inewsize*sizeof(double)))==NULL) {
            fprintf(E->trace.fpt,"ERROR(expand tracer arrays )-no memory 78 (%d)\n",kk);
            fflush(E->trace.fpt);
            exit(10);
        }
    }


    fprintf(E->trace.fpt,"Expanding physical memory of ielement, basicq, and extraq to %d from %d\n",
            inewsize,E->trace.max_ntracers[j]);

    E->trace.max_ntracers[j]=inewsize;

    return;
}




/****** REDUCE  TRACER ARRAYS *****************************************/

static void reduce_tracer_arrays(struct All_variables *E)
{

    int inewsize;
    int kk;
    int iempty_space;
    int j;

    int icushion=100;

    for (j=1;j<=E->sphere.caps_per_proc;j++) {


        /* if physical size is double tracer size, reduce it */

        iempty_space=(E->trace.max_ntracers[j]-E->trace.ntracers[j]);

        if (iempty_space>(E->trace.ntracers[j]+icushion)) {


            inewsize=E->trace.ntracers[j]+E->trace.ntracers[j]/4+icushion;

            if (inewsize<1) {
                fprintf(E->trace.fpt,"Error(reduce tracer arrays)-something up (hdf3)\n");
                fflush(E->trace.fpt);
                exit(10);
            }


            if ((E->trace.ielement[j]=(int *)realloc(E->trace.ielement[j],inewsize*sizeof(int)))==NULL) {
                fprintf(E->trace.fpt,"ERROR(reduce tracer arrays )-no memory (ielement)\n");
                fflush(E->trace.fpt);
                exit(10);
            }


            for (kk=0;kk<=((E->trace.number_of_basic_quantities)-1);kk++) {
                if ((E->trace.basicq[j][kk]=(double *)realloc(E->trace.basicq[j][kk],inewsize*sizeof(double)))==NULL) {
                    fprintf(E->trace.fpt,"AKM(reduce tracer arrays )-no memory (%d)\n",kk);
                    fflush(E->trace.fpt);
                    exit(10);
                }
            }

            for (kk=0;kk<=((E->trace.number_of_extra_quantities)-1);kk++) {
                if ((E->trace.extraq[j][kk]=(double *)realloc(E->trace.extraq[j][kk],inewsize*sizeof(double)))==NULL) {
                    fprintf(E->trace.fpt,"AKM(reduce tracer arrays )-no memory 783 (%d)\n",kk);
                    fflush(E->trace.fpt);
                    exit(10);
                }
            }


            fprintf(E->trace.fpt,"Reducing physical memory of ielement, basicq, and extraq to %d from %d\n",
                    E->trace.max_ntracers[j],inewsize);

            E->trace.max_ntracers[j]=inewsize;

        } /* end if */

    } /* end j */

    return;
}


/********** PUT AWAY LATER ************************************/
/*                                             */
/* rlater has a similar structure to basicq     */
/* ilatersize is the physical memory and       */
/* ilater is the number of tracers             */

static void put_away_later(struct All_variables *E, int j, int it)
{
    int kk;
    void expand_later_array();


    /* The first tracer in initiates memory allocation. */
    /* Memory is freed after parallel communications    */

    if (E->trace.ilatersize[j]==0) {

        E->trace.ilatersize[j]=E->trace.max_ntracers[j]/5;

        for (kk=0;kk<=((E->trace.number_of_tracer_quantities)-1);kk++) {
            if ((E->trace.rlater[j][kk]=(double *)malloc(E->trace.ilatersize[j]*sizeof(double)))==NULL) {
                fprintf(E->trace.fpt,"AKM(put_away_later)-no memory (%d)\n",kk);
                fflush(E->trace.fpt);
                exit(10);
            }
        }
    } /* end first particle initiating memory allocation */


    /* Put tracer in later array */

    E->trace.ilater[j]++;

    if (E->trace.ilater[j] >= (E->trace.ilatersize[j]-5)) expand_later_array(E,j);

    /* stack basic and extra quantities together (basic first) */

    for (kk=0;kk<=((E->trace.number_of_basic_quantities)-1);kk++)
        E->trace.rlater[j][kk][E->trace.ilater[j]]=E->trace.basicq[j][kk][it];

    for (kk=0;kk<=((E->trace.number_of_extra_quantities)-1);kk++)
        E->trace.rlater[j][E->trace.number_of_basic_quantities+kk][E->trace.ilater[j]]=E->trace.extraq[j][kk][it];


    return;
}


/****** EXPAND LATER ARRAY *****************************************/

void expand_later_array(struct All_variables *E, int j)
{

    int inewsize;
    int kk;
    int icushion;

    /* expand rlater by 20% */

    icushion=100;

    inewsize=E->trace.ilatersize[j]+E->trace.ilatersize[j]/5+icushion;

    for (kk=0;kk<=((E->trace.number_of_tracer_quantities)-1);kk++) {
        if ((E->trace.rlater[j][kk]=(double *)realloc(E->trace.rlater[j][kk],inewsize*sizeof(double)))==NULL) {
            fprintf(E->trace.fpt,"AKM(expand later array )-no memory (%d)\n",kk);
            fflush(E->trace.fpt);
            exit(10);
        }
    }


    fprintf(E->trace.fpt,"Expanding physical memory of rlater to %d from %d\n",
            inewsize,E->trace.ilatersize[j]);

    E->trace.ilatersize[j]=inewsize;

    return;
}


/***** EJECT TRACER ************************************************/

static void eject_tracer(struct All_variables *E, int j, int it)
{

    int ilast_tracer;
    int kk;


    ilast_tracer=E->trace.ntracers[j];

    /* put last tracer in ejected tracer position */

    E->trace.ielement[j][it]=E->trace.ielement[j][ilast_tracer];

    for (kk=0;kk<=((E->trace.number_of_basic_quantities)-1);kk++)
        E->trace.basicq[j][kk][it]=E->trace.basicq[j][kk][ilast_tracer];

    for (kk=0;kk<=((E->trace.number_of_extra_quantities)-1);kk++)
        E->trace.extraq[j][kk][it]=E->trace.extraq[j][kk][ilast_tracer];



    E->trace.ntracers[j]--;

    return;
}



/********** ICHECK PROCESSOR SHELL *************/
/* returns -99 if rad is below current shell  */
/* returns 0 if rad is above current shell    */
/* returns 1 if rad is within current shell   */
/*                                            */
/* Shell, here, refers to processor shell     */
/*                                            */
/* shell is defined as bottom boundary up to  */
/* and not including the top boundary unless  */
/* the shell in question is the top shell     */

int icheck_processor_shell(struct All_variables *E,
                           int j, double rad)
{

    const int noz = E->lmesh.noz;
    const int nprocz = E->parallel.nprocz;
    double top_r, bottom_r;

    if (nprocz==1) return 1;

    top_r = E->sx[j][3][noz];
    bottom_r = E->sx[j][3][1];

    /* First check bottom */

    if (rad<bottom_r) return -99;


    /* Check top */

    if (rad<top_r) return 1;

    /* top processor */

    if ( (rad<=top_r) && (E->parallel.me_loc[3]==nprocz-1) ) return 1;

    /* If here, means point is above processor */
    return 0;
}


/********* ICHECK THAT PROCESSOR SHELL ********/
/*                                            */
/* Checks whether a given radius is within    */
/* a given processors radial domain.          */
/* Returns 0 if not, 1 if so.                 */
/* The domain is defined as including the bottom */
/* radius, but excluding the top radius unless   */
/* we the processor domain is the one that       */
/* is at the surface (then both boundaries are   */
/* included).                                    */

int icheck_that_processor_shell(struct All_variables *E,
                                int j, int nprocessor, double rad)
{
    int icheck_processor_shell();
    int me = E->parallel.me;

    /* nprocessor is right on top of me */
    if (nprocessor == me+1) {
        if (icheck_processor_shell(E, j, rad) == 0) return 1;
        else return 0;
    }

    /* nprocessor is right on bottom of me */
    if (nprocessor == me-1) {
        if (icheck_processor_shell(E, j, rad) == -99) return 1;
        else return 0;
    }

    /* Shouldn't be here */
    fprintf(E->trace.fpt, "Should not be here\n");
    fprintf(E->trace.fpt, "Error(check_shell) nprocessor: %d, radius: %f\n",
            nprocessor, rad);
    fflush(E->trace.fpt);
    exit(10);

    return 0;
}
/*lhy 20181011
 * read in composition from a file*/
void read_chemicals_from_files(struct All_variables *E, float **comp_e, int *nf_comp)
{
	const int fncomp = 4;	
	const int ncomp = E->composition.ncomp;	
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];
	int i,j,k,el,n,m,ll,mm,fj;
  	float v1, v2, v3, T, g, tt, volume1, integral1;
	float **comp_d;
	float *C;
  	char output_file[255], input_s[1000];
  	FILE *fp;
	int Is0 = (E->parallel.me==0);
/*read in composition from file*/
	if(Is0){
		fprintf(stderr,"read composition from file\n");
	}
	sprintf(output_file,"%s.coord.%d.inp.dat",E->control.old_P_file,E->parallel.me);
	comp_d = (float**)malloc((ncomp+1)*sizeof(float*));
	for(i=0;i<=ncomp;i++){
		comp_d[i] = (float*)malloc((E->lmesh.nno+1)*sizeof(float));
	}
  	C = (float*)malloc(fncomp*sizeof(float));
  	fp = fopen(output_file,"r");
  	if (fp == NULL) {
    	fprintf(stderr,"(Tracer_setup.c #1) Cannot open %s\n",output_file);
	    parallel_process_termination();
  	}
  	if (E->parallel.me==0){
    	fprintf(stderr,"Reading %s for initial composition\n",output_file);
	}
    for(i=1;i<=E->lmesh.nno;i++)  {
    	fgets(input_s,1000,fp);
      	if(sscanf(input_s,"%g %g %f %f %f %f %f %f %f",&(v1),&(v1),&(T),
			&(T),&(T),&(C[0]),&(C[1]),&(C[2]),&(C[3])) != 9) {
        	fprintf(stderr,"Error while reading file '%s'\n", output_file);
        	exit(8);
      	}
		for(j=0;j<=ncomp;j++){
			fj = nf_comp[j];
			comp_d[j][i] = C[fj];
		}
    }
	if(E->parallel.me==0){
	for(i=1;i<=10;i++){
			fprintf(stderr,"comp_d[%d] : %.4e, %.4e\n",i,comp_d[0][i],comp_d[1][i]); //debug
		}
	}
  	fclose (fp);
/*get element value from node value*/
	if(Is0){
		fprintf(stderr,"get element compostion value from node value\n");
	}
	for(el=1;el<=E->lmesh.nel;el++)
	for(k=0;k<=ncomp;k++){
		volume1 = 0.0;
		integral1 = 0.0;
		for(j=1;j<=vpts;j++)
		for(i=1;i<=ends;i++) {
			n = E->ien[1][el].node[i];
			volume1 += E->N.vpt[GNVINDEX(i,j)] * E->gDA[1][el].vpt[j];
			integral1 += comp_d[k][n] * E->N.vpt[GNVINDEX(i,j)] * E->gDA[1][el].vpt[j];
        }
		comp_e[k][el] = integral1/volume1;
	}
  
	for(i=0;i<=ncomp;i++){
		free(comp_d[i]);
	}
	free(comp_d);
  	free(C);
}


void calculate_melting_flux(E)
   struct All_variables *E;
{

const float year=365*24*3600;
const int nflavors=E->trace.nflavors;
int Is0=(E->parallel.me==0);
int number_of_tracers;
int itracer,flavor;
int ipos,intpos;
int ispecie;
int crust_specie;
double z0,x,z,theta,phi,rad;
float Dcrust;
float *aveT,T300;
int i,j,k,kk,Tnode,node,iel,m,mm,nz,comp;
static int been_here=0;
int el,el_above;
float R=E->data.radius/1e3;
float latent_heat;

float P,T,depth;
double Tnd;
float degree_of_dry_melting(float T, float P, float Mcpx, float *Ts);
float katz_wet_melting(float T, float P, float Mcpx, float water);
float kelley_wetmelting(float T, float P, float water, float Mcpx, int fertility);
float depth_pressure();
double get_temperature(struct All_variables *E,
                       int j, int nelem,
                       double theta, double phi, double rad);
double get_element_area();
float F_tracer,F0_tracer,F1_tracer,Fout_tracer;
float Tsol,Fout,F_inf,*velocity;
float **Fout_el_tracer_flavored,**dFdt_el_tracer_flavored,**MFsurf_tracer_flavored;
float *Fout_el_tracer,*Fout_el_tracer1,*CO2_el_tracer,*dFdt_el_tracer,*dCO2dt_el_tracer,*dFdt_el_euler;
float *MFtotal_tracer_flavored, *MFglobal_tracer_flavored;
float *MFsurf_tracer,*CO2_tracer,*MFsurf_euler,MFtotal_tracer,MFtotal_euler,CO2_total,CO2_global;
double *Ts,*Tl;
int *NT;
int N;
int ielx,iely,ielz,ielxy;
float dt,dz,dr;
float velocity_scale;
float time_scale,Tave;
float *MFave,*MFmax,MFave_total,MFmax_total;
float Fmax,Fave,velo_deepest;
float volume;
double area;
float dx;
float MFglobal_tracer,MFglobal_euler,MFglobal_Fave,MFglobal_Fmax;
float rad1,rad2;
float U_dot_grad_F();
float F[9];
double **comp_melt;

FILE *fp;
char tempstring[256];

parallel_process_sync(E);
/*parameters and allocate space*/
latent_heat=E->lunar.latent_heat;
#ifdef MELT_DEBUG
if(E->parallel.me==0) //debug
	fprintf(stderr,"calculating melt flux via Dr LiMingming's implement\n");
#endif

m=E->sphere.caps_per_proc;

if(m!=1)
{
	fprintf(stderr,"caps per proc must be 1\n");
	parallel_process_termination();
}

time_scale=pow(E->data.radius,2.0)/E->data.therm_diff/year; //radius^2/kappa in year
velocity_scale=E->data.therm_diff/E->data.radius*year; //kappa/radius in m/year
#ifdef MELT_DEBUG
if(Is0){
	fprintf(stderr,"time_scale=%.4e, velocity_scale=%.4e\n",time_scale,velocity_scale);
} //lhy debug
#endif
dt=time_scale*E->advection.dt_before;
F_inf=E->control.F_inf;
if(E->control.melting_model==4){
	Ts=E->sol[1];
	Tl=E->liq[1];
}
Fout_el_tracer_flavored=(float**)malloc((nflavors)*sizeof(float*));
for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
	Fout_el_tracer_flavored[flavor]=(float*)malloc((E->lmesh.nel+1)*sizeof(float)); }
Fout_el_tracer=malloc((E->lmesh.nel+1)*sizeof(float));
Fout_el_tracer1=malloc((E->lmesh.nel+1)*sizeof(float));
CO2_el_tracer=malloc((E->lmesh.nel+1)*sizeof(float));
dFdt_el_tracer_flavored=(float**)malloc((nflavors)*sizeof(float*));
for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
	dFdt_el_tracer_flavored[flavor]=(float*)malloc((E->lmesh.nel+1)*sizeof(float)); 
}
dFdt_el_tracer=malloc((E->lmesh.nel+1)*sizeof(float));
dCO2dt_el_tracer=malloc((E->lmesh.nel+1)*sizeof(float));
dFdt_el_euler=malloc((E->lmesh.nel+1)*sizeof(float));
velocity=malloc((E->lmesh.nel+1)*sizeof(float));
NT=malloc((E->lmesh.nel+1)*sizeof(int));
MFsurf_tracer_flavored=(float**)malloc((nflavors)*sizeof(float*));
for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
	MFsurf_tracer_flavored[flavor]=(float*)malloc((E->lmesh.elx*E->lmesh.ely+1)*sizeof(float)); 
}
MFsurf_tracer=malloc((E->lmesh.elx*E->lmesh.ely+1)*sizeof(float));
CO2_tracer=malloc((E->lmesh.elx*E->lmesh.ely+1)*sizeof(float));
MFsurf_euler=malloc((E->lmesh.elx*E->lmesh.ely+1)*sizeof(float));
MFave=malloc((E->lmesh.elx*E->lmesh.ely+1)*sizeof(float));
MFmax=malloc((E->lmesh.elx*E->lmesh.ely+1)*sizeof(float));
MFtotal_tracer_flavored=(float*)malloc((nflavors)*sizeof(float));
MFglobal_tracer_flavored=(float*)malloc((nflavors)*sizeof(float));

for(iel=1;iel<=E->lmesh.nel;iel++)
{
	for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
		Fout_el_tracer_flavored[flavor][iel]=0.0;
	}
	Fout_el_tracer[iel]=0.0;
	Fout_el_tracer1[iel]=0.0;
	CO2_el_tracer[iel]=0.0;
	for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
		dFdt_el_tracer_flavored[flavor][iel]=0.0;
	}
	dFdt_el_tracer[iel]=0.0;
	dCO2dt_el_tracer[iel]=0.0;
	dFdt_el_euler[iel]=0.0;
	NT[iel]=0;
	velocity[iel]=0.0;
	for(i=1;i<=8;i++)
	{
		velocity[iel]+=E->sphere.cap[m].V[3][E->ien[m][iel].node[i]]*(float)E->N.ppt[GNPINDEX(i,1)];
	}
	velocity[iel] = velocity[iel]*velocity_scale;
}
for(i=1;i<=E->lmesh.elx*E->lmesh.ely;i++)
{
	for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
		MFsurf_tracer_flavored[flavor][i]=0.0; 
	}
	MFsurf_tracer[i]=MFsurf_euler[i]=0.0;
	CO2_tracer[i]=0.0;
	MFave[i]=MFmax[i]=0.0;
}

number_of_tracers=E->trace.ntracers[m];

comp_melt=(double**)malloc((E->sphere.caps_per_proc+1)*sizeof(double*));
for(mm=0;mm<=E->sphere.caps_per_proc;mm++) {
	comp_melt[mm]=(double*)malloc((E->lmesh.nno+1)*sizeof(double));
}

#ifdef MELT_DEBUG
if(Is0){
	fprintf(stderr,"calculate_melting_flux local initialization done\n"); //lhy debug
}
#endif


/*initialization*/
if(been_here==0)
{
	if(E->control.restart==0 && E->trace.ic_method!=2)
	{
		for(itracer=1;itracer<=number_of_tracers;itracer++)
		{
			//ip_F is 2; ip_Fout is 3; fco2 is 4
			E->trace.extraq[m][2][itracer]=0.0;
		}
	}

	
	if(E->control.restart==0)
	{
		for(itracer=1;itracer<=number_of_tracers;itracer++)
			E->trace.extraq[m][4][itracer]=E->control.fco2;
	}

	E->F0_el=malloc((E->lmesh.nel+1)*sizeof(float));
	E->F_el=malloc((E->lmesh.nel+1)*sizeof(float));
	E->UdotgradF=malloc((E->lmesh.nel+1)*sizeof(float));
	for(i=1;i<=E->lmesh.nel;i++)
	{
		E->UdotgradF[i]=0.0;
		E->F0_el[i]=E->F_el[i]=0.0;
		E->melt_latent_heat[i]=0.0;
	}

	E->Fnode=malloc((E->lmesh.nno+1)*sizeof(float));
	for(i=1;i<=E->lmesh.nno;i++)
	{
		E->Fnode[i]=0.0;
	}


}

#ifdef MELT_DEBUG
if(Is0){
	fprintf(stderr,"calculate_melting_flux initialization done\n"); //lhy debug
}
#endif

/*get melting tracerwise*/
for(itracer=1;itracer<=number_of_tracers;itracer++)
{
	el=E->trace.ielement[m][itracer];
	rad=0.0;

	for(node=1;node<=8;node++)
	{
		rad+=E->sx[m][3][E->ien[m][el].node[node]]*E->N.ppt[GNPINDEX(node,1)];
	}
	rad=E->trace.basicq[m][2][itracer];

	theta=E->trace.basicq[m][0][itracer];
	phi=E->trace.basicq[m][1][itracer];

	F0_tracer=E->trace.extraq[m][2][itracer];

	depth=R*(1.0-rad);
	Tnd=get_temperature(E,m,el,theta,phi,rad);

	T=E->data.ref_temperature*Tnd+E->data.Tsurf; //I think here should add on temperature at the top

	//T=1873.0;P=8.0;
	if(E->control.melting_model==1)
	{
		P=depth_pressure(rad);
		F_tracer=degree_of_dry_melting(T,P,E->control.Mcpx,&Tsol);
	}
	else if(E->control.melting_model==2)
	{
		P=depth_pressure(rad);
		F_tracer=kelley_wetmelting(T,P,E->control.water, E->control.Mcpx,E->control.fertility);
	}
	else if(E->control.melting_model==3)
	{
		P=depth_pressure(rad);
		F_tracer=katz_wet_melting(T,P,E->control.Mcpx,E->control.water);
	}
	else if(E->control.melting_model==4)
	{
		F_tracer=linear_melting(E,Tnd,rad);
	}
	else
	{
		fprintf(stderr,"Wrong melting model\n");
		exit(10);
	}
	if(F_tracer>F0_tracer && F_tracer>F_inf)
	{
		E->trace.extraq[m][2][itracer]=F_tracer;
		E->trace.extraq[m][3][itracer]=F_tracer-F0_tracer; //lhy test
		E->trace.extraq[m][5][itracer]=F_tracer-F0_tracer; //lhy tes
	}
	else if(F_tracer>F0_tracer)
	{
		E->trace.extraq[m][2][itracer]=F_tracer;
		E->trace.extraq[m][5][itracer]=F_tracer-F0_tracer; //lhy tes
	}
	else
	{
		E->trace.extraq[m][2][itracer]=F0_tracer;
		E->trace.extraq[m][3][itracer]=0.0;
		E->trace.extraq[m][5][itracer]=0.0; //lhy test
	}
	NT[el]=NT[el]+1;
	flavor=E->trace.extraq[m][0][itracer];
	Fout_el_tracer_flavored[flavor][el]=Fout_el_tracer_flavored[flavor][el]+E->trace.extraq[m][3][itracer];
	Fout_el_tracer[el]=Fout_el_tracer[el]+E->trace.extraq[m][3][itracer];
	Fout_el_tracer1[el]=Fout_el_tracer1[el]+E->trace.extraq[m][5][itracer];
	if(F_tracer>F0_tracer && F_tracer>F_inf)
	{
		CO2_el_tracer[el]=CO2_el_tracer[el]+E->trace.extraq[m][4][itracer]*E->eco[m][el].area*R*R*R*3300.0*1e9; //Here I cannot figure out how 3300 comes to be. But it isn't important after all
		E->trace.extraq[m][4][itracer]=0.0;
	}
}

#ifdef MELT_DEBUG
if(Is0){
	fprintf(stderr,"calculate_melting_flux melt by tracer done\n"); //lhy debug
}
#endif

for(i=1;i<=E->lmesh.nel;i++)
{
	E->melt_latent_heat[i]=0.0;
}
/*take average in every element*/
for(iely=1;iely<=E->lmesh.ely;iely++)
for(ielx=1;ielx<=E->lmesh.elx;ielx++)
for(ielz=1;ielz<=E->lmesh.elz;ielz++)
{
	el=ielz+(ielx-1)*E->lmesh.elz+(iely-1)*E->lmesh.elx*E->lmesh.elz;
	if(NT[el]==0)
	{
		for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
			Fout_el_tracer_flavored[flavor][el]=0.0;
		}
	Fout_el_tracer[el]=0.0;
	Fout_el_tracer1[el]=0.0;
		CO2_el_tracer[el]=0.0; }
	else
	{
		for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
			Fout_el_tracer_flavored[flavor][el]=Fout_el_tracer_flavored[flavor][el]/NT[el];
		}
		Fout_el_tracer[el]=Fout_el_tracer[el]/NT[el];
		Fout_el_tracer1[el]=Fout_el_tracer1[el]/NT[el];
		CO2_el_tracer[el]=CO2_el_tracer[el]/NT[el];
	}
//	if(E->parallel.me==0)
//	printf("el=%d %.6e\n",el,Fout_el_tracer[el]);
	if(been_here==0)
	{
		for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
			dFdt_el_tracer_flavored[flavor][el]=0.0;
		}
		dFdt_el_tracer[el]=0.0;
		dCO2dt_el_tracer[el]=0.0;
	}
	else
	{	
		for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
			dFdt_el_tracer_flavored[flavor][el]=Fout_el_tracer_flavored[flavor][el]/dt;
		}
		dFdt_el_tracer[el]=Fout_el_tracer[el]/dt;
		dCO2dt_el_tracer[el]=CO2_el_tracer[el]/dt;
	}
	//E->melt_latent_heat[el]=latent_heat*dFdt_el_tracer[el]*time_scale/(E->data.Cp*E->control.refT);
}

#ifdef MELT_DEBUG
if(Is0){
	fprintf(stderr,"calculate_melting_flux calculate melt average in element\n"); //lhy debug
}
#endif

/*calculate melt fraction in node*/
for(j=1;j<=E->lmesh.noy;j++)
for(i=1;i<=E->lmesh.nox;i++)
for(k=1;k<=E->lmesh.noz;k++)
{
	node=k+(i-1)*E->lmesh.noz+(j-1)*E->lmesh.noz*E->lmesh.nox;

	rad=E->sx[m][3][node];
	/*
	if(rad<0.9529)
	{
		continue;
	}
	*/
	Tnd=E->T[m][node];
	depth=R*(1.0-rad);
	T=E->data.ref_temperature*Tnd+E->data.Tsurf; //I think here should add on temperature at the top

	//T=1873.0;P=8.0;
	if(E->control.melting_model==1)
	{
		P=depth_pressure(rad);
		E->Fnode[node]=degree_of_dry_melting(T,P,E->control.Mcpx,&Tsol);
	}
	else if(E->control.melting_model==2)
	{
		P=depth_pressure(rad);
		E->Fnode[node]=kelley_wetmelting(T,P,E->control.water, E->control.Mcpx,E->control.fertility);
	}
	else if(E->control.melting_model==3)
	{
		P=depth_pressure(rad);
		E->Fnode[node]=katz_wet_melting(T,P,E->control.Mcpx,E->control.water);
	}
	else if(E->control.melting_model==4)
	{
		nz = k+E->parallel.me_loc[3]*(E->lmesh.noz-1);
		E->Fnode[node] = max((Tnd-Ts[nz])/(Tl[nz]-Ts[nz]),0.0);
#ifdef MELT_DEBUG
		if(Is0&&i==1&&j==1){
			fprintf(stderr,"lnz=%d, Ts=%.4e, Tl=%.4e\n",k,Ts[nz],Tl[nz]);
		}//lhy debug
#endif
	}
	else
	{
		fprintf(stderr,"Wrong melting model\n");
		parallel_process_termination();
	}
	E->melting[1][node]=E->Fnode[node];
}
E->melting_volume = return_bulk_value_d(E,E->melting,0); //calculate melt volume here
for(comp=0;comp<E->composition.ncomp;comp++){
for(node=1;node<=E->lmesh.nno;node++){
	comp_melt[1][node] = E->melting[1][node]*E->composition.comp_node[1][comp][node];
}
E->comp_melting_volume[comp] = return_bulk_value_d(E,comp_melt,0);
}



for(iely=1;iely<=E->lmesh.ely;iely++)
for(ielx=1;ielx<=E->lmesh.elx;ielx++)
for(ielz=1;ielz<=E->lmesh.elz;ielz++)
{
	el=ielz+(ielx-1)*E->lmesh.elz+(iely-1)*E->lmesh.elx*E->lmesh.elz;
	E->F0_el[el]=max(E->F_el[el],F_inf);
	E->F_el[el]=0.0;
	for(i=1;i<=8;i++)
	{
		E->F_el[el]+=E->Fnode[E->ien[m][el].node[i]]*0.125;
	}
	if(E->F_el[el]>F_inf)
	{
		if(been_here==0)
		{
			dFdt_el_euler[el]=E->UdotgradF[el];
		}
		else
		{
			dFdt_el_euler[el]=(E->F_el[el]-E->F0_el[el])/dt+E->UdotgradF[el];
		}
/*
		if(E->parallel.me==0)
		{
			printf("%.6e %.6e\n",(E->F_el[el]-E->F0_el[el])/dt,fabs(E->UdotgradF[el]));
		}
*/
	}

	for(i=1;i<=8;i++)
	{
		F[i]=E->Fnode[E->ien[m][el].node[i]];
	}
	E->UdotgradF[el]=U_dot_grad_F(E,el,F);
}

for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
	MFtotal_tracer_flavored[flavor]=0.0;
}
MFtotal_tracer=MFtotal_euler=0.0;
MFave_total=MFmax_total=0.0;
CO2_total=0.0;


i=0;
for(iely=1;iely<=E->lmesh.ely;iely++)
for(ielx=1;ielx<=E->lmesh.elx;ielx++)
{
	i++;
	Fmax=Fave=0.0;
	N=0;
	for(ielz=1;ielz<=E->lmesh.elz;ielz++)
	{
		el=ielz+(ielx-1)*E->lmesh.elz+(iely-1)*E->lmesh.elx*E->lmesh.elz;
		volume=E->eco[m][el].area*R*R*R;
		if(dFdt_el_euler[el]>0.0)
		{
			MFsurf_euler[i]=MFsurf_euler[i]+dFdt_el_euler[el]*volume;
			MFtotal_euler=MFtotal_euler+dFdt_el_euler[el]*volume;
		}
		for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
			MFsurf_tracer_flavored[flavor][i]=MFsurf_tracer_flavored[flavor][i]+dFdt_el_tracer_flavored[flavor][el]*volume; 
		}
		MFsurf_tracer[i]=MFsurf_tracer[i]+dFdt_el_tracer[el]*volume;
		if(Is0&&i==1){
			fprintf(stderr,"el = %d, dFdt = %.4e, volume = %.4e\n",el, dFdt_el_tracer[el], volume);
		}
		for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
			MFtotal_tracer_flavored[flavor]=MFtotal_tracer_flavored[flavor]+dFdt_el_tracer_flavored[flavor][el]*volume;
		}
		MFtotal_tracer=MFtotal_tracer+dFdt_el_tracer[el]*volume;
		CO2_tracer[i]=CO2_tracer[i]+dCO2dt_el_tracer[el];
		CO2_total+=dCO2dt_el_tracer[el];

		if(E->F_el[el]>F_inf)
		{
			N++;
			Fave=Fave+E->F_el[el];
			if(E->F_el[el]>Fmax)
			{
				Fmax=E->F_el[el];
			}
		}
	}
	if(N>0)
	{
		Fave=Fave/N;
	}
	else
	{
		Fave=0.0;
	}

	velo_deepest=0.0;
	area=0.0;
	for(ielz=1;ielz<=E->lmesh.elz;ielz++)
	{
		el=ielz+(ielx-1)*E->lmesh.elz+(iely-1)*E->lmesh.elx*E->lmesh.elz;
		if(E->F_el[el]>F_inf && velocity[el]>0.0)
		{
			velo_deepest=velocity[el];
			area=get_element_area(E,m,el);
			break;
		}
	}
	MFave[i]=velo_deepest*Fave*area*R*R;
	MFmax[i]=velo_deepest*Fmax*area*R*R;
	MFave_total=MFave_total+MFave[i];
	MFmax_total=MFmax_total+MFmax[i];
}

if(E->control.melt_latent_heating==1)
{
	for(el=1;el<=E->lmesh.nel;el++)
	{
		if(dFdt_el_euler[el]>0.0)
		{
			E->melt_latent_heat[el]=latent_heat*dFdt_el_euler[el]*time_scale/(E->data.Cp*E->data.ref_temperature);
		}
	}
}
else if(E->control.melt_latent_heating==2)
{
	if(Is0){
		for(el=1;el<E->lmesh.elz;el++)
		{
			fprintf(stderr,"el = %d, melt_rate = %.4e, melt_rate1 = %.4e\n",el,Fout_el_tracer[el]/dt,Fout_el_tracer1[el]/dt);//lhy debug
		}
	}
	for(el=1;el<=E->lmesh.nel;el++)
	{
		E->melt_latent_heat[el]=latent_heat*Fout_el_tracer1[el]/dt*time_scale/(E->data.Cp*E->data.ref_temperature);
	}
}


MPI_Allreduce(&MFave_total, &MFglobal_Fave, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(&MFmax_total, &MFglobal_Fmax, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
	MPI_Allreduce(MFtotal_tracer_flavored+flavor, MFglobal_tracer_flavored+flavor, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
}
MPI_Allreduce(&MFtotal_tracer, &MFglobal_tracer, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(&MFtotal_euler, &MFglobal_euler, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(&CO2_total, &CO2_global, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

//printf("%d %f %f %f %f\n",E->parallel.me,MFave_total,MFmax_total,MFtotal_tracer,MFtotal_euler);

//fprintf(stderr,"%07d %.4e %.4e %.4e %.4e %.4e\n",E->advection.timesteps,E->monitor.elapsed_time,MFglobal_tracer,MFglobal_euler,MFglobal_Fave,MFglobal_Fmax);

if(E->parallel.me==0)
{
	fprintf(stderr,"%07d %.4e %.4e %.4e %.4e %.4e",E->advection.timesteps,dt,E->monitor.elapsed_time*time_scale,MFglobal_tracer,MFglobal_euler,CO2_global);
	for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
		fprintf(stderr," %.4e",MFglobal_tracer_flavored[flavor]);
	}
	fprintf(stderr,"\n");
	fprintf(E->fp_mf,"%07d %.4e %.4e %.4e %.4e %.4e",E->advection.timesteps,dt,E->monitor.elapsed_time*time_scale,MFglobal_tracer,MFglobal_euler,CO2_global);
	for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
		fprintf(E->fp_mf," %.4e",MFglobal_tracer_flavored[flavor]);
	}
	fprintf(E->fp_mf,"\n");
	fflush(E->fp_mf);
}


//if(E->monitor.solution_cycles % E->control.record_every == 0)
if(E->monitor.solution_cycles % E->control.MF_save_step == 0)
{
	sprintf(tempstring,"%s.MF.%d.%d",E->control.data_file,E->parallel.me,E->monitor.solution_cycles);
	if(access(tempstring,F_OK)==-1)
	{
		fp=fopen(tempstring,"w");
		for(i=1;i<=E->lmesh.elx*E->lmesh.ely;i++)
		{
			fprintf(fp,"%.6e %.6e",MFsurf_tracer[i],MFsurf_euler[i]);
			for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
				fprintf(fp," %.6e",MFsurf_tracer_flavored[flavor][i]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
	}
	else
	{
		fprintf(stderr,"%s already exist\n",tempstring);
	}

/*
	sprintf(tempstring,"%s/CO2/CO2.%d.%d",E->control.data_dir,E->parallel.me,E->monitor.solution_cycles);
	if(access(tempstring,F_OK)==-1)
	{
		fp=fopen(tempstring,"w");
		for(i=1;i<=E->lmesh.elx*E->lmesh.ely;i++)
		{
			fprintf(fp,"%.6e\n",CO2_tracer[i]);
		}
		fclose(fp);
	}
	else
	{
		fprintf(stderr,"%s already exist\n",tempstring);
	}
*/


/*	sprintf(tempstring,"%s.dFdt.%d.%d",E->control.data_file,E->parallel.me,E->monitor.solution_cycles);
	if(access(tempstring,F_OK)==-1)
	{
		fp=fopen(tempstring,"w");
		for(iely=1;iely<=E->lmesh.ely;iely++)
		for(ielx=1;ielx<=E->lmesh.elx;ielx++)
		for(ielz=1;ielz<=E->lmesh.elz;ielz++)
		{
			el=ielz+(ielx-1)*E->lmesh.elz+(iely-1)*E->lmesh.elx*E->lmesh.elz;
			if(dFdt_el_euler[el]>0.0 || dFdt_el_tracer[el]>0.0)
			{
//				fprintf(fp,"%d %.6e\n",el,dFdt_el_euler[el]);
				fprintf(fp,"%d %.6e %.6e",el,dFdt_el_tracer[el],dFdt_el_euler[el]);
				for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
					fprintf(fp," %.6e",dFdt_el_tracer_flavored[flavor][el]);
				}
				fprintf(fp,"\n");
			}
		}
		fclose(fp);
	}
	else
	{
		fprintf(stderr,"%s already exist\n",tempstring);
	}
	*/
}

for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
	free(Fout_el_tracer_flavored[flavor]); 
}
free(Fout_el_tracer_flavored);
Fout_el_tracer=malloc((E->lmesh.nel+1)*sizeof(float));
for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
	free(dFdt_el_tracer_flavored[flavor]); 
}
free(dFdt_el_tracer_flavored);
for (flavor=0; flavor<nflavors; flavor++) { //modified by lhy
	free(MFsurf_tracer_flavored[flavor]); 
}
free(MFsurf_tracer_flavored);
free(MFtotal_tracer_flavored);
free(MFglobal_tracer_flavored);
for(mm=0;mm<=E->sphere.caps_per_proc;mm++) {
	free(comp_melt[mm]);
}
free(comp_melt);
free(Fout_el_tracer);
free(CO2_el_tracer);
free(dCO2dt_el_tracer);
free(NT);
free(dFdt_el_tracer);
free(dFdt_el_euler);
free(velocity);
free(MFsurf_tracer);
free(CO2_tracer);
free(MFsurf_euler);
free(MFave);
free(MFmax);

been_here++;
return;
}

void calculate_melting_flux_postp(E)
   struct All_variables *E;
{
	if (E->parallel.me == 0)
		fprintf(stderr, "calculate_melting_flux_postp\n");
	float *MFsurf;	// Initialization
	const float year = 365 * 24 * 3600;
	int i, j, k;
	int el, ielx, iely, ielz, node, nz;
	float volume, R, MFTotal, MFGlobal, temp, time_scale, Tnd;
	double *Ts, *Tl;
	static int been_here = 0;
	char filename[500];
	FILE *fp;
	R = E->data.radius/1e3;
	time_scale = pow(E->data.radius, 2.0) / E->data.therm_diff / year; //radius^2/kappa in year
	MFsurf = (float*)malloc((E->lmesh.elx * E->lmesh.ely + 1) * sizeof(float));
	for(i = 0; i <= E->lmesh.elx * E->lmesh.ely; i++)
		MFsurf[i] = 0.0;
	MFTotal = 0.0;
	if(been_here == 0){
	    E->Fnode = (float*)malloc((E->lmesh.nno+1)*sizeof(float));
	    for(i = 0; i <= E->lmesh.nno; i++)
		    E->Fnode[i]=0.0;
	    E->F_el = (float*)malloc((E->lmesh.nel+1)*sizeof(float));
	    for(i = 0; i <= E->lmesh.nel; i++)
	    	E->F_el[i] = 0.0;
    }
	Ts = E->sol[1];
	Tl = E->liq[1];
	for(j = 1; j <= E->lmesh.noy; j++) // calculate melt
	for(i = 1; i <= E->lmesh.nox; i++)
	for(k = 1; k <= E->lmesh.noz; k++)
	{
		node = k + (i - 1) * E->lmesh.noz +
			  (j - 1) * E->lmesh.noz * E->lmesh.nox;
		Tnd = E->T[1][node];
		if(E->control.melting_model == 4)
		{
			nz = k+E->parallel.me_loc[3]*(E->lmesh.noz-1);
			E->Fnode[node] = max((Tnd-Ts[nz])/(Tl[nz]-Ts[nz]),0.0);
	#ifdef MELT_DEBUG
			if(Is0&&i==1&&j==1){
				fprintf(stderr,"lnz=%d, Ts=%.4e, Tl=%.4e\n",k,Ts[nz],Tl[nz]);
			}//lhy debug
	#endif
		}
		else
		{
			fprintf(stderr,"Wrong melting model\n");
			parallel_process_termination();
		}
		E->melting[1][node] = E->Fnode[node];
	}
	for(el = 1; el <= E->lmesh.nel; el++)
	{
		for(i = 1; i <= 8; i++)
		{
			E->F_el[el] += E->Fnode[E->ien[1][el].node[i]] * 0.125;
		}
	}
	i=0;
	for(iely = 1; iely <= E->lmesh.ely; iely++)
	for(ielx = 1; ielx <= E->lmesh.elx; ielx++)
	{
		i++;
		for(ielz=1; ielz <= E->lmesh.elz; ielz++)
		{
			el = ielz + (ielx - 1) * E->lmesh.elz +
				(iely - 1) * E->lmesh.elx * E->lmesh.elz;
			volume = E->eco[1][el].area * R * R * R;
			temp = E->F_el[el] * volume;
			MFsurf[i] += temp;	// MFsurf in km^3
			MFTotal += temp;
		}
	}
	MPI_Allreduce(&MFTotal, &MFGlobal, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	if(E->parallel.me == 0){    // Output
		fprintf(E->fp_mf, "%07d %.4e %.4e",E->advection.timesteps,
				E->monitor.elapsed_time * time_scale, MFGlobal);
		fprintf(E->fp_mf,"\n");
		fflush(E->fp_mf);
	}
	if(E->monitor.solution_cycles % E->control.MF_save_step == 0) // Output
	{
		sprintf(filename, "%s.MFp.%d.%d", E->control.data_file,
				E->parallel.me, E->monitor.solution_cycles);
        fp = fopen(filename, "w");
        if (fp == NULL) {
            fprintf(stderr,"(Initial_temperature.c #1) Cannot open %s\n", filename);
            parallel_process_termination();
        }
		for(i = 1; i <= E->lmesh.elx * E->lmesh.ely; i++)
		{
			fprintf(fp, "%.6e", MFsurf[i]);
			fprintf(fp, "\n");
		}
		fclose(fp);
	}
	free(MFsurf);	// free space
	been_here++;
	return;
}


static float depth_pressure(depth_citcom)
        double depth_citcom;
{
float density[43];
static float depth[43],pressure[43];
float depth_earth;
float P;
int i;
char tempstring[200],tempstring1[200];
FILE *fpin_prem;
static int been_here=0;

if(fabs(depth_citcom-1.0)<1e-9)
{
        return 0.0;
}


        depth_earth=6371.0*(1.0-depth_citcom);

if(been_here==0)
{
        sprintf(tempstring1,"../Reg/prem.input");

        if ((fpin_prem=fopen(tempstring1,"r"))==NULL)
        {
        	fprintf(stderr,"ERROR-file does not existdfdfdfd.\n");
			parallel_process_termination();
        }

        i=1;
        while(fgets(tempstring,200,fpin_prem)!=NULL)
        {
        sscanf(tempstring,"%f %f %f",&depth[i],&density[i],&pressure[i]);
        i=i+1;
        }
        for(i=1;i<42;i++)
        {
        if(depth[i]<depth_earth && depth_earth<depth[i+1])
                {
                P=pressure[i]+((depth_earth-depth[i])/(depth[i+1]-depth[i]))*(pressure[i+1]-pressure[i]);
                break;
                }
        }
        fclose(fpin_prem);
	been_here++;
	return P;
}
else
{
        for(i=1;i<42;i++)
        {
   	     if(depth[i]<depth_earth && depth_earth<depth[i+1])
                {
                	P=pressure[i]+((depth_earth-depth[i])/(depth[i+1]-depth[i]))*(pressure[i+1]-pressure[i]);
	                break;
                }
        }
	return P;
}
}

static double get_temperature(struct All_variables *E,
                       int j, int nelem,
                       double theta, double phi, double rad)
{
    int iwedge;
    const int sphere_key = 0;

    double shape[9];
    double temperature,T[7];
    int i;

    void velo_from_element_d();

    full_get_shape_functions(E, shape, nelem, theta, phi, rad);

	iwedge=shape[0];

	if(iwedge==1)
	{
		T[1]=E->T[j][E->ien[j][nelem].node[1]];
		T[2]=E->T[j][E->ien[j][nelem].node[2]];
		T[3]=E->T[j][E->ien[j][nelem].node[3]];
		T[4]=E->T[j][E->ien[j][nelem].node[5]];
		T[5]=E->T[j][E->ien[j][nelem].node[6]];
		T[6]=E->T[j][E->ien[j][nelem].node[7]];
	}
	else if(iwedge==2)
	{
		T[1]=E->T[j][E->ien[j][nelem].node[1]];
		T[2]=E->T[j][E->ien[j][nelem].node[3]];
		T[3]=E->T[j][E->ien[j][nelem].node[4]];
		T[4]=E->T[j][E->ien[j][nelem].node[5]];
		T[5]=E->T[j][E->ien[j][nelem].node[7]];
		T[6]=E->T[j][E->ien[j][nelem].node[8]];
	}

	temperature=0.0;
	for(i=1;i<=6;i++)
	{
		temperature+=T[i]*shape[i];
	}

    return temperature;
}

static float degree_of_dry_melting(float T, float P, float Mcpx, float *Ts)
{
  const float T0 = 0.0;
  const float A1 = 1085.7;  // in C
  const float A2 =  132.9;  // in C GPa-1
  const float A3 =   -5.1;  // in C GPa-2
  const float B1 = 1475.0;  // in C
  const float B2 =   80.0;  // in C GPa-1
  const float B3 =   -3.2;  // in C GPa-2
  const float C1 = 1780.0;  // in C
  const float C2 =   45.0;  // in C GPa-1
  const float C3 =   -2.0;  // in C GPa-2
  const float r0 =   0.50;  // in cpx/melt
  const float r1 =   0.08;  // in cpx/melt/GPa
  const float beta1 = 1.50; // exponent
  const float beta2 = 1.50; // exponent
  static int been_here = 0;
  FILE *fp;
  int i;
  float pp,Tsol,Tlhl,Tliq,Tprime,Rcpx,Fcpxout,Tcpxout,F;
  if (been_here++ == 0) {
    fp = fopen("katz_fig1.dat","w");
    for (i=0;i<=160;i++) {
       pp = i/20.;
       Tsol = A1+pp*(A2+pp*A3);
       Tlhl = B1+pp*(B2+pp*B3);
       Tliq = C1+pp*(C2+pp*C3);
       fprintf(fp,"%.02f %.01f %.01f %.01f\n",pp,Tsol,Tlhl,Tliq);
    }
    fclose(fp);
  }

  Tsol = A1+P*(A2+P*A3) + T0; // in K
  Tlhl = B1+P*(B2+P*B3) + T0;
  Tliq = C1+P*(C2+P*C3) + T0;
  *Ts = Tsol;

  Tprime = (T-Tsol)/(Tlhl-Tsol);
  if (T < Tsol) return(0.0);
  if (T > Tliq) return(1.0);
  Rcpx = r0 + r1*P;
if(P>3.5)
{
        Rcpx=r0+r1*3.5;
}
  Fcpxout = Mcpx/Rcpx;

  F = pow(Tprime,beta1);

  if (F > Fcpxout) {
     Tcpxout = pow(Fcpxout,1.0/beta1)*(Tlhl-Tsol)+Tsol;
     Tprime = (T-Tcpxout)/(Tliq-Tcpxout);
     F = Fcpxout + (1.0-Fcpxout)*pow(Tprime,beta2);
  }

  return F;
}

static float katz_wet_melting(float T, float P, float Mcpx, float water)
{
const float T0 = 0.0;
const float A1 = 1085.7;  // in C
const float A2 =  132.9;  // in C GPa-1
const float A3 =   -5.1;  // in C GPa-2
const float B1 = 1475.0;  // in C
const float B2 =   80.0;  // in C GPa-1
const float B3 =   -3.2;  // in C GPa-2
const float C1 = 1780.0;  // in C
const float C2 =   45.0;  // in C GPa-1
const float C3 =   -2.0;  // in C GPa-2
const float r0 =   0.50;  // in cpx/melt
const float r1 =   0.08;  // in cpx/melt/GPa
const float beta1 = 1.50; // exponent
const float beta2 = 1.50; // exponent
const float  K=43.0;
const float  gama=0.75;
const float  Dh2o=0.01;
const float  zeta1=12.00;
const float  zeta2=1.00;
const float  ramda=0.60;
static int been_here = 0;
FILE *fp;
int i;
float pp,Tsol,Tlhl,Tliq,Tprime,Rcpx,Fcpxout,Tcpxout,F;
float Xh2o,Xh2o_bulk,Xh2o_sat;
float startF,endF,error;
float deltaT;
float F_cal;
float Fopx,Fopx_cal;
char ch;
Tsol = A1+P*(A2+P*A3); // in oC 
Tlhl = B1+P*(B2+P*B3);
Tliq = C1+P*(C2+P*C3);

Xh2o_bulk=water;

Rcpx = r0 + r1*P;
if(P>3.5)
{
        Rcpx=r0+r1*3.5;
}
Fcpxout = Mcpx/Rcpx;


startF=0.0;
endF=1.0;
F=0.0;
do
{
        Xh2o=Xh2o_bulk/(Dh2o+F*(1.0-Dh2o));
        Xh2o_sat=zeta1*pow(P,ramda)+zeta2*P;
        if(Xh2o>Xh2o_sat)
        {
                Xh2o=Xh2o_sat;
        }
        deltaT=K*pow(Xh2o,gama);
        F_cal=(T-(Tsol-deltaT))/(Tlhl-Tsol);
        if(F_cal<0.0)
        {
                F_cal=0.0;
        }
        else
        {
                F_cal=pow(F_cal,beta1);
        }
        error=F_cal-F;
        if(error>0)
        {
                startF=F;
        }
        else
        {
                endF=F;
        }
        F=0.5*(startF+endF);
}while(fabs(error)>1e-7 && fabs(startF-endF)>1e-7);

if (F > Fcpxout)
{
        startF=0.0;
        endF=1.0;
        F=0.0;
        do
        {
                Xh2o=Xh2o_bulk/(Dh2o+F*(1.0-Dh2o));
                Xh2o_sat=zeta1*pow(P,ramda)+zeta2*P;
                if(Xh2o>Xh2o_sat)
                {
                        Xh2o=Xh2o_sat;
                }
                deltaT=K*pow(Xh2o,gama);
                Tcpxout = pow(Fcpxout,1.0/beta1)*(Tlhl-Tsol)+Tsol-deltaT;
                Tprime = (T-Tcpxout)/(Tliq-Tcpxout);
                if(Tprime<0.0)
                {
                        F_cal=Fcpxout;
                }
                else
                {
                        F_cal = Fcpxout + (1.0-Fcpxout)*pow(Tprime,beta2);
                }
                error=F_cal-F;
                if(error>0)
                {
                        startF=F;
                }
                else
                {
                        endF=F;
                }
                F=0.5*(startF+endF);
        }while(fabs(error)>1e-7 && fabs(startF-endF)>1e-7);
}
return F;
}

static float kelley_wetmelting(float T, float P, float water, float Mcpx, int fertility)
{
        float a,b,c,x,y;
        float water_cal,endF,startF,F,error;
        char ch;
        float term1,term2,term3;
        float Fcpxout,Rcpx;
        float Xh2o,Xh2o_bulk,Xh2o_sat;
        float Tprime;
        float deltaT;
        float F_cal;
        const float T0 = 0.0;
        const float A1 = 1085.7;  // in C
        const float A2 =  132.9;  // in C GPa-1
        const float A3 =   -5.1;  // in C GPa-2
        const float B1 = 1475.0;  // in C
        const float B2 =   80.0;  // in C GPa-1
        const float B3 =   -3.2;  // in C GPa-2
        const float C1 = 1780.0;  // in C
        const float C2 =   45.0;  // in C GPa-1
        const float C3 =   -2.0;  // in C GPa-2
        const float r0 =   0.50;  // in cpx/melt
        const float r1 =   0.08;  // in cpx/melt/GPa
        const float beta1 = 1.50; // exponent
        const float beta2 = 1.50; // exponent
        const float  K=43.0;
        const float  gama=0.75;
        const float  Dh2o=0.01;
        const float  zeta1=12.00;
        const float  zeta2=1.00;
        const float  ramda=0.60;
        float Tsol,Tlhl,Tliq,Tcpxout;

        if(fertility==1)
        {
                a=-5.1404654;
                b = 132.899012;
                c = 1120.66061;
                x = -221.34;
                y = 536.86;
        }
        else
        {
                a = -5.1404654;
                b = 132.899012;
                c = 1159.66061;
                x = -136.88;
                y = 332.01;
        }
        startF=0.0;
        endF=1.0;
        F=0.0;
        do
        {
                term1=Dh2o*(1.0-F)+F;
                term2=(T-(a*P*P+b*P+c)-(x*log(P)+y)*F)/(-60.0);
                if(term2<0.0)
                {
                        startF=F;
                        F=0.5*(startF+endF);
                        error=1.0;
                }
                else
                {
                        term3=pow(term2,1.85);
                        water_cal=term1*term3;
                        error=water_cal-water;
                        if(error<0.0)
                        {
                                startF=F;
                        }
                        else
                        {
                                endF=F;
                        }
                }
                F=0.5*(startF+endF);
        }while(fabs(error)>1e-7 && fabs(startF-endF)>1e-7);

        Rcpx = r0 + r1*P;
        if(P>3.5)
        {
                Rcpx=r0+r1*3.5;
        }
        Fcpxout = Mcpx/Rcpx;

        if (F > Fcpxout)
        {

                Tsol = A1+P*(A2+P*A3); // in oC 
                Tlhl = B1+P*(B2+P*B3);
                Tliq = C1+P*(C2+P*C3);

                Xh2o_bulk=water;

                startF=0.0;
                endF=1.0;
                F=0.0;
                do
                {
                        Xh2o=Xh2o_bulk/(Dh2o+F*(1.0-Dh2o));
                        Xh2o_sat=zeta1*pow(P,ramda)+zeta2*P;
                        if(Xh2o>Xh2o_sat)
                        {
                                Xh2o=Xh2o_sat;
                        }
                        deltaT=K*pow(Xh2o,gama);
                        Tcpxout = pow(Fcpxout,1.0/beta1)*(Tlhl-Tsol)+Tsol-deltaT;
                        Tprime = (T-Tcpxout)/(Tliq-Tcpxout);
                        if(Tprime<0.0)
                        {
                                F_cal=Fcpxout;
                        }
                        else
                        {
                                F_cal = Fcpxout + (1.0-Fcpxout)*pow(Tprime,beta2);
                        }
                        error=F_cal-F;
                        if(error>0)
                        {
                                startF=F;
                        }
                        else
                        {
                                endF=F;
                        }
                        F=0.5*(startF+endF);
                }while(fabs(error)>1e-7 && fabs(startF-endF)>1e-7);
        }        
return F;
}

static float U_dot_grad_F(struct All_variables *E,int e,float *F)
{
const float year=365*24*3600;
float UDF[4][9],rad[9],theta[9],phi[9];
float vtheta,vphi,vrad;
float UdotgradF;
int i,j;
int m=1;
const int dims=E->mesh.nsd;
const int vpts=vpoints[dims];
const int ppts=ppoints[dims];
const int ends=enodes[dims];
const int sphere_key=1;
float VV[4][9],XX[4][9],area;
float velocity_scale=E->data.therm_diff/E->data.radius*year; //kappa/radius in m/year
float R=E->data.radius; //test
void velo_from_element();
float temp1,temp2,temp3;


velo_from_element(E,VV,m,e,sphere_key);

for(j=1;j<=ends;j++)
{
	XX[1][j]=E->sx[m][1][E->ien[m][e].node[j]];
	XX[2][j]=E->sx[m][2][E->ien[m][e].node[j]];
	XX[3][j]=E->sx[m][3][E->ien[m][e].node[j]];
}

vtheta=vphi=vrad=0.0;
for(i=1;i<=8;i++)
{
	vtheta+=VV[1][i]*0.125;
	vphi+=VV[2][i]*0.125;
	vrad+=VV[3][i]*0.125;
}
vtheta*=velocity_scale;
vphi*=velocity_scale;
vrad*=velocity_scale;

for(i=1;i<=vpts;i++)
{
//	vtheta[i]=vphi[i]=vrad[i]=0.0;
	theta[i]=phi[i]=rad[i]=0.0;
	for(j=1;j<=ends;j++)
	{
//		vtheta[i]+=VV[1][j]*E->N.vpt[GNVINDEX(j,i)];
//		vphi[i]+=VV[2][j]*E->N.vpt[GNVINDEX(j,i)];
//		vrad[i]+=VV[3][j]*E->N.vpt[GNVINDEX(j,i)];
		theta[i]+=XX[1][j]*E->N.vpt[GNVINDEX(j,i)];
		phi[i]+=XX[2][j]*E->N.vpt[GNVINDEX(j,i)];
		rad[i]+=XX[3][j]*E->N.vpt[GNVINDEX(j,i)];
	}
}

/*
for(i=1;i<=vpts;i++)
{
	vtheta[i]=vtheta[i]/velocity_scale;
	vphi[i]=vphi[i]/velocity_scale;
	vrad[i]=vrad[i]/velocity_scale;
//	rad[i]=rad[i]*R;
}

*/
/*
for(i=1;i<=vpts;i++)
{
	UDF[i]=0.0;
	for(j=1;j<=ends;j++)
	{
		UDF[i]=UDF[i]+vrad[i]*F[j]*E->gNX[m][e].vpt[GNVXINDEX(2,j,i)]
			+1.0/rad[i]*vtheta[i]*F[j]*E->gNX[m][e].vpt[GNVXINDEX(0,j,i)]
			+1.0/(rad[i]*sin(theta[i]))*vphi[i]*F[j]*E->gNX[m][e].vpt[GNVXINDEX(1,j,i)];
	}
	UDF[i]/=R;
}
*/

for(i=1;i<=vpts;i++)
{
	UDF[1][i]=0.0;
	UDF[2][i]=0.0;
	UDF[3][i]=0.0;
	for(j=1;j<=ends;j++)
	{
		UDF[1][i]+=1.0/rad[i]*F[j]*E->gNX[m][e].vpt[GNVXINDEX(0,j,i)];
		UDF[2][i]+=1.0/(rad[i]*sin(theta[i]))*F[j]*E->gNX[m][e].vpt[GNVXINDEX(1,j,i)];
		UDF[3][i]+=1.0*F[j]*E->gNX[m][e].vpt[GNVXINDEX(2,j,i)];
	}
}

temp1=temp2=temp3=0.0;
for(i=1;i<=vpts;i++)
{
	temp1+=UDF[1][i]*E->gDA[m][e].vpt[i];
	temp2+=UDF[2][i]*E->gDA[m][e].vpt[i];
	temp3+=UDF[3][i]*E->gDA[m][e].vpt[i];
}
temp1/=E->eco[m][e].area;
temp2/=E->eco[m][e].area;
temp3/=E->eco[m][e].area;

//printf("%.6e %.6e %.6e %.6e %.6e %.6e\n",vtheta,temp1,vphi,temp2,vrad,temp3);
//terminate();
UdotgradF=(vtheta*temp1+vphi*temp2+vrad*temp3)/R;

/*
UdotgradF=0.0;
for(i=1;i<=vpts;i++)
{
	UdotgradF+=UDF[i]*E->gDA[m][e].vpt[i];
}

UdotgradF /= E->eco[m][e].area;
*/

return UdotgradF;
}

static double get_element_area(struct All_variables *E, int m, int iel)
{

double area;
int es,i,j,ii,ia[5],lev;
double aa,y1[4],y2[4],angle[6],xx[4][5],area_sphere_cap();
void get_angle_sphere_cap();
float surf_area;
double temp;

        ia[1]=E->ien[m][iel].node[5];
        ia[2]=E->ien[m][iel].node[6];
        ia[3]=E->ien[m][iel].node[7];
        ia[4]=E->ien[m][iel].node[8];

        for(i=1;i<=4;i++)
        {
                xx[1][i] = E->x[m][1][ia[i]]/E->sx[m][3][ia[1]];
                xx[2][i] = E->x[m][2][ia[i]]/E->sx[m][3][ia[1]];
                xx[3][i] = E->x[m][3][ia[i]]/E->sx[m][3][ia[1]];
        }

        get_angle_sphere_cap(xx,angle);

        area = area_sphere_cap(angle);

return area;
}

static float linear_melting(struct All_variables *E,float T, float rad){
	const float tiny=1e-6;
	int Is0=(E->parallel.me==0);
	double *Ts=E->sol[1];
	double *Tl=E->liq[1];
	int j,nz;
	double F,Tsol,Tliq,rad1,rad2;
	T = T;
	for(j=1;j<E->lmesh.noz;j++)  {
		rad1=E->sx[1][3][j];
		rad2=E->sx[1][3][j+1];
		if(rad>rad1-tiny&&rad<rad2+tiny){
			nz = j + E->parallel.me_loc[3] * (E->lmesh.noz - 1);
			Tsol = Ts[nz]*(rad-rad2)/(rad1-rad2)+Ts[nz+1]*(rad-rad1)/(rad2-rad1);
			Tliq = Tl[nz]*(rad-rad2)/(rad1-rad2)+Tl[nz+1]*(rad-rad1)/(rad2-rad1);
			F = (T-Tsol)/(Tliq-Tsol);
			break;
		}
	}
	return F;
}

static float linear_melting_postp(struct All_variables *E,float T, float rad){
	const float tiny=1e-6;
	int Is0=(E->parallel.me==0);
	double *Ts=E->sol[1];
	double *Tl=E->liq[1];
	int j,nz;
	double F,Tsol,Tliq,rad1,rad2;
	double LL = E->lunar.latent_heat/(E->data.Cp*E->data.ref_temperature);
	if(Is0)
		fprintf(stderr,"LL: %.4e\n",LL);
	for(j=1;j<E->lmesh.noz;j++)  {
		rad1=E->sx[1][3][j];
		rad2=E->sx[1][3][j+1];
		if(rad>rad1-tiny&&rad<rad2+tiny){
			nz = j + E->parallel.me_loc[3] * (E->lmesh.noz - 1);
			Tsol = Ts[nz]*(rad-rad2)/(rad1-rad2)+Ts[nz+1]*(rad-rad1)/(rad2-rad1);
			Tliq = Tl[nz]*(rad-rad2)/(rad1-rad2)+Tl[nz+1]*(rad-rad1)/(rad2-rad1);
			F = (T-Tsol)/(Tliq-Tsol+LL);
			break;
		}
	}
	return F;
}

static void correct_ic_sol(struct All_variables *E)
{
	const int number_of_tracers = E->trace.ntracers[1];
	const double t_scale = pow(E->data.radius,2.0)/E->data.therm_diff/(365*24*3600*1e6);
	const int ic_flavor = E->trace.ic_flavor;
	const int Is1 = (E->parallel.me==1);
	const int nflavors = E->trace.nflavors;
	const double ic_inter =  E->chemical.interface[ic_flavor][1];
	const double ic_thick = E->chemical.interface[ic_flavor][2] 
								- E->chemical.interface[ic_flavor][1];
	int kk;
	double t_model;
	t_model =  E->monitor.elapsed_time;
	double ic_sol_inter = ic_inter + ic_thick * t_model/(E->control.t_ic_sol/t_scale);
	if(Is1){
		fprintf(stderr, "correct_ic_sol");
	}
	for (kk=1;kk<=number_of_tracers;kk++) {
		if((E->trace.extraq[1][0][kk]==nflavors-1)&&(E->trace.extraq[1][6][kk]<t_model)){
//		if((E->trace.extraq[1][0][kk]==nflavors-1)&&
//			((E->trace.extraq[1][6][kk]<t_model)||
//			(E->trace.basicq[1][2][kk]<ic_sol_inter))){
			E->trace.extraq[1][0][kk] = ic_flavor;
			/*
			if(Is1){
				fprintf(stderr,"correct_ic_sol, rad = %.4e, t_sol = %.4e\n",
						E->trace.basicq[1][2][kk],E->trace.extraq[1][6][kk]);
			}
			*/
		}
	}
}
