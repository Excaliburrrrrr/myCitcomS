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

#include <mpi.h>

#include <math.h>
#include <sys/types.h>

#include "element_definitions.h"
#include "global_defs.h"
#include "citcom_init.h"
#include "interuption.h"
#include "output.h"
#include "parallel_related.h"
#include "checkpoints.h"

extern int Emergency_stop;

void solver_init(struct All_variables *E);
void calculate_melting_flux_postp(struct All_variables *E);

int main(argc,argv)
     int argc;
     char **argv;

{	/* Functions called by main*/
  void general_stokes_solver();
  void general_stokes_solver_pseudo_surf();
  void global_default_values();
  void read_instructions();
  void initial_setup();
  void initial_conditions();
  void post_processing();
  void read_velocity_boundary_from_file();
  void read_rayleigh_from_file();
  void read_mat_from_file();
  void read_temperature_boundary_from_file();

  void output_finalize();
  void tracer_advection();
  void heat_flux();

  void melting_initialize(); //lhy solidus_liquidus
  void read_sol_liq_from_file();

  void update_cmb_T(); //lhy update_cmb_T
  void boundary_write_q_files(); //lhy alter_write_q_file
  void read_vtk(); //post_process
  //void lunar_test_functions();
  void get_melting_status();	
  void output_coord(struct All_variables *E);

  float cpu_time_on_vp_it;

  int cpu_total_seconds,k,need_init_sol;
  double CPU_time0(),time,initial_time,start_time;

  struct All_variables *E;
  MPI_Comm world;

  MPI_Init(&argc,&argv); /* added here to allow command-line input */

  if (argc < 2)   {
    fprintf(stderr,"Usage: %s PARAMETERFILE\n", argv[0]);
    parallel_process_termination();
  }



  /* this section reads input, allocates memory, and set some initial values;
   * replaced by CitcomS.Controller.initialize() and
   * CitcomS.Solver.initialize() in Pyre. */
  world = MPI_COMM_WORLD;
  E = citcom_init(&world); /* allocate global E and do initializaion here */

  /* define common aliases for full/regional functions */
  solver_init(E);

  start_time = time = CPU_time0();

  /* Global interuption handling routine defined once here */
  set_signal();

  /* default values for various parameters */
  global_default_values(E);

  /* read input parameters from file */
  int Is0 = (E->parallel.me == 0);
  if(Is0){
  	  fprintf(stderr,"start read_instructions\n");
  } //debug
  read_instructions(E, argv[1]);

  /* create mesh, setup solvers etc. */
  if(Is0){
  	  fprintf(stderr,"start initial_setup\n");
  } //debug
  initial_setup(E);

  cpu_time_on_vp_it = CPU_time0();
  initial_time = cpu_time_on_vp_it - time;
  if (E->parallel.me == 0)  {
    fprintf(stderr,"Input parameters taken from file '%s'\n",argv[1]);
    fprintf(stderr,"Initialization complete after %g seconds\n\n",initial_time);
    fprintf(E->fp,"Initialization complete after %g seconds\n\n",initial_time);
    fflush(E->fp);
  }
  /*
  if (Is0){
  	lunar_test_functions(E);
  }*/ //debug models in Lunar_model
  


  /* write all config parameters to a file named pidXXXXXXXXX */
  print_all_config_parameters(E);

  /* this section sets the initial condition;
   * replaced by CitcomS.Controller.launch() ->
   * CitcomS.Solver.launch() in Pyre. */
  if (E->control.post_p) {
      /* the initial condition is from previous checkpoint */
      //read_checkpoint(E);
	  //read_vtk(E,0);
	  //heat_flux(E); //lhy alter_write_q_file
	  if(E->convection.sol_liq){
		if (E->parallel.me == 0){
			fprintf(stderr,"restarting solidus_liquidus\n");
		}
		melting_initialize(E);
		read_sol_liq_from_file(E);
  	    //get_melting_status(E);	
		//boundary_write_q_files(E);
	  }
      need_init_sol = 0;	/*  */
      E->monitor.solution_cycles = E->advection.min_timesteps-1;

      /* the program will finish after post_processing */
//      post_processing(E);
//      (E->problem_output)(E, E->monitor.solution_cycles);

//      citcom_finalize(E, 0);
	  output_coord(E); //output coordinate for use of ploting surface topography
	  goto begin_iterate;
  }
  else if (E->control.restart) {
      /* the initial condition is from previous checkpoint */
      read_checkpoint(E);
      if(E->control.tracer && (E->trace.ic_method_for_flavors == 99)){
	/* 
	   if ggrd tracer input is selected, this will override
	   existing tracers, or allow addition of tracers after
	   restart from a thermal model
	 
	*/
	initialize_tracers(E);
        if (E->composition.on)
	  init_composition(E);
	need_init_sol = 1;
      }else
	need_init_sol = 0;
	
	if(E->convection.sol_liq){
		if (E->parallel.me == 0)
			fprintf(stderr,"restarting solidus_liquidus\n");
		melting_initialize(E);
		read_sol_liq_from_file(E);
	  	if(E->convection.sol_liq){
  	    	get_melting_status(E);	
		}
	}
	heat_flux(E); //lhy alter_write_q_file
	boundary_write_q_files(E);
  }
  else {
      /* regular init, or read T from file only */
	  if(Is0){
		  fprintf(stderr,"start initial_condition\n");
	  }//debug
      initial_conditions(E);
      need_init_sol = E->control.need_init_sol;	/*  */
  }
  if(need_init_sol){
    /* find first solution */
      if(E->control.pseudo_free_surf) {
          if(E->mesh.topvbc == 2)
              general_stokes_solver_pseudo_surf(E);
          else
              assert(0);
      }
      else
          general_stokes_solver(E);
  }

  /* stop the computation if only computes stokes' problem */
  if (E->control.stokes&&(!E->control.post_p))  {
    if(E->control.tracer==1)
      tracer_advection(E);
	heat_flux(E); //lhy alter_write_q_file
	if(E->convection.sol_liq){
  		get_melting_status(E);	
	}
	boundary_write_q_files(E);
    (E->problem_output)(E, E->monitor.solution_cycles);
    citcom_finalize(E, 0);
  }


  (E->problem_output)(E, E->monitor.solution_cycles);

  /* information about simulation time and wall clock time */
  if(!E->control.post_p){
	output_time(E, E->monitor.solution_cycles);

  	if(!E->control.restart)	/* if we have not restarted, print new
				   checkpoint, else leave as is to
				   allow reusing directories */
    	output_checkpoint(E);
   }


 
begin_iterate:;
  /* this section advances the time step;
   * replaced by CitcomS.Controller.march() in Pyre. */
  while ( E->control.keep_going   &&  (Emergency_stop == 0) ) {

    /* The next few lines of code were replaced by
     * pyCitcom_PG_timestep_solve() in Pyre version.
     * If you modify here, make sure its Pyre counterpart
     * is modified as well */

  	if(!E->control.post_p){
    	E->monitor.solution_cycles++;
	}
	else{
		int keep_going = 1;
		char filename[1000];
		while((E->monitor.solution_cycles<E->advection.max_timesteps-1)&&keep_going){
    		E->monitor.solution_cycles++;
			sprintf(filename,"%s/%s.proc0.%d.vts",E->control.data_dir_old,E->control.data_prefix_old,E->monitor.solution_cycles);
			if(!access(filename,0)){
				keep_going = 0;
				if(Is0){
					fprintf(stderr,"file found: %s\n",filename);
				}
			}
		}
	}
	
    if(E->monitor.solution_cycles>E->control.print_convergence)
      E->control.print_convergence=1;

    
    if(!E->control.post_p){
    	(E->next_buoyancy_field)(E);
    }
    /* */


    if(((E->advection.total_timesteps < E->advection.max_total_timesteps) &&
	(E->advection.timesteps < E->advection.max_timesteps)) ||
       (E->advection.total_timesteps < E->advection.min_timesteps) )
      E->control.keep_going = 1;
    else{
	if (Is0){
      		fprintf(E->fp,"max step reached\n");
	}
      	E->control.keep_going = 0;
    }

	if ((E->control.post_p)&&(E->monitor.solution_cycles >= E->advection.max_timesteps-1)){
		if (Is0){
      			fprintf(E->fp,"max time of post-processing reached\n");
		}
		E->control.keep_going = 0;
	}

    cpu_total_seconds = CPU_time0()-start_time;
    if (cpu_total_seconds > E->control.record_all_until)  {
	if (Is0){
      		fprintf(E->fp,"limit of cpu total seconds reached\n");
	}
      E->control.keep_going = 0;
    }

    if (E->monitor.T_interior > E->monitor.T_interior_max_for_exit)  {
      fprintf(E->fp,"quit due to maxT = %.4e sub_iteration%d\n",E->monitor.T_interior,E->advection.last_sub_iterations);
      parallel_process_termination();
    }

	if(!E->control.post_p){
    	if(E->control.tracer==1)
      		tracer_advection(E);
    	general_stokes_solver(E);
	}
	else{	// post process
		read_vtk(E,E->monitor.solution_cycles);
		calculate_melting_flux_postp(E); // melt calculation
	}
	if(Is0){
		fprintf(stderr,"Next heat_flux\n");
	}
	heat_flux(E);
	if(E->convection.sol_liq){
		get_melting_status(E);	
	}
	//fprintf(stderr,"proc:%d E.heat_flux: %.4e\n",E->parallel.me,E->sphere.cap[1].heat_flux); //lhy debug
	if(!E->control.post_p){
		if(E->parallel.me_loc[3]==0&&E->control.coreT)
			update_cmb_T(E); //lhy update_cmb_T
	}
    if(E->output.write_q_files)
      if ((E->monitor.solution_cycles % E->output.write_q_files)==0)
	  	boundary_write_q_files(E);
    
	if(!E->control.post_p){
		if ((E->monitor.solution_cycles % E->control.record_every)==0) {
			(E->problem_output)(E, E->monitor.solution_cycles);
    	}

    /* information about simulation time and wall clock time */
    	output_time(E, E->monitor.solution_cycles);
	

    /* print checkpoint every checkpoint_frequency, unless we have restarted,
       then, we would like to avoid overwriting 
    */
    if ( ((E->monitor.solution_cycles % E->control.checkpoint_frequency)==0) &&
	 ((!E->control.restart) || (E->monitor.solution_cycles != E->monitor.solution_cycles_init))){
	output_checkpoint(E);
    }
    /* updating time-dependent material group
     * if mat_control is 0, the material group has already been
     * initialized in initial_conditions() */
    if(E->control.mat_control==1)
      read_mat_from_file(E);

#ifdef USE_GGRD
    /* updating local rayleigh number (based on Netcdf grds, the
       rayleigh number may be modified laterally in the surface
       layers) */
    /* no counterpart in pyre */
    if(E->control.ggrd.ray_control)
      read_rayleigh_from_file(E);
#endif

    /* updating plate velocity boundary condition */
    if(E->control.vbcs_file==1)
      read_velocity_boundary_from_file(E);

    /* updating plate temperature boundary condition */
    if(E->control.tbcs_file)
      read_temperature_boundary_from_file(E);   
	}
	else{
		(E->problem_output)(E, E->monitor.solution_cycles);
	}


    if (E->parallel.me == 0)  {
      fprintf(E->fp,"CPU total = %g & CPU = %g for step %d time = %.4e dt = %.4e  maxT = %.4e sub_iteration%d\n",CPU_time0()-start_time,CPU_time0()-time,E->monitor.solution_cycles,E->monitor.elapsed_time,E->advection.timestep,E->monitor.T_interior,E->advection.last_sub_iterations);

      time = CPU_time0();
    }

  }


  /* this section prints time accounting;
   * no counterpart in pyre */
  if (E->parallel.me == 0)  {
    fprintf(stderr,"cycles=%d\n",E->monitor.solution_cycles);
    cpu_time_on_vp_it=CPU_time0()-cpu_time_on_vp_it;
    fprintf(stderr,"Average cpu time taken for velocity step = %f\n",
	    cpu_time_on_vp_it/((float)(E->monitor.solution_cycles-E->control.restart)));
    fprintf(E->fp,"Initialization overhead = %f\n",initial_time);
    fprintf(E->fp,"Average cpu time taken for velocity step = %f\n",
	    cpu_time_on_vp_it/((float)(E->monitor.solution_cycles-E->control.restart)));
  }
  citcom_finalize(E, 0);
  return(0);

}
