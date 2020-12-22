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


struct CONVECTION { /* information controlling convection problems */

    int tic_method;
    float half_space_age;

	int tic_max_T ;//lhy 20180104
	float tic_max_T_value;

	int sol_liq; //lhy solidus_liquidus

#define PERTURB_MAX_LAYERS 255
    int number_of_perturbations;
    int perturb_ll[PERTURB_MAX_LAYERS];
    int perturb_mm[PERTURB_MAX_LAYERS];
    int load_depth[PERTURB_MAX_LAYERS];
    float perturb_mag[PERTURB_MAX_LAYERS];

  float blob_center[3];
  float blob_radius;
  float blob_dT;
  int blob_bc_persist;

    struct SOURCES {
	    int number;
	    float t_offset;
	    float *Q;
		float *QER;
		float *QER1;
	    float *lambda;
	}  heat_sources;




} convection;


