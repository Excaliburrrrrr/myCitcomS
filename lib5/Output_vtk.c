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

/* Routine to process the output of the finite element cycles
   and to turn them into a coherent suite of files  */


#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "output.h"
#include "string.h"

#ifdef USE_GZDIR
#include "zlib.h"
#endif

#define CHUNK 16384

void  visc_from_nodes_to_gint();
static void write_binary_array(int nn, float* array, FILE * f);
static void write_ascii_array(int nn, int perLine, float *array, FILE *fp);
static void read_vtk_single(FILE *fp, double *data, int inode);
static void read_vtk_single_float(FILE *fp, float *data, int inode);
static void read_vtk_triple_float(FILE *fp, float **data, int inode);
static void type_of_data(char *type, const char *line);
static void cart_sph_project(float vx,float vy,float vz,float theta,float fi,float *vth,float *vfi,float *vr);

static void vts_file_header(struct All_variables *E, FILE *fp)
{

    const char format[] =
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"StructuredGrid\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
        "  <StructuredGrid WholeExtent=\"%s\">\n"
        "    <Piece Extent=\"%s\">\n";

    char extent[64], header[1024];

    snprintf(extent, 64, "%d %d %d %d %d %d",
             E->lmesh.ezs, E->lmesh.ezs + E->lmesh.elz,
             E->lmesh.exs, E->lmesh.exs + E->lmesh.elx,
             E->lmesh.eys, E->lmesh.eys + E->lmesh.ely);

    snprintf(header, 1024, format, extent, extent);

    fputs(header, fp);

    return;
}


static void vts_file_trailer(struct All_variables *E, FILE *fp)
{
    const char trailer[] =
        "    </Piece>\n"
        "  </StructuredGrid>\n"
        "</VTKFile>\n";

    fputs(trailer, fp);

    return;
}


static void vtk_point_data_header(struct All_variables *E, FILE *fp)
{
    fputs("      <PointData Scalars=\"temperature\" Vectors=\"velocity\">\n", fp);
    return;
}


static void vtk_point_data_trailer(struct All_variables *E, FILE *fp)
{
    fputs("      </PointData>\n", fp);
    return;
}


static void vtk_cell_data_header(struct All_variables *E, FILE *fp)
{
    fputs("      <CellData>\n", fp);
    return;
}


static void vtk_cell_data_trailer(struct All_variables *E, FILE *fp)
{
    fputs("      </CellData>\n", fp);
    return;
}


static void vtk_output_temp(struct All_variables *E, FILE *fp)
{
    int i;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    float* floattemp = malloc(nodes*sizeof(float));

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"temperature\" format=\"%s\">\n", E->output.vtk_format);

    for(i=0;i < nodes;i++)
        floattemp[i] =  (float) *(E->T[1]+i+1);

    if (strcmp(E->output.vtk_format,"binary") == 0) {
        write_binary_array(nodes,floattemp,fp);
    } else {
        write_ascii_array(nodes,1,floattemp,fp);
    }
    fputs("        </DataArray>\n", fp);
    free(floattemp);
    return;
}

static void vtk_output_diff_temp(struct All_variables *E, FILE *fp)
{
    int i, noz;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    float* floattemp = malloc(nodes*sizeof(float));

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"diff_temperature\" format=\"%s\">\n", E->output.vtk_format);
    for(i=0;i < nodes;i++){
		noz = i%E->lmesh.noz+1;
        floattemp[i] =  (float) *(E->T[1]+i+1) - E->Have.T[noz];
	}
    if (strcmp(E->output.vtk_format,"binary") == 0) {
        write_binary_array(nodes,floattemp,fp);
    } else {
        write_ascii_array(nodes,1,floattemp,fp);
    }
    fputs("        </DataArray>\n", fp);
    free(floattemp);
    return;
}


/*lhy solidus_liquidus
 * melting percent output
 */
static void vtk_output_melting(struct All_variables *E, FILE *fp)
{
    int i;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    float* melting = malloc(nodes*sizeof(float));

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"melting\" format=\"%s\">\n", E->output.vtk_format);
	//for(i=1;i <= nodes;i++) E->melting[1][i] = E->T[1][i];
    for(i=0;i < nodes;i++)		//melting[i] =  (float) *(E->T[1]+i+1);
        melting[i] =  (float) *(E->melting[1]+i+1);

    if (strcmp(E->output.vtk_format,"binary") == 0) {
        write_binary_array(nodes,melting,fp);
    } else {
        write_ascii_array(nodes,1,melting,fp);
    }
    fputs("        </DataArray>\n", fp);
    free(melting);
    return;
}

static void vtk_output_velo(struct All_variables *E, FILE *fp)
{
    int i, j;
    int nodes=E->sphere.caps_per_proc*E->lmesh.nno;
    double sint, sinf, cost, cosf;
    float *V[4];
    const int lev = E->mesh.levmax;
    float* floatvel = malloc(nodes*3*sizeof(float));

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"%s\">\n", E->output.vtk_format);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        V[1] = E->sphere.cap[j].V[1];
        V[2] = E->sphere.cap[j].V[2];
        V[3] = E->sphere.cap[j].V[3];

        for(i=1; i<=E->lmesh.nno; i++) {
            sint = E->SinCos[lev][j][0][i];
            sinf = E->SinCos[lev][j][1][i];
            cost = E->SinCos[lev][j][2][i];
            cosf = E->SinCos[lev][j][3][i];

            floatvel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3+0] = (float)(V[1][i]*cost*cosf - V[2][i]*sinf + V[3][i]*sint*cosf);
            floatvel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3+1] = (float)(V[1][i]*cost*sinf + V[2][i]*cosf + V[3][i]*sint*sinf);
            floatvel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3+2] = (float)(-V[1][i]*sint + V[3][i]*cost);
        }
    }

    if (strcmp(E->output.vtk_format, "binary") == 0)
        write_binary_array(nodes*3,floatvel,fp);
    else
        write_ascii_array(nodes*3,3,floatvel,fp);
    fputs("        </DataArray>\n", fp);

    free(floatvel);
    return;
}

static void vtk_output_tide(struct All_variables *E, FILE *fp)
{
    int i, j;
    int nodes=E->sphere.caps_per_proc*E->lmesh.nno;
    double sint, sinf, cost, cosf;
    const int lev = E->mesh.levmax;
    float* floatvel = malloc(nodes*3*sizeof(float));


    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"tidal_force\" NumberOfComponents=\"3\" format=\"%s\">\n", E->output.vtk_format);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->lmesh.nno; i++) {
            sint = E->SinCos[lev][j][0][i];
            sinf = E->SinCos[lev][j][1][i];
            cost = E->SinCos[lev][j][2][i];
            cosf = E->SinCos[lev][j][3][i];

            floatvel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3+0] = (float)(E->tide_force[j][i]*cost*cosf);
            floatvel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3+1] = (float)(E->tide_force[j][i]*cost*sinf);
            floatvel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3+2] = (float)(-E->tide_force[j][i]*sint);
        }
    }

    if (strcmp(E->output.vtk_format, "binary") == 0)
        write_binary_array(nodes*3,floatvel,fp);
    else
        write_ascii_array(nodes*3,3,floatvel,fp);
    fputs("        </DataArray>\n", fp);

    free(floatvel);
    return;
}


static void vtk_output_visc(struct All_variables *E, FILE *fp)
{
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    int lev = E->mesh.levmax;

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"viscosity\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0) {
            write_binary_array(nodes,&E->VI[lev][1][1],fp);
        } else {
            write_ascii_array(nodes,1,&E->VI[lev][1][1],fp);
        }

    fputs("        </DataArray>\n", fp);
    return;
}

static void vtk_output_visc_extra(struct All_variables *E, FILE *fp) //lhy new_smoothing_method
{
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    int lev = E->mesh.levmax;

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"viscosity_extra\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0) {
            write_binary_array(nodes,&E->VIextra[1][1],fp);
        } else {
            write_ascii_array(nodes,1,&E->VIextra[1][1],fp);
        }

    fputs("        </DataArray>\n", fp);
    return;
}


static void vtk_output_coord(struct All_variables *E, FILE *fp)
{
    /* Output Cartesian coordinates as most VTK visualization softwares
       assume it. */
    int i, j;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    float* floatpos = malloc(nodes*3*sizeof(float));

    fputs("      <Points>\n", fp);
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"%s\">\n", E->output.vtk_format);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->lmesh.nno; i++){
                floatpos[((j-1)*E->lmesh.nno+i-1)*3] = (float)(E->x[j][1][i]);
	        floatpos[((j-1)*E->lmesh.nno+i-1)*3+1]=(float)(E->x[j][2][i]);
	        floatpos[((j-1)*E->lmesh.nno+i-1)*3+2]=(float)(E->x[j][3][i]);
        }
    }

    if (strcmp(E->output.vtk_format, "binary") == 0)
        write_binary_array(nodes*3,floatpos,fp);
    else
        write_ascii_array(nodes*3,3,floatpos,fp);
    fputs("        </DataArray>\n", fp);
    fputs("      </Points>\n", fp);
    free(floatpos);
    return;
}

static void vtk_output_stress(struct All_variables *E, FILE *fp)
{
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
   /* for stress computation */
    void allocate_STD_mem();
    void compute_nodal_stress();
    void free_STD_mem();
    float *SXX[NCS],*SYY[NCS],*SXY[NCS],*SXZ[NCS],*SZY[NCS],*SZZ[NCS];
    float *divv[NCS],*vorv[NCS];

    /* those are sorted like stt spp srr stp str srp  */
    allocate_STD_mem(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);
    compute_nodal_stress(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);
    free_STD_mem(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"stress\" NumberOfComponents=\"6\" format=\"%s\">\n", E->output.vtk_format);

    if (strcmp(E->output.vtk_format, "binary") == 0) {
        write_binary_array(nodes*6,&E->gstress[1][1],fp);
    } else {
        write_ascii_array(nodes*6,6,&E->gstress[1][1],fp);
    }

    fputs("        </DataArray>\n", fp);
    return;
}

static void vtk_output_comp_nd(struct All_variables *E, FILE *fp)
{
    int i, j, k;
    char name[255];
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    float* floatcompo = malloc (nodes*sizeof(float));
	/*compute and output composition0*/
    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->lmesh.nno; i++) {
        	floatcompo[(j-1)*E->lmesh.nno+i-1] = 1.0;
	    }
    }
	for(k=0;k<E->composition.ncomp;k++) {
        for(j=1; j<=E->sphere.caps_per_proc; j++) {
            for(i=1; i<=E->lmesh.nno; i++) {
                floatcompo[(j-1)*E->lmesh.nno+i-1] -= (float) (E->composition.comp_node[j][k][i]);
	    	}
        }
	}
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"composition0\" format=\"%s\">\n", E->output.vtk_format);
	if (strcmp(E->output.vtk_format, "binary") == 0)
		write_binary_array(nodes,floatcompo,fp);
    else
        write_ascii_array(nodes,1,floatcompo,fp);
    fputs("        </DataArray>\n", fp);
	/*output other compositions*/
    for(k=0;k<E->composition.ncomp;k++) {
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"composition%d\" format=\"%s\">\n", k+1, E->output.vtk_format);

        for(j=1; j<=E->sphere.caps_per_proc; j++) {
            for(i=1; i<=E->lmesh.nno; i++) {
                floatcompo[(j-1)*E->lmesh.nno+i-1] = (float) (E->composition.comp_node[j][k][i]);
	    }
        }

        if (strcmp(E->output.vtk_format, "binary") == 0)
            write_binary_array(nodes,floatcompo,fp);
        else
            write_ascii_array(nodes,1,floatcompo,fp);
        fputs("        </DataArray>\n", fp);
    }
    free(floatcompo);
    return;
}


static void vtk_output_surf(struct All_variables *E,  FILE *fp, int cycles)
{
    int i, j, k;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    char output_file[255];
    float* floattopo = malloc (nodes*sizeof(float));

    if((E->output.write_q_files == 0) || (cycles == 0) ||
      (cycles % E->output.write_q_files)!=0)
        heat_flux(E);
  /* else, the heat flux will have been computed already */

    if(E->control.use_cbf_topo){
        get_CBF_topo(E,E->slice.tpg,E->slice.tpgb);
    }
    else{
        get_STD_topo(E,E->slice.tpg,E->slice.tpgb,E->slice.divg,E->slice.vort,cycles);
    }

    fprintf(fp,"        <DataArray type=\"Float32\" Name=\"surface\" format=\"%s\">\n", E->output.vtk_format);

    for(j=1;j<=E->sphere.caps_per_proc;j++){
        for(i=1;i<=E->lmesh.nsf;i++){
            for(k=1;k<=E->lmesh.noz;k++){
                floattopo[(j-1)*E->lmesh.nno + (i-1)*E->lmesh.noz + k-1] = 0.0;
            }

            if (E->parallel.me_loc[3]==E->parallel.nprocz-1) {

                /* choose either STD topo or pseudo-free-surf topo */
                if(E->control.pseudo_free_surf)
                floattopo[(j-1)*E->lmesh.nno + i*E->lmesh.noz-1] = E->slice.freesurf[j][i];
                else
                floattopo[(j-1)*E->lmesh.nno + i*E->lmesh.noz-1] = E->slice.tpg[j][i];

            }
        }
    }

    if (strcmp(E->output.vtk_format, "binary") == 0)
        write_binary_array(nodes,floattopo,fp);
    else
        write_ascii_array(nodes,1,floattopo,fp);

    fputs("        </DataArray>\n", fp);
  return;
}


static void write_vtm(struct All_variables *E, int cycles)
{
    FILE *fp;
    char vtm_file[255];
    int n;

    const char header[] =
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
        "  <vtkMultiBlockDataSet>\n";

    snprintf(vtm_file, 255, "%s.%d.vtm",
             E->control.data_file, cycles);
    fp = output_open(vtm_file, "w");
    fputs(header, fp);

    for(n=0; n<E->parallel.nproc; n++) {
        fprintf(fp, "    <DataSet index=\"%d\" file=\"%s.proc%d.%d.vts\"/>\n",
                n, E->control.data_prefix, n, cycles);
    }
    fputs("  </vtkMultiBlockDataSet>\n",fp);
    fputs("</VTKFile>",fp);

    fclose(fp);
}

static void write_visit(struct All_variables *E, int cycles)
{
    FILE *fp;
    char visit_file[255];
    int n;

    const char header[] = "!NBLOCKS %d\n";

    snprintf(visit_file, 255, "%s.%d.visit",
             E->control.data_file, cycles);
    fp = output_open(visit_file, "w");
    fprintf(fp, header, E->parallel.nproc);

    for(n=0; n<E->parallel.nproc; n++) {
        fprintf(fp, "%s.proc%d.%d.vts\n",
                E->control.data_prefix, n, cycles);
    }
    fclose(fp);
}

static void write_pvts(struct All_variables *E, int cycles)
{
    FILE *fp;
    char pvts_file[255];
    int i,j,k;
    snprintf(pvts_file, 255, "%s.%d.pvts",
             E->control.data_file,cycles);
    fp = output_open(pvts_file, "w");

    const char format[] =
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
        "  <PStructuredGrid WholeExtent=\"%s\" GhostLevel=\"#\">\n"
        "    <PPointData Scalars=\"temperature\" Vectors=\"velocity\">\n"
        "      <DataArray type=\"Float32\" Name=\"temperature\" format=\"%s\"/>\n"
        "      <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"%s\"/>\n"
        "      <DataArray type=\"Float32\" Name=\"viscosity\" format=\"%s\"/>\n";

    char extent[64], header[1024];

    snprintf(extent, 64, "%d %d %d %d %d %d",
        E->lmesh.ezs, E->lmesh.ezs + E->lmesh.elz*E->parallel.nprocz,
        E->lmesh.exs, E->lmesh.exs + E->lmesh.elx*E->parallel.nprocx,
        E->lmesh.eys, E->lmesh.eys + E->lmesh.ely*E->parallel.nprocy);

    snprintf(header, 1024, format, extent, E->output.vtk_format,
             E->output.vtk_format, E->output.vtk_format);
    fputs(header, fp);

    if (E->output.stress){
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"stress\" NumberOfComponents=\"6\" format=\"%s\"/>\n", E->output.vtk_format);
    }
    if (E->output.comp_nd && E->composition.on){
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"composition1\" format=\"%s\"/>\n", E->output.vtk_format);
    }
    if (E->output.surf){
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"surface\" format=\"%s\"/>\n", E->output.vtk_format);
    }

    fputs("    </PPointData>\n \n"
    "    <PCellData>\n"
    "    </PCellData>\n \n"
    "    <PPoints>\n"
    "      <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"binary\" />\n"
    "    </PPoints>\n", fp);

    for(i=0; i < E->parallel.nprocy;i++){
        for(j=0; j < E->parallel.nprocx;j++){
            for(k=0; k < E->parallel.nprocz;k++){
                fprintf(fp, "    <Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s.proc%d.%d.vts\"/>\n",
                    (k%E->parallel.nprocz)*E->lmesh.elz,
                    (k%E->parallel.nprocz+1)*E->lmesh.elz,
                    (j%E->parallel.nprocx)*E->lmesh.elx, (j%E->parallel.nprocx+1)*E->lmesh.elx,
                    (i%E->parallel.nprocy)*E->lmesh.ely, (i%E->parallel.nprocy+1)*E->lmesh.ely,
                    E->control.data_prefix,
                    i*E->parallel.nprocx*E->parallel.nprocz+j*E->parallel.nprocz+k, cycles);
            }
        }
    }

    fputs("  </PStructuredGrid>\n",fp);
    fputs("</VTKFile>",fp);

    fclose(fp);
}

static void write_ascii_array(int nn, int perLine, float *array, FILE *fp)
{
    int i;

    switch (perLine) {
    case 1:
        for(i=0; i<nn; i++)
            fprintf(fp, "%.4e\n", array[i]);
        break;
    case 3:
        for(i=0; i < nn/3; i++)
            fprintf(fp,"%.4e %.4e %.4e\n",array[3*i],array[3*i+1],array[3*i+2]);
        break;
    case 6:
        for(i=0; i < nn/6; i++)
            fprintf(fp,"%.4e %.4e %.4e %.4e %.4e %.4e\n",
                    array[6*i],array[6*i+1],array[6*i+2],
                    array[6*i+3],array[6*i+4],array[6*i+5]);
        break;
    }
    return;
}

static void FloatToUnsignedChar(float * floatarray, int nn, unsigned char * chararray)
{
    /* simple float to unsigned chararray routine via union
    nn=length(intarray) chararray has to be BIG ENOUGH! */
    int i;
    union FloatToUnsignedChars
        {
            float input;
            unsigned char output[4];
        } floattransform;

    for (i=0; i<nn; i++){
        floattransform.input=floatarray[i];
        chararray[4*i]=floattransform.output[0];
        chararray[4*i+1]=floattransform.output[1];
        chararray[4*i+2]=floattransform.output[2];
        chararray[4*i+3]=floattransform.output[3];
    }
    return;
}

static void IntToUnsignedChar(int * intarray, int nn, unsigned char * chararray)
{
    /* simple int - to unsigned chararray routine via union
    nn=length(intarray) chararray has to be BIG ENOUGH! */
    int i;
    union IntToUnsignedChars
        {
            int input;
            unsigned char output[4];
        } inttransform;

    for (i=0; i<nn; i++){
        inttransform.input=intarray[i];
        chararray[4*i]=inttransform.output[0];
        chararray[4*i+1]=inttransform.output[1];
        chararray[4*i+2]=inttransform.output[2];
        chararray[4*i+3]=inttransform.output[3];
        }
}


static void zlibcompress(unsigned char* in, int nn, unsigned char** out, int *nn2)
/* function to compress "in" to "out" reducing size from nn to nn2 */
{
#ifdef USE_GZDIR
    int ntemp=0;

    /* in and out of z-stream */
    unsigned char inz[CHUNK];
    unsigned char outz[CHUNK];

    /* compression level */
    int level = Z_DEFAULT_COMPRESSION;
    int ret,flush;
    int i,j,k;

    /* zlib compression stream */
    z_stream strm;

    /* hope compressed data will be <= uncompressed */
    *out = malloc(sizeof(unsigned char)*nn);

    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;

    /* zlib init */
    ret = deflateInit(&strm, level);
    if (ret == Z_OK){
        i=0;     // position in "in" array
        do{
            j=0; // position in "inz"
            do{
                inz[j++]=in[i++];
            } while((j<CHUNK) && (i<nn)); // stopps if "inz"-buffer is full or "in" array empty
            strm.avail_in=j;              // set number of input chars

            flush = (i==nn) ? Z_FINISH : Z_NO_FLUSH; // done?
            strm.next_in = inz;           // set input buffer

            do{
                strm.avail_out = CHUNK;   // set number of max output chars
                strm.next_out = outz;     // set output buffer

                /* zlib compress */
                ret = deflate(&strm, flush);
                assert(ret != Z_STREAM_ERROR);

                /* zlib changed strm.avail_out=CHUNK
                 to the number of chars we can NOT use
                 in outz */

                for (k=0;k<CHUNK-strm.avail_out;k++){
                    (*out)[ntemp+k]=outz[k];
                }

                /* increase position in "out" */
                ntemp+=(CHUNK-strm.avail_out);
            }while(strm.avail_out==0);
            assert(strm.avail_in == 0);

        }while (flush != Z_FINISH);
    }
    else{fprintf(stderr,"Error during compression init\n");}

    // now we know how short "out" should be!
    *nn2=ntemp;
    *out = realloc(*out,sizeof(unsigned char)*ntemp);

    (void)deflateEnd(&strm);
#endif

    return;
}

static void base64(unsigned char * in, int nn, unsigned char* out)
{
    /*takes *in*-array and "in"-length-"nn" and fills "out"-array
    with base64(in) "out" needs to be big enough!!!
    length(out) >= 4* |^ nn/3.0 ^| */
    char cb64[]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    int len;
    int i;

    for (i=0; i < nn; i+=3){

        len = (3 < nn-i ? 3 : nn-i);
        if (len >= 3){
        /* normal base64 encoding */
            out[i/3*4+0] = cb64[ in[i] >> 2 ];
            out[i/3*4+1] = cb64[ ((in[i] & 0x03) << 4) | ((in[i+1] & 0xf0) >> 4) ];
            out[i/3*4+2] = cb64[ ((in[i+1] & 0x0f) << 2) | ((in[i+2] & 0xc0) >> 6)];
            out[i/3*4+3] = cb64[ in[i+2] & 0x3f ];
        } else if (len == 2){
        /* at the end of array fill up with '=' */
            out[i/3*4+0] = cb64[ in[i] >> 2 ];
            out[i/3*4+1] = cb64[ ((in[i] & 0x03) << 4) | ((in[i+1] & 0xf0) >> 4) ];
            out[i/3*4+2] = cb64[((in[i+1] & 0x0f) << 2)];
            out[i/3*4+3] = (unsigned char) '=';
        } else if (len == 1){
        /* at the end of array fill up with '=' */
            out[i/3*4+0] = cb64[ in[i] >> 2 ];
            out[i/3*4+1] = cb64[ ((in[i] & 0x03) << 4) ];
            out[i/3*4+2] = (unsigned char) '=';
            out[i/3*4+3] = (unsigned char) '=';
        }
    }
}


static void base64plushead(unsigned char * in, int nn, int orinn, unsigned char* out)
{
    /* writing vtk compatible zlib compressed base64 encoded data to "out" */
    int i;
    unsigned char * b64head;
    int b64bodylength;
    unsigned char * b64body;
    /* header of data */
    unsigned char * charhead = malloc(sizeof(unsigned char)*16);
    /* - consists of "1" (number of pieces) */
    /* - original datalength in byte */
    /* - original datalength in byte */
    /* - new datalength after z-lib compression */
    int * headInts= malloc(sizeof(int)*4);
    headInts[0]=1;
    headInts[1]=orinn;
    headInts[2]=orinn;
    headInts[3]=nn;
    // transform to unsigned char
    IntToUnsignedChar(headInts,4,charhead);

    // base64: 16byte -> 24byte
    b64head =  malloc(sizeof(unsigned char)*24);
    // fills b64head
    base64(charhead, 16, b64head);

    // base64 data
    b64bodylength = 4*ceil((double) nn/3.0);
    b64body = malloc(sizeof(unsigned char)*b64bodylength);
    // writes base64 data to b64body
    base64(in,nn,b64body);

    // combines header and body
    for (i=0; i<24 ; i++){
        out[i]=b64head[i];
    }

    for (i=0; i<b64bodylength ; i++){
        out[24+i]=b64body[i];
    }

    if(b64body){free(b64body);}
    if(b64head){free(b64head);}
    if(headInts){free(headInts);}
    if(charhead){free(charhead);}
}

static void write_binary_array(int nn, float* array, FILE * f)
{
    /* writes vtk-data array of floats and performs zip and base64 encoding */
    int chararraylength=4*nn;	/* nn floats -> 4*nn unsigned chars */
    unsigned char * chararray = malloc (chararraylength * sizeof(unsigned char));
    int compressedarraylength = 0;
    unsigned char * compressedarray;
    unsigned char ** pointertocompressedarray= &compressedarray;
    int base64plusheadlength;
    unsigned char * base64plusheadarray;

    FloatToUnsignedChar(array,nn,chararray);

    /* compression routine */
    zlibcompress(chararray,chararraylength,pointertocompressedarray,&compressedarraylength);

    /* special header for zip compressed and bas64 encoded data
    header needs 4 int32 = 16 byte -> 24 byte due to base64 (4*16/3) */
    base64plusheadlength = 24 + 4*ceil((double) compressedarraylength/3.0);
    base64plusheadarray = malloc(sizeof(unsigned char)* base64plusheadlength);

    /* fills base64plusheadarray with everything ready for simple writing */
    base64plushead(compressedarray,compressedarraylength, chararraylength, base64plusheadarray);

    fwrite(base64plusheadarray,sizeof(unsigned char),base64plusheadlength,f);
    fprintf(f,"\n");
    free(chararray);
    free(base64plusheadarray);
    free(compressedarray);
}

/**********************************************************************/

static void cart_sph_project(float vx,float vy,float vz,float theta,float fi,float *vth,float *vfi,float *vr)
{
	*vth = cos(theta)*cos(fi)*vx+cos(theta)*sin(fi)*vy-sin(theta)*vz;
	*vfi = -sin(fi)*vx+cos(fi)*vy;
	*vr = sin(theta)*cos(fi)*vx+sin(theta)*sin(fi)*vy+cos(theta)*vz;
	return;
}

void vtk_output(struct All_variables *E, int cycles)
{
	void output_total_melting(struct All_variables *E);
    char output_file[255];
    FILE *fp;

	output_horiz_avg(E, cycles);
  if(cycles ==0)     //zwb 20200914
    output_coord(E);

    snprintf(output_file, 255, "%s.proc%d.%d.vts",
             E->control.data_file, E->parallel.me, cycles);
    fp = output_open(output_file, "w");

    /* first, write volume data to vts file */
    vts_file_header(E, fp);

    /* write node-based field */
    vtk_point_data_header(E, fp);

    vtk_output_temp(E, fp);

    vtk_output_diff_temp(E, fp);

	if (E->convection.sol_liq){
		vtk_output_melting(E, fp);
	}

    vtk_output_velo(E, fp);

    vtk_output_visc(E, fp);
	if(E->viscosity.SMOOTH&&E->viscosity.SMOOTH_test)
 		vtk_output_visc_extra(E, fp); //lhy new_smoothing_method

	/*if (E->convection.sol_liq){
		vtk_output_melting(E, fp);
	}*/

    if (E->output.stress)
        vtk_output_stress(E, fp);

    if (E->output.comp_nd && E->composition.on)
        vtk_output_comp_nd(E, fp);

    if (E->output.surf)
        vtk_output_surf(E, fp, cycles);

	if (E->output.tide)
		vtk_output_tide(E,fp);

    vtk_point_data_trailer(E, fp);

    /* write element-based field */
    vtk_cell_data_header(E, fp);
    /* TODO: comp_el, heating */
    vtk_cell_data_trailer(E, fp);

    /* write coordinate */
    vtk_output_coord(E, fp);

    vts_file_trailer(E, fp);
    fclose(fp);

	output_surf_botm(E, cycles);
	if (E->output.geoid)      /* this needs to be called after the surface
		and bottom topo has been computed! */
		output_geoid(E, cycles);
	if (E->output.field) //20170630 lhy spectrum
		output_field(E,cycles);

    if (E->parallel.me == 0) {
        if (E->sphere.caps == 12) {
            /* VTK multiblock file */
            write_vtm(E, cycles);
            /* LLNL VisIt multiblock file */
            write_visit(E, cycles);
        }
        else
            write_pvts(E, cycles);
    }
    return;
}
/*post_process lhy 20180215
 * reading vtk files
 * determine types of data by indicators in vtk files and then read in
 * data that follows*/
void read_vtk(struct All_variables *E, int cycles){
	int nno, ii, i, j;
    int lev,find;
	float vth, vfi, vr, vx, vy, vz, theta, fi;
	char input_file[255], inputs_s[1000], type[100], name[100];
	FILE *fp;
	int Is0 = (E->parallel.me == 0);

	if (E->parallel.me == 0) {
		fprintf(stderr,"read_vtk, step: %d\n", E->monitor.solution_cycles);
	}
    nno = E->sphere.caps_per_proc*E->lmesh.nno;
    lev = E->mesh.levmax;
	ii = E->monitor.solution_cycles;

	sprintf(input_file,"%s.proc%d.%d.vts",E->control.old_P_file,E->parallel.me,ii);
  	fp=fopen(input_file,"r");
  	if (fp == NULL) {
    	fprintf(E->fp,"(void read_vtk #1) Cannot open %s\n",input_file);
    	fprintf(stderr,"(void read_vtk #1) Cannot open %s\n",input_file);
    	parallel_process_termination();
  	}
	for (i=0;i<6;i++){
		fgets(inputs_s,1000,fp);
	}
	while (!feof(fp)){
		type_of_data(type, inputs_s);
		find=0; //if find this type, set to 1
	/*	if (E->parallel.me == 0) {
			fprintf(stderr,"Type: %s\n",type); //lhy debug
		}*/
		if (!strcmp("temperature",type)){
			read_vtk_single(fp,E->T[1],nno);
			find=1;
		}
		else if (!strcmp("viscosity",type)){
			read_vtk_single_float(fp,E->VI[lev][1],nno); //read viscosity
      		visc_from_nodes_to_gint(E,E->VI[lev],E->EVI[lev],E->mesh.levmax);//project to g-point
			find=1;
		}
		else if (!strcmp("velocity",type)){
			read_vtk_triple_float(fp,E->sphere.cap[1].V,nno); //read velocity
			for (j=1;j<=nno;j++){
				//convert from cartisian to sphere
				theta = (float)E->sx[1][1][j];
				fi = (float)E->sx[1][2][j];
				vx = E->sphere.cap[1].V[1][j];
				vy = E->sphere.cap[1].V[2][j];
				vz = E->sphere.cap[1].V[3][j];
				cart_sph_project(vx,vy,vz,theta,fi,&vth,&vfi,&vr);
				E->sphere.cap[1].V[1][j] = vth;
				E->sphere.cap[1].V[2][j] = vfi;
				E->sphere.cap[1].V[3][j] = vr;
			}
			find=1;
		}
		else {
			/*read in chemical composition*/
			for (i=0;i<E->composition.ncomp;i++){
				sprintf(name,"composition%d",i+1);
				if (!strcmp(name,type)){
					read_vtk_single(fp,E->composition.comp_node[1][i],nno);
					find=1;
					break;
				}
			}
		}
		if (!(find||(type[0] == '\0'))){
			for (i=0;i<nno;i++){
				fgets(inputs_s,1000,fp);
			}
		}
		fgets(inputs_s,1000,fp);
	}
	fclose(fp);
	return;
}
/*post_process lhy 20180105 read vtk data that is a value
 * read vtk data that is a float value
 * data should be a pointer that has dimension (inode+1)*/
static void read_vtk_single(FILE *fp, double *data, int inode){
	int ii;
	float tempV;
	char input_s[1000];
	for (ii=1;ii<=inode;ii++){
		fgets(input_s,1000,fp);
		if(sscanf(input_s,"%g",&(tempV)) != 1) {
			fprintf(stderr,"Error while reading scalar data from file\n");
			exit(8);
		}
		data[ii]=tempV;
	}
	return;
}
/*post_process lhy 20180105
 * read vtk data that is a float value
 * data should be a pointer that has dimension (inode+1)*/
static void read_vtk_single_float(FILE *fp, float *data, int inode){
	int ii;
	float tempV;
	char input_s[1000];
	for (ii=1;ii<=inode;ii++){
		fgets(input_s,1000,fp);
		if(sscanf(input_s,"%g",&(tempV)) != 1) {
			fprintf(stderr,"Error while reading scalar data from file\n");
			exit(8);
		}
		data[ii]=tempV;
	}
	return;
}
/*post_process lhy 20180105
 * read vtk data that is a 3-float-vector
 * data should be a double pointer that has dimention 4*(inode+1)*/
static void read_vtk_triple_float(FILE *fp, float **data, int inode){
	int ii,jj;
	float tempV[3];
	char input_s[1000];
	for (ii=1;ii<=inode;ii++){
		fgets(input_s,1000,fp);
		if(sscanf(input_s,"%g %g %g",&(tempV[0]),&(tempV[1]),&(tempV[2])) != 3) {
			fprintf(stderr,"Error while reading vector data from file\n");
			exit(8);
		}
		data[1][ii] = tempV[0];
		data[2][ii] = tempV[1];
		data[3][ii] = tempV[2];
	}
	return;
}
/*post_process lhy 20180215
 * find location of different types of data
 * name of types include "temperature", "viscosity", "composition" and "velocity"
 * data follows location of type identification will be read by "void read_vtk" later on*/
static void type_of_data(char *type, const char *line){
	int i,j,k1,k2;
	char tempc[5];
	//find the location of indifier "Name"
	for (i=0;i<strlen(line)-4;i++){
		for (j=0;j<4;j++){
			tempc[j] = line[i+j];
		}
		tempc[4]='\0';
		if (!strcmp("Name", tempc)){
			break;
		}
	}
	//save string between two '"' in type and return
	for (k1=i;k1<strlen(line);k1++){
		if (line[k1] == '"'){
			break;
		}
	}
	for (k2=k1+1;k2<strlen(line);k2++){
		if (line[k2] == '"'){
			break;
		}
	}
	for (i=0;i<k2-k1-1;i++){
		type[i] = line[k1+i+1];
	}
	type[k2-k1-1] = '\0';
	return;
}
