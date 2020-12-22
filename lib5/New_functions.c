/*functions difined by haoyaunli goes here*/

#include <stdlib.h>
#include <stdint.h>
#include <mpi.h>
#include <math.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "parsing.h"
#include "output.h"
void test_functions(struct All_variables *E, int e){
	fprintf(stderr,"test started\n");
	void global_coordinate_from_local(double *X, struct All_variables *E, double* epzl, int e);
	int check_points_in_element(struct All_variables *E, double* X, int e);
	void local_coordinate_from_global(double *epzl, struct All_variables *E, double *X, int e, int ndivide);
	double *XX, *epzl, *epzl1;
	int ele,ival;
//initiate variables
	XX = (double*)malloc(4*sizeof(double));
	epzl = (double*)malloc(4*sizeof(double));
	epzl1 = (double*)malloc(4*sizeof(double));
//assign nature coordinate
	epzl[1] = 0.5;
	epzl[2] = 0.5;
	epzl[3] = 0.5;
//derive global coordinate
	global_coordinate_from_local(XX, E, epzl, e);
//	fprintf(stderr,"points: x = %.4e, y = %.4e, z = %.4e\n",XX[1],XX[2],XX[3]); //debug
//find which element is this point in
	for (ele=1;ele<E->lmesh.nel;ele++){
		ival = check_points_in_element(E,XX,ele);
		if(ival){
			fprintf(stderr,"point in element: %d\n",ele);
			break;
		}	
	}
//derive local coordinate invertly
	local_coordinate_from_global(epzl1, E, XX, e, 10);
	fprintf(stderr,"epzl1 = %.4e, epzl2 = %.4e, epzl3 = %.4e\n", epzl1[1], epzl1[2], epzl1[3]);//stderr
//free space
	free(XX);
	free(epzl);
	free(epzl1);
	fprintf(stderr,"test ended\n");
}

/*get global coordinate of a point from its nature coordinate withing an element*/
void global_coordinate_from_local(double *X, struct All_variables *E, double* epzl, int e){
	void global_coordinate_of_a_point();
	double **Xa;
	int i,nno;
//initiate variables
	Xa = (double**)malloc(9*sizeof(double*));
	for (i=0;i<=8;i++){
		Xa[i] = (double*)malloc(4*sizeof(double));
	}	
//get node coordinate
	for (i=1;i<=8;i++){
		nno = E->ien[1][e].node[i];
		Xa[i][1] = E->x[1][1][nno] ;
		Xa[i][2] = E->x[1][2][nno] ;
		Xa[i][3] = E->x[1][3][nno] ;
	}
	global_coordinate_of_a_point(X,E,Xa,epzl);
//froe space
	for (i=0;i<=8;i++){
		free(Xa[i]);
	}
	free(Xa);
}

/*get nature coordinate for a element by its global coordinate*/
void local_coordinate_from_global(double *epzl, struct All_variables *E, double *X, int e, int ndivide)
{
	void global_coordinate_from_local(double *X, struct All_variables *E, double* epzl, int e);
	int sphere_check_element(struct All_variables *E, double* X, double** Xn);
	double **Xa, **Xn;
	double *epzltemp;
	double epzlbound[5],epzlbound1[5];
	int i,nno,level,ein,ival,inode;
	double rad1,rad2,rad;
	//fprintf(stderr,"local_coordinate_from_global started\n"); //debug
//initialize
	Xa = (double**)malloc(9*sizeof(double*));
	for (i=0;i<=8;i++){
		Xa[i] = (double*)malloc(4*sizeof(double));
	}
	Xn = (double**)malloc(5*sizeof(double*));
	for (i=0;i<=4;i++){
		Xn[i] = (double*)malloc(4*sizeof(double));
	}
	epzltemp = (double*)malloc(4*sizeof(double));
//get nodes coordinates
	for (i=1;i<=8;i++){
		nno = E->ien[1][e].node[i];
		Xa[i][1] = E->x[1][1][nno] ;
		Xa[i][2] = E->x[1][2][nno] ;
		Xa[i][3] = E->x[1][3][nno] ;
		//fprintf(stderr,"inode = %d\nX = %.4e, Y = %.4e, Z = %.4e\n",i,Xa[i][1],Xa[i][2],Xa[i][3]); //debug
	}	
//get local coordinate on the r dimension
	nno = E->ien[1][e].node[1];
	rad1 = E->sx[1][3][nno];
	nno = E->ien[1][e].node[5];
	rad2 = E->sx[1][3][nno];
	rad = sqrt(X[1]*X[1]+X[2]*X[2]+X[3]*X[3]);
	epzl[3] = (-1.0)*(rad-rad2)/(rad1-rad2)+(rad-rad1)/(rad2-rad1);
	//fprintf(stderr,"rad1 = %.8e, rad2 = %.8e, rad = %.8e\n",rad1, rad2, rad);
	//fprintf(stderr,"epzl3: %.4e\n",epzl[3]); //debug
//get local coordinate on the sphere
	//assign intial value for epzlbound
	epzlbound[1] = -1.0;
	epzlbound[2] = 1.0;
	epzlbound[3] = -1.0;
	epzlbound[4] = 1.0;
	for (level=1;level<=ndivide;level++){
		//loop for four inside elements
		//fprintf(stderr,"level = %d\n",level);//debug
		for (ein=1;ein<=4;ein++){
			//get nature coordinate bound for this inside element
			//fprintf(stderr,"ein = %d\n",ein); //stderr
			if ((ein==1)||(ein==4)){
				epzlbound1[1] = epzlbound[1];
				epzlbound1[2] = (epzlbound[1]+epzlbound[2])/2.0;
			}
			else{
				epzlbound1[1] = (epzlbound[1]+epzlbound[2])/2.0;
				epzlbound1[2] = epzlbound[2];
			}
			if ((ein==1)||(ein==2)){
				epzlbound1[3] = epzlbound[3];
				epzlbound1[4] = (epzlbound[3]+epzlbound[4])/2.0;
			}
			else{
				epzlbound1[3] = (epzlbound[3]+epzlbound[4])/2.0;
				epzlbound1[4] = epzlbound[4];
			}
			//get global coordinate for nodes of this inside element
			for (inode=1;inode<=4;inode++){
				if((inode==1)||(inode==4)){
					epzltemp[1] = epzlbound1[1];
				}
				else{
					epzltemp[1] = epzlbound1[2];	
				}
				if((inode==1)||(inode==2)){
					epzltemp[2] = epzlbound1[3];
				}
				else{
					epzltemp[2] = epzlbound1[4];	
				}
				epzltemp[3] = epzl[3];
				global_coordinate_from_local(Xn[inode], E, epzltemp, e);
				//fprintf(stderr,"inode = %d\nX = %.4e, Y = %.4e, Z = %.4e\n",inode,Xn[inode][1],Xn[inode][2],Xn[inode][3]); //debug
			}
			//determine if points in this inside element
			ival = sphere_check_element(E,X,Xn);
			/*
			if (ival){
				fprintf(stderr,"Yes\n");
			}
			else{
				fprintf(stderr,"No\n");
			}*///debug
			//get natrue boundary for next step
			if (ival){
				for (i=1;i<=4;i++){
					epzlbound[i] = epzlbound1[i];
				}
				break;
			}
			else if((!ival)&&ein==4){
				fprintf(stderr,"something wrong while get global coordinate::\n");
				fprintf(stderr,"points: X = %.4e, Y = %.4e, Z = %.4e\n", X[1],X[2],X[3]);
				fprintf(stderr,"element: %d\n", e);
				exit(10);
			}
		}
	}
//assign average of boundary to epzl1 and epzl2
	epzl[1] = (epzlbound[1]+epzlbound[2])/2.0;
	epzl[2] = (epzlbound[3]+epzlbound[4])/2.0;
//free space
	free(epzltemp);
	for (i=0;i<=8;i++){
		free(Xa[i]);
	}
	free(Xa);	
	for (i=0;i<=4;i++){
		free(Xn[i]);
	}
	free(Xn);	
	//fprintf(stderr,"local_coordinate_from_global ended\n"); //debug
}

/*check if a point of given cartisen coordinate with element of number e*/
int check_points_in_element(struct All_variables *E, double* X, int e)
{
	int sphere_check_element();
	double **Xa, rad1, rad2, rad;
	int i,nno,iath,ival,irad;
	double tiny = 1e-10;
//initiate variables
	Xa = (double**)malloc(9*sizeof(double*));
	for (i=0;i<=8;i++){
		Xa[i] = (double*)malloc(4*sizeof(double));
	}
//get node coordinates
//	fprintf(stderr,"ele = %d\n",e); //debug
	for (i=1;i<=8;i++){
		nno = E->ien[1][e].node[i];
		Xa[i][1] = E->x[1][1][nno] ;
		Xa[i][2] = E->x[1][2][nno] ;
		Xa[i][3] = E->x[1][3][nno] ;
/*
		fprintf(stderr,"node %d:\n",i);
		fprintf(stderr,"x = %.4e, y = %.4e, z = %.4e\n",Xa[i][1],Xa[i][2],Xa[i][3]); //debug
*/
	}
//get two radius of upper and lower boundary of node
	nno = E->ien[1][e].node[1];
	rad1 = E->sx[1][3][nno]; 
	nno = E->ien[1][e].node[5];
	rad2 = E->sx[1][3][nno];
	//fprintf(stderr,"rad1 = %.4e, rad2 = %.4e\n",rad1,rad2); //debug
//check if withing this range
	rad = sqrt(X[1]*X[1]+X[2]*X[2]+X[3]*X[3]);
	if ((rad>(rad1-tiny))&&(rad<(rad2+tiny))){
		irad = 1;
	}
	else{
		irad = 0;
	}	
	//fprintf(stderr,"irad = %d\n",irad); //debug
//check if points is in the element
	iath = sphere_check_element(E, X, Xa);
	ival = (iath&&irad);
/*
	if (ival){
		fprintf(stderr,"Yes\n");
	}
	else{
		fprintf(stderr,"No\n");
	}
	//debug
*/

//free space
	for (i=0;i<=8;i++){
		free(Xa[i]);
	}
	free(Xa);
	return ival;
}
/*get globale coordinate of a point with given nature coordinate in a given element of 8 nodes*/
void global_coordinate_of_a_point(double *Xc,struct All_variables *E,double **Xa,double *epzl){
	/* Xc: return value of global coordinate
	 * Xa: global coordinates of nodes
	 * epzl: nature coordiantes of point*/
	double lpoly();
	const int enodes = 8;
	const int nsd = 3;
	int i,d;
	double vpt[9];
	for(i=1; i<=enodes; i++){
		vpt[i] = 1.0;
		for(d=1; d<=nsd; d++){
			vpt[i] *= lpoly(bb[d-1][i],epzl[d]);
		}
	}
	for(d=1; d<=nsd; d++){
			Xc[d] = 0.0;
		for(i=1; i<=enodes; i++){
			Xc[d] += vpt[i]*Xa[i][d];
		}
	}
}

/*check if the profject of a point on the sphere is in a sphere element*/
int sphere_check_element(struct All_variables *E, double* X, double** Xn){
	/*X is a 3-vector of point coordinate
	 * Xn is a 4*3-vector of node coordiate*/
	void full_coordinate_information_of_a_point();
	int icheck_bounds();
	double **rnode,*point;
	double rad;
	int i,ival;
//initialize
	rnode = (double**)malloc(5*sizeof(double*));
	for (i=0;i<=4;i++){
		rnode[i] = (double*)malloc(10*sizeof(double));
	}
	point = (double*)malloc(4*sizeof(double));
//get input vectors for icheck_bounds
	for (i=1;i<=4;i++){
		full_coordinate_information_of_a_point(rnode[i], E, Xn[i]);
		//fprintf(stderr,"rnode %d\n",i); //debug
	}
	rad = sqrt(X[1]*X[1]+X[2]*X[2]+X[3]*X[3]);
	point[1] = X[1]/rad;
	point[2] = X[2]/rad;
	point[3] = X[3]/rad;
//	fprintf(stderr,"points: %.4e, %.4e, %.4e\n",point[1],point[2],point[3]);//debug
//determine if is in the element
	ival=icheck_bounds(E,point,rnode[1],rnode[2],rnode[3],rnode[4]);
	//fprintf(stderr,"ival = %d\n",ival);//debug
//free space
	for (i=0;i<=4;i++){
		free(rnode[i]);
		}
	free(point);
	free(rnode);

	return ival;	
} 
/*this is to generate input vector for icheck_bound_input()*/
void full_coordinate_information_of_a_point(double *XX,struct All_variables *E, double *X)
{
	//fprintf(stderr,"full_coordinate_information_of_a_point started\n"); //debug
	void cart_to_sphere();
	double temp;
	XX[1] = X[1];
	XX[2] = X[2];
	XX[3] = X[3];
	cart_to_sphere(E,X[1],X[2],X[3],&XX[4],&XX[5],&temp);
	temp = sqrt(X[1]*X[1]+X[2]*X[2]);
	XX[6] = cos(XX[4]); 
	XX[7] = sin(XX[4]); 
	XX[8] = cos(XX[5]);
	XX[9] = sin(XX[5]);	
	//fprintf(stderr,"full_coordinate_information_of_a_point ended\n"); //debug
}
