#include <R.h>
#include <Rmath.h>
#include <stdlib.h> 
#include <Rdefines.h>
#include <time.h>


double gauss(double, double );
double gauss2(double, double, double);
double Fgauss(double, double );
double max(double, double );
void q_sort(double[], int, int);
void calcSum(double[], double[], double[], int, int, double);


/* calcPDF */

void calcPDF(double *res, double *r, int *len, 
         int *layer_start, int *layer_end, int *layers_num,
         double *np, int *nrow, double *scale, 
	     double *a1, double *b1, 
		 double *a2, double *b2, 
		 double *a3, double *b3, 
		 double *a4, double *b4,
	     double *a5, double *b5,
		 double *c,  double *scatterLengths, 
		 int *type, double *Qmin,
	     double *dr, double *minR) {

  int i, j, k, l, p;

  double ** nanop = (double **)R_alloc( *nrow, sizeof( double * ));
  for( int i = 0; i < *nrow; i++ )
    nanop[i] = (double *)R_alloc(3, sizeof( double ));
  int * count =  (int *)R_alloc( *len, sizeof( int ));


  double ff1, ff2;
  double q4;
  double dist;

  j = 0;
  for (i=0; i < *nrow; i++) {
    nanop[i][0] = np[j];
    nanop[i][1] = np[j+1];
    nanop[i][2] = np[j+2];
    j+=3;
  }

  q4 = pow(*Qmin/(4*M_PI),2);
  if(*type) { // x-ray scattering 
    ff1 = a1[0]*exp(-b1[0]*q4) + a2[0]*exp(-b2[0]*q4) +
	       a3[0]*exp(-b3[0]*q4) + a4[0]*exp(-b4[0]*q4) + 
		   a5[0]*exp(-b5[0]*q4)+c[0]; 
    ff2 = a1[1]*exp(-b1[1]*q4) + a2[1]*exp(-b2[1]*q4) +
	       a3[1]*exp(-b3[1]*q4) + a4[1]*exp(-b4[1]*q4) + 
		   a5[1]*exp(-b5[1]*q4)+c[1]; 
  }
  else { // neutron scattering 
    ff1 = scatterLengths[0];
	ff2 = scatterLengths[1];
  }
  

  for (l=0; l<*layers_num; l++){
    for (i=1; i < *len; i++)
      count[i] = 0;
    	
    for (j=layer_start[l]-1; j < layer_end[l]; j++) {
      for (k=layer_start[l]-1; k < j; k++) { 
        dist = sqrt(pow(nanop[j][0]-nanop[k][0],2)+
                pow(nanop[j][1]-nanop[k][1],2)+
		        pow(nanop[j][2]-nanop[k][2],2));
        i = round( (dist-*minR)/ *dr);
        if(i >= 0 && i < *len)
	      count[i] += 1;
      }
    }
	
//	p = layer_end[l] - layer_start[l] + 1;
	for (i=1; i < *len; i++)
	  res[i] += count[i] * scale[l]; //nrow_tot[l] / pow(p*ff1 + (nrow_tot[l]-p)*ff2, 2); 
  }
   
  res[0] = 0;  //prevent div by 0
  for (i=1; i < *len; i++) {
    res[i] = 2*ff1*ff1* res[i] / (4.0*M_PI*pow(r[i],2)* *dr)/ *layers_num;
  }

}

/* calcPDF_CS */

void calcPDF_CS(int *layer_start, int *layer_end, 
        int *layerS_start, int *layerS_end, int *layers_num, 
        double *res, double *r, int *len, double *scale, 
        double *np_mu, int *nrow_mu, double *np_nu, int *nrow_nu,
	    double *a1, double *b1, 
		double *a2, double *b2, 
		double *a3, double *b3, 
		double *a4, double *b4,
	    double *a5, double *b5,
		double *c,  double *scatterLengths, 
		int *type, double *Qmin,
		double *dr, double *minR) {

  int i, j, k, l, p1, p2;
  double ** nanop_mu = (double **)R_alloc( *nrow_mu, sizeof( double * ));
  for( int i = 0; i < *nrow_mu; i++ )
    nanop_mu[i] = (double *)R_alloc(3, sizeof( double ));
  double ** nanop_nu = (double **)R_alloc( *nrow_nu, sizeof( double * ));
  for( int i = 0; i < *nrow_nu; i++ )
    nanop_nu[i] = (double *)R_alloc(3, sizeof( double ));
  int * count =  (int *)R_alloc( *len, sizeof( int ));
  
  double ff1, ff2;
  double q4;

  double dist;
  j = 0;
  for (i=0; i < *nrow_mu; i++) {
    nanop_mu[i][0] = np_mu[j];
    nanop_mu[i][1] = np_mu[j+1];
    nanop_mu[i][2] = np_mu[j+2];
    j+=3;
  }
  j = 0;
  for (i=0; i < *nrow_nu; i++) {
    nanop_nu[i][0] = np_nu[j];
    nanop_nu[i][1] = np_nu[j+1];
    nanop_nu[i][2] = np_nu[j+2];
    j+=3;
  }

  if(*type) { // x-ray scattering 
    q4 = pow(*Qmin/(4*M_PI),2);
    ff1 = a1[0]*exp(-b1[0]*q4) + a2[0]*exp(-b2[0]*q4) +
	       a3[0]*exp(-b3[0]*q4) + a4[0]*exp(-b4[0]*q4) + 
		   a5[0]*exp(-b5[0]*q4)+c[0]; 
    ff2 = a1[1]*exp(-b1[1]*q4) + a2[1]*exp(-b2[1]*q4) +
	       a3[1]*exp(-b3[1]*q4) + a4[1]*exp(-b4[1]*q4) + 
		   a5[1]*exp(-b5[1]*q4)+c[1];   
  }
  else { // neutron scattering 
    ff1 = scatterLengths[0]; 
	ff2 = scatterLengths[1];
  }

  for (l=0; l<*layers_num; l++){
    for (i=1; i < *len; i++)
      count[i] = 0;
    	
    for (j=layerS_start[l]-1; j < layerS_end[l]; j++) {
      for (k=layer_start[l]-1; k < layer_end[l]; k++) {
        dist = sqrt(pow(nanop_nu[j][0]-nanop_mu[k][0],2)+
	   	    pow(nanop_nu[j][1]-nanop_mu[k][1],2)+
	   	    pow(nanop_nu[j][2]-nanop_mu[k][2],2));	
        
		i = round( (dist-*minR)/ *dr);
        if(i >= 0 && i < *len)
		  count[i] += 1;
      }
    }
//    p1 = layer_end[l] - layer_start[l] + 1;
//    p2 = layerS_end[l] - layerS_start[l] + 1;
    for (i=0; i<*len; i++)
	  res[i] += count[i] * scale[l]; //(p1+p2) / pow(p1*ff1 + p2*ff2, 2);
  
  }
  
  // mult by two since do not loop by matrix twice
  res[0] = 0;  //prevent div by 0

  for (i=1; i < *len; i++) {
    res[i] = 2 * ff1*ff2 *res[i] / (4.0*M_PI*pow(r[i],2)* *dr ) / *layers_num;
  }  
  
}

/* calcSc */

void calcSc(double *res, double *Q, int *len,  int *nrow,
	    int *atomType, int *natomTypes,
		    double *a1,
		    double *b1, double *a2, double *b2, double *a3,
		    double *b3, double *a4, double *b4,
		    double *a5, double *b5, double *c,
		    double *scatterLengths, int *type) {

  int i, j, k;

  double ** ffM =  (double **)R_alloc( *len, sizeof( double * ));
  for( int i = 0; i < *len; i++ )
    ffM[i] = (double *)R_alloc(*natomTypes, sizeof( double ));

  double ff, ff1, ff2, ffAv;

  double ntmp, q4;

  if(*type) { // x-ray scattering
    Rprintf("Calculating X-ray scattering\n");
    for (i=0; i < *len; i++) {
      q4 = Q[i]/(4*M_PI);
      for (j=0; j < *natomTypes; j++) {
	ffM[i][j] = (a1[atomType[j]])*exp(-(b1[atomType[j]])*q4)+
	  (a2[atomType[j]])*exp(-(b2[atomType[j]])*q4)+
	  (a3[atomType[j]])*exp(-(b3[atomType[j]])*q4)+
	  (a4[atomType[j]])*exp(-(b4[atomType[j]])*q4)+
	  (a5[atomType[j]])*exp(-(b5[atomType[j]])*q4)+(c[atomType[j]]);

      }
    }
  }
  else { // neutron scattering
    Rprintf("Calculating neutron scattering\n");
    q4 = Q[0]/(4*M_PI);
    for (i=0; i < *len; i++) {
      for (j=0; j < *natomTypes; j++) {
	ffM[i][j] = scatterLengths[atomType[j]];
      }
    }
  }

  for (i=0; i < *len; i++) {
    ntmp = 0;
    ffAv = 0;
    for (j=0; j < *nrow; j++) {
      for (k=0; k < *nrow; k++) {

	q4 = Q[i]/(4*M_PI);

	ff1 = ffM[i][atomType[j]];
	ff2 = ffM[i][atomType[k]];

	ff = ff2*ff1;
	if(j==0)
	  ffAv += ff2;

	ntmp = ntmp + ff;
      }
    }

    res[i] =  1/((double)*nrow*pow( (ffAv/(double)*nrow),2))* ntmp;
  }
}

/* calcSc_CS */
void calcSc_CS(double *res,
		       double *Q, int *len, double *minQ,
		       double *dQ, double *npC, double *npS,
		       int *nrowC, int *nrowS,
		       int *atomTypeC, int *atomTypeS,
		       int *natomTypesC, int *natomTypesS,
		       double *a1,
		       double *b1, double *a2, double *b2, double *a3,
		       double *b3, double *a4, double *b4,
		       double *a5, double *b5,double *c,
		       double *scatterLengths, int *type) {

  int i, j, k;

  double ** ffMC =  (double **)R_alloc( *len, sizeof( double * ));
  for( int i = 0; i < *len; i++ )
    ffMC[i] = (double *)R_alloc(*natomTypesC, sizeof( double ));
  double ** ffMS =  (double **)R_alloc( *len, sizeof( double * ));
  for( int i = 0; i < *len; i++ )
    ffMS[i] = (double *)R_alloc(*natomTypesS, sizeof( double ));


  double ff, ff1, ff2;
  double ffAvCS,ffAvC,ffAvS = 0;

  double ntmpC,ntmpS,ntmpCS, q4;

  if(*type) { // x-ray scattering
    Rprintf("Calculating X-ray scattering\n");
    for (i=0; i < *len; i++) {
      q4 = Q[i]/(4*M_PI);
      for (j=0; j < *natomTypesC; j++) {
	ffMC[i][j] = (a1[atomTypeC[j]])*exp(-(b1[atomTypeC[j]])*q4)+
	  (a2[atomTypeC[j]])*exp(-(b2[atomTypeC[j]])*q4)+
	  (a3[atomTypeC[j]])*exp(-(b3[atomTypeC[j]])*q4)+
	  (a4[atomTypeC[j]])*exp(-(b4[atomTypeC[j]])*q4)+
	  (a5[atomTypeC[j]])*exp(-(b5[atomTypeC[j]])*q4)+(c[atomTypeC[j]]);
      }
      for (j=0; j < *natomTypesS; j++) {
	  ffMS[i][j] = (a1[atomTypeS[j]])*exp(-(b1[atomTypeS[j]])*q4)+
	    (a2[atomTypeS[j]])*exp(-(b2[atomTypeS[j]])*q4)+
	    (a3[atomTypeS[j]])*exp(-(b3[atomTypeS[j]])*q4)+
	    (a4[atomTypeS[j]])*exp(-(b4[atomTypeS[j]])*q4)+
	    (a5[atomTypeS[j]])*exp(-(b5[atomTypeS[j]])*q4)+(c[atomTypeS[j]]);
      }
    }
  }
  else { // neutron scattering
    Rprintf("Calculating neutron scattering\n");
    q4 = Q[0]/(4*M_PI);
    for (i=0; i < *len; i++) {
      for (j=0; j < *natomTypesC; j++) {
	ffMC[i][j] = scatterLengths[atomTypeC[j]];
      }
      for (j=0; j < *natomTypesS; j++) {
	ffMS[i][j] = scatterLengths[atomTypeS[j]];
      }
    }
  }
  ffAvCS = 0;
  for (j=0; j < *nrowC; j++) {
    ffAvCS += scatterLengths[atomTypeC[j]];
  }
  for (j=0; j < *nrowS; j++) {
    ffAvCS += scatterLengths[atomTypeS[j]];
  }
  ffAvCS = ffAvCS/(double)(*nrowC+*nrowS);

  for (i=0; i < *len; i++) {
    ntmpC = ntmpS = ntmpCS = 0;
    ffAvC=ffAvS=0;
    for (j=0; j < *nrowC; j++) {
      for (k=0; k < *nrowC; k++) {

	ff1 = ffMC[i][atomTypeC[j]];
	ff2 = ffMC[i][atomTypeC[k]];

	ff = ff2*ff1;

	if(j==0)
	  ffAvC += ff2;
	ntmpC = ntmpC + ff;

      }
    }
    for (j=0; j < *nrowS; j++) {
      for (k=0; k < *nrowS; k++) {

	ff1 = ffMS[i][atomTypeS[j]];
	ff2 = ffMS[i][atomTypeS[k]];

	ff = ff2*ff1;

	if(j==0)
	  ffAvS += ff2;
	ntmpS = ntmpS + ff;
      }
    }

    for (j=0; j < *nrowC; j++) {
      for (k=0; k < *nrowS; k++) {

	ff1 = ffMC[i][atomTypeC[j]];
	ff2 = ffMS[i][atomTypeS[k]];

	ff = ff2*ff1;

	ntmpCS = ntmpCS + ff;
      }
    }

    Rprintf("** %f *  %f \n", ffAvCS, (ffAvC + ffAvS)/  (double)(*nrowC+*nrowS) );
    //ntmpC =   (1/((double)(*nrowC) *pow( (ffAvC/(double)*nrowC),2)))* ntmpC *(double)(*nrowC)/(double)(*nrowC+*nrowS);
    //ntmpS =   (1/((double)(*nrowS) *pow( (ffAvS/(double)*nrowS),2)))* ntmpS *(double)(*nrowS)/(double)(*nrowC+*nrowS);
    //ntmpCS =  (1/((double)(*nrowC+*nrowS)*pow(ffAvCS,2)))*ntmpCS*2;

    ntmpC =   (1/((double)(*nrowC) * pow(ffAvCS,2)))* ntmpC *(double)(*nrowC)/(double)(*nrowC+*nrowS);
    ntmpS =   (1/((double)(*nrowS) * pow(ffAvCS,2)))* ntmpS *(double)(*nrowS)/(double)(*nrowC+*nrowS);
    ntmpCS =  (1/((double)(*nrowC+*nrowS)*pow(ffAvCS,2)))*ntmpCS*2;

    res[i] =  ntmpC + ntmpS + ntmpCS;
  }
  Rprintf("S: %d C: %d \n", *nrowS, *nrowC);
}

/* calcQDepPDF */

void calcQDepPDF(double *res,
		 double *Q, double *r,
		 int *len, 
		 double *totalScattVec) {

  int i, j, k;
    
  for (i=0; i < *len; i++) 
    res[i] = totalScattVec[i];

  for (i=0; i < *len; i++) 
    res[i] = res[i] * sin(Q[i]* *r) * Q[i];

}

/* calcQDepPDF */

void calcQDepPDF_smoothSQ(double *res,
			  double *Q, double *r,
			  int *len, double *np, int *nrow,
			  double *a1, double *b1,
			  double *a2, double *b2,
			  double *a3, double *b3,
			  double *a4, double *b4,
			  double *a5, double *b5,
			  double *c, int *preComp,
		 double *totalScattVec) {


  int i, j, k;
  double *ff =  (double *)R_alloc(*len, sizeof( double ));
  double q4;
  for (i=0; i < *len; i++) {
    q4 = Q[i]/(4*M_PI);
    ff[i] = (*a1)*exp(-(*b1)*q4)+
      (*a2)*exp(-(*b2)*q4)+
      (*a3)*exp(-(*b3)*q4)+
      (*a4)*exp(-(*b4)*q4)+
      (*a5)*exp(-(*b5)*q4)+(*c);
    ff[i] = pow(ff[i],2);
  }


  if(*preComp == 0) {
    double ** nanop = (double **)R_alloc( *nrow, sizeof( double * ));
    for( int i = 0; i < *nrow; i++ )
      nanop[i] = (double *)R_alloc(3, sizeof( double ));
    double ** dist =  (double **)R_alloc( *nrow, sizeof( double * ));
    for( int i = 0; i < *nrow; i++ )
      dist[i] = (double *)R_alloc(*nrow, sizeof( double ));

    double ntmp;

    j = 0;
    for (i=0; i < *nrow; i++) {
      nanop[i][0] = np[j];
      nanop[i][1] = np[j+1];
      nanop[i][2] = np[j+2];
      j+=3;
    }

    for (j=0; j < *nrow; j++) {
      for (k=0; k < *nrow; k++) {
	dist[j][k] = sqrt(pow(nanop[j][0]-nanop[k][0],2)+
			  pow(nanop[j][1]-nanop[k][1],2)+
			  pow(nanop[j][2]-nanop[k][2],2));

      }
    }
    for (i=0; i < *len; i++) {
      for (j=0; j < *nrow; j++) {
	for (k=0; k < *nrow; k++) {
	  if(dist[j][k]!=0)
	    res[i] = res[i] + (ff[i] *(sin(Q[i]*dist[j][k]))/(Q[i]*dist[j][k]));
	}
      }
    }
  }
  else {
    for (i=0; i < *len; i++) {
      res[i] = totalScattVec[i];
    }
  }
  for (i=0; i < *len; i++) {
    if(*preComp == 0)
      res[i] = 1/((double)*nrow*ff[i]) * res[i] * sin(Q[i]* *r) * Q[i];
    else
      res[i] = res[i] * sin(Q[i]* *r); // assume res represents
                                       // smoothed [S(Q)-1] Q
  }
}
/* calcRedTotalScatt */
void calcRedTotalScatt(double *res, double *Q, int *len, double *minQ,
		       double *dQ, double *np, int *nrow, double *a1,
		       double *b1, double *a2, double *b2, double *a3,
		       double *b3, double *a4, double *b4,
		       double *a5, double *b5,double *c)
		       {

  int i, j, k;
  double ** nanop = (double **)R_alloc( *nrow, sizeof( double * ));
  for( int i = 0; i < *nrow; i++ )
    nanop[i] = (double *)R_alloc(3, sizeof( double ));
  double ** dist =  (double **)R_alloc( *nrow, sizeof( double * ));
  for( int i = 0; i < *nrow; i++ )
    dist[i] = (double *)R_alloc(*nrow, sizeof( double ));

  double *ff =  (double *)R_alloc(*len, sizeof( double ));

  double ntmp, q4;


  j = 0;
  for (i=0; i < *nrow; i++) {
    nanop[i][0] = np[j];
    nanop[i][1] = np[j+1];
    nanop[i][2] = np[j+2];
    j+=3;
  }
  for (i=0; i < *len; i++) {
    q4 = Q[i]/(4*M_PI);
    ff[i] = (*a1)*exp(-(*b1)*q4)+
      (*a2)*exp(-(*b2)*q4)+
      (*a3)*exp(-(*b3)*q4)+
      (*a4)*exp(-(*b4)*q4)+
      (*a5)*exp(-(*b5)*q4)+(*c);

    ff[i] = pow(ff[i],2);
  }
  for (j=0; j < *nrow; j++) {
    for (k=0; k < *nrow; k++) {
      dist[j][k] = sqrt(pow(nanop[j][0]-nanop[k][0],2)+
			pow(nanop[j][1]-nanop[k][1],2)+
			pow(nanop[j][2]-nanop[k][2],2));
    }
  }
  for (i=0; i < *len; i++) {
    for (j=0; j < *nrow; j++) {
      for (k=0; k < *nrow; k++) {
	if(dist[j][k]!=0)
	  res[i] = res[i] + (ff[i] * (sin(Q[i]*dist[j][k]))/(dist[j][k]));
      }
    }
    res[i] = 1/((double)*nrow*ff[i] * res[i]);
  }

}


/* simPart */

void simPart(double *res, 
	     double *a1, 
	     double *a2, 
	     double *a3,  
	     int *a, int *b, int *c) {

  double sh;
  int i, j, k, h;
  int cnt = 0;
  
  for (i=-*a; i < *a+1; i++) { 
     for (j=-*b; j < *b+1; j++) { 
       for (k=-*c; k < *c+1; k++) { 	    
		 res[cnt] =   ((a1[0]*(double)i) + (a2[0]*(double)j) + (a3[0]*(double)k));
		 res[cnt+1] = ((a1[1]*(double)i) + (a2[1]*(double)j) + (a3[1]*(double)k));
		 res[cnt+2] = ((a1[2]*(double)i) + (a2[2]*(double)j) + (a3[2]*(double)k));
         cnt +=3;  
       }
     }
  }
}

/* **************************************
  CUBIC PACKING WITH STACKING FAULTS
***************************************** */
void simPartStackCub(double *res, int *stacks, double *a, int *Nx, int *Ny, int *Nz, int *nStacks) 
{
  double  b, by, bz;
  int i, j, k, plane=0;
  int cnt = 0;

  b = *a/sqrt(2);
  by = *a*sqrt(3)/(2*sqrt(2));
  bz = *a*sqrt(3)/3;
  *nStacks = 0;
  for (i=-*Nz; i < *Nz+1; i++) {  //z
	if (i == stacks[*nStacks]){
	  *nStacks += 1;
	  plane += 1;
	  plane = plane%3;
	}
    for (j=-*Ny; j < *Ny+1; j++) {  //t
      for (k=-*Nx; k < *Nx+1; k++) {  //x
	    res[cnt] = ((double)k + 0.5*(double)j + 0.5*(double)plane) * b;  //x
		res[cnt+1] =((double)j + (double)plane/3) * by;  //t
		res[cnt+2] =  (double)i* bz ;  //z
	    cnt +=3;
       }
    }
	plane += 1; 
    plane = plane%3;
  }
  
}


/* **************************************
  HEXAGONAL PACKING WITH STACKING FAULTS
***************************************** */
void simPartStackHex(double *res, int *stacks, double *a, double *c, int *Nx, int *Ny, int *Nz, int *nStacks) 
{
  double by;
  int i, j, j2, k;
  int Lp = 0, Lc = 1, Lf = 2, Ltmp;
  int cnt = 0;

  by = *a*sqrt(3)/2;
  *nStacks = 0;
  

  for (i=-*Nz; i < *Nz+1; i++) {  //z
 	if (i == stacks[*nStacks]){
	  *nStacks += 1;
	  Ltmp = Lc;
      Lc = Lf;
      Lf = Lp;	  
      Lp = Ltmp;
    }
	else{
	  Ltmp = Lp;
	  Lp = Lc;
	  Lc = Ltmp;	
	}
	
    for (j=-*Ny; j < *Ny+1; j++) {  //t
	  j2 = j%2;
      for (k=-*Nx; k < *Nx+1; k++) {  //x
	    res[cnt] = ( (double)k + 0.5*((double)j2 + (double)Lc) )* *a;  //x
		res[cnt+1] =((double)j + (double)Lc/3)*by;  //t
		res[cnt+2] =  (double)i* *c/2 ;  //z
	    cnt +=3;
      }
    }
  }  
}









/* termRip */

void termRip(double *res, double *pdf, int *len, double *qmax,
	     double *deltar, double *rmax, int *lenRip) {

  int i, j;

  double zpi = 8.00*atan(1.00);
  double sincut = 0.025;
  double rcut   = 1.0 / (*qmax*sincut);
  double z      = *deltar* *qmax;

  int pdf_bin  = round((*rmax+rcut) / *deltar);
  double *pdf_sinc = (double *)R_alloc(2*pdf_bin+2, sizeof( double ));
  double *ppp = (double *)R_alloc(pdf_bin, sizeof( double ));
  double *pdf_calc = (double *)R_alloc(pdf_bin, sizeof( double ));
  int n = round(rcut / *deltar) + 100;

  pdf_sinc[0] = 0.0;
  for (i=1; i < (2*pdf_bin+2); i++)
    if(i < n)
      pdf_sinc[i] = (sin(z*(double)i))/(*deltar*(double)i);
    else
      pdf_sinc[i] = 0.0;

  for (i=0; i < (pdf_bin); i++) {
    ppp[i] = 0;
    if(i < (*len))
      pdf_calc[i] = pdf[i];
    else
      pdf_calc[i] = 0.0;
  }

  for (i=1; i < *lenRip; i++) {
    ppp[i] = pdf_calc[i]*(*qmax -pdf_sinc[2*i]);
    for (j=1; j <= (i-1); j++)
      ppp[i] = ppp[i] + pdf_calc[j]*(pdf_sinc[i-j]-pdf_sinc[i+j]);
    for (j=(i+1); j < pdf_bin; j++)
      ppp[i] = ppp[i] + pdf_calc[j]*(pdf_sinc[j-i]-pdf_sinc[j+i]);
  }
  for (i=0; i < *lenRip; i++)
    pdf_calc[i] = ppp[i]* *deltar/zpi*2.0;

  for (i=0; i < *lenRip; i++)
    res[i] = pdf_calc[i];
  for (i=*lenRip; i < (*len); i++)
    res[i] = pdf[i];
	
	
	

}


/* gaussConvol */

void gaussConvol(double *SQ, double *Q, int *M, int *len, double *Qdamp,  double *dQ){

  int i, j;
  double *I2 = (double *)R_alloc(*len, sizeof( double ));

  for (i=*M; i < (*len-*M); i++){
    I2[i] = 0;
	for(j=-*M; j < (*M+1); j++)
	  I2[i] = I2[i] + SQ[i+j]* *dQ*1/sqrt(2*3.1415926* pow(*Qdamp,2) )* exp(- pow( Q[i+j]-Q[i],2)/(2* pow(*Qdamp,2)));
	I2[i] = I2[i] - 1;
  }
  for(i=0; i < *M; i++)
    I2[i] = SQ[i]-1;
  for(i=(*len-*M); i<*len; i++ )
    I2[i] = SQ[i]-1	;

  for(i=0; i < *len; i++)
    SQ[i] = I2[i];
}



/* calcIn */
void calcIn(double *res, double *Q, int *len, double *minQ,
		    double *dQ, double *np, int *nrow,
		    int *atomType, int *natomTypes,
		    double *a1,
		    double *b1, double *a2, double *b2, double *a3,
		    double *b3, double *a4, double *b4,
		    double *a5, double *b5, double *c,
		    double *scatterLengths, int *type,
		    double *sigma,
		    int *useN, double *n, double *delta) {

  int i, j, k;

  double ** nanop = (double **)R_alloc( *nrow, sizeof( double * ));
  for( int i = 0; i < *nrow; i++ )
    nanop[i] = (double *)R_alloc(3, sizeof( double ));
  double ** dist =  (double **)R_alloc( *nrow, sizeof( double * ));
  for( int i = 0; i < *nrow; i++ )
    dist[i] = (double *)R_alloc(*nrow, sizeof( double ));
  double ** ffM =  (double **)R_alloc( *len, sizeof( double * ));
  for( int i = 0; i < *len; i++ )
    ffM[i] = (double *)R_alloc(*natomTypes, sizeof( double ));


  double ff, ff1, ff2;
  double * ffAv = (double *)R_alloc(*len, sizeof( double ));
  double * fA = (double *)R_alloc(*len, sizeof( double ));

  double ntmp, q4;

  j = 0;
  for (i=0; i < *nrow; i++) {
    nanop[i][0] = np[j];
    nanop[i][1] = np[j+1];
    nanop[i][2] = np[j+2];
    j+=3;
  }

  for (j=0; j < *nrow; j++) {
    for (k=0; k < *nrow; k++) {
      dist[j][k] = sqrt(pow(nanop[j][0]-nanop[k][0],2)+
			pow(nanop[j][1]-nanop[k][1],2)+
			pow(nanop[j][2]-nanop[k][2],2));

    }
  }
  if(*type) { // x-ray scattering
    Rprintf("Calculating X-ray scattering\n");
    for (i=0; i < *len; i++) {
      q4 = Q[i]/(4*M_PI);
      for (j=0; j < *natomTypes; j++) {
	ffM[i][j] = (a1[atomType[j]])*exp(-(b1[atomType[j]])*q4)+
	  (a2[atomType[j]])*exp(-(b2[atomType[j]])*q4)+
	  (a3[atomType[j]])*exp(-(b3[atomType[j]])*q4)+
	  (a4[atomType[j]])*exp(-(b4[atomType[j]])*q4)+
	  (a5[atomType[j]])*exp(-(b5[atomType[j]])*q4)+(c[atomType[j]]);

      }
    }
  }
  else { // neutron scattering
    Rprintf("Calculating neutron scattering\n");
    q4 = Q[0]/(4*M_PI);
    for (i=0; i < *len; i++) {
      for (j=0; j < *natomTypes; j++) {
	ffM[i][j] = scatterLengths[atomType[j]];
      }
    }
  }
  for (i=0; i < *len; i++) {
    for (j=0; j < *nrow; j++) {
      for (k=0; k < *nrow; k++) {

	q4 = Q[i]/(4*M_PI);

	ff1 = ffM[i][atomType[j]];
	ff2 = ffM[i][atomType[k]];

	ffAv[i] += ff2*ff1;

	if(j==0)
	  fA[i] += ff1;
      }
    }
    ffAv[i] = ffAv[i]/(pow((double)*nrow,2)) ;

    fA[i] = fA[i]/(double)*nrow;
  }


  for (i=0; i < *len; i++) {
    ntmp = 0;
    for (j=0; j < *nrow; j++) {
      for (k=0; k < *nrow; k++) {

	q4 = Q[i]/(4*M_PI);

	ff1 = ffM[i][atomType[j]];
	ff2 = ffM[i][atomType[k]];

	ff = ff2*ff1;

	ntmp += (1+res[i])*((double)*nrow/ffAv[i]) + pow((ff1 - fA[i]),2);
      }

    }

    res[i] = ntmp;
  }
}

/* calcIn_CS */
void calcIn_CS(double *res,
		       double *Q, int *len, double *minQ,
		       double *dQ, double *npC, double *npS,
		       int *nrowC, int *nrowS,
		       int *atomTypeC, int *atomTypeS,
		       int *natomTypesC, int *natomTypesS,
		       double *a1,
		       double *b1, double *a2, double *b2, double *a3,
		       double *b3, double *a4, double *b4,
		       double *a5, double *b5,double *c,
		       double *scatterLengths, int *type,
		       int *useN, double *n, double *delta,
		       double *sigmaC, double *sigmaS) {

  int i, j, k;

  double ** nanopC = (double **)R_alloc( *nrowC, sizeof( double * ));
  for( int i = 0; i < *nrowC; i++ )
    nanopC[i] = (double *)R_alloc(3, sizeof( double ));

  double ** nanopS = (double **)R_alloc( *nrowS, sizeof( double * ));
  for( int i = 0; i < *nrowS; i++ )
    nanopS[i] = (double *)R_alloc(3, sizeof( double ));

  double ** distC =  (double **)R_alloc( *nrowC, sizeof( double * ));
  double ** distS =  (double **)R_alloc( *nrowS, sizeof( double * ));
  double ** distCS =  (double **)R_alloc( *nrowC, sizeof( double * ));

  for( int i = 0; i < *nrowC; i++ )
    distC[i] = (double *)R_alloc(*nrowC, sizeof( double ));
  for( int i = 0; i < *nrowS; i++ )
    distS[i] = (double *)R_alloc(*nrowS, sizeof( double ));
  for( int i = 0; i < *nrowC; i++ )
    distCS[i] = (double *)R_alloc(*nrowS, sizeof( double ));

  double ** ffMC =  (double **)R_alloc( *len, sizeof( double * ));
  for( int i = 0; i < *len; i++ )
    ffMC[i] = (double *)R_alloc(*natomTypesC, sizeof( double ));
  double ** ffMS =  (double **)R_alloc( *len, sizeof( double * ));
  for( int i = 0; i < *len; i++ )
    ffMS[i] = (double *)R_alloc(*natomTypesS, sizeof( double ));


  double ff, ff1, ff2;

  double * ffAvCS = (double *)R_alloc(*len, sizeof( double ));
  double * ffAvC = (double *)R_alloc(*len, sizeof( double ));
  double * ffAvS = (double *)R_alloc(*len, sizeof( double ));

  double * fACS = (double *)R_alloc(*len, sizeof( double ));
  double * fAC = (double *)R_alloc(*len, sizeof( double ));
  double * fAS = (double *)R_alloc(*len, sizeof( double ));
  double ntmpC,ntmpS,ntmpCS, q4;

  j = 0;
  for (i=0; i < *nrowC; i++) {
    nanopC[i][0] = npC[j];
    nanopC[i][1] = npC[j+1];
    nanopC[i][2] = npC[j+2];
    j+=3;
  }
  j = 0;
  for (i=0; i < *nrowS; i++) {
    nanopS[i][0] = npS[j];
    nanopS[i][1] = npS[j+1];
    nanopS[i][2] = npS[j+2];
    j+=3;
  }

  for (j=0; j < *nrowC; j++) {
    for (k=0; k < *nrowC; k++) {
      distC[j][k] = sqrt(pow(nanopC[j][0]-nanopC[k][0],2)+
			 pow(nanopC[j][1]-nanopC[k][1],2)+
			 pow(nanopC[j][2]-nanopC[k][2],2));
    }
  }
  for (j=0; j < *nrowS; j++) {
    for (k=0; k < *nrowS; k++) {
      distS[j][k] = sqrt(pow(nanopS[j][0]-nanopS[k][0],2)+
			 pow(nanopS[j][1]-nanopS[k][1],2)+
			 pow(nanopS[j][2]-nanopS[k][2],2));
    }
  }
  for (j=0; j < *nrowC; j++) {
    for (k=0; k < *nrowS; k++) {
      distCS[j][k] = sqrt(pow(nanopC[j][0]-nanopS[k][0],2)+
			  pow(nanopC[j][1]-nanopS[k][1],2)+
			  pow(nanopC[j][2]-nanopS[k][2],2));
    }
  }
  if(*type) { // x-ray scattering
    Rprintf("Calculating X-ray scattering\n");
    for (i=0; i < *len; i++) {
      q4 = Q[i]/(4*M_PI);
      for (j=0; j < *natomTypesC; j++) {
	ffMC[i][j] = (a1[atomTypeC[j]])*exp(-(b1[atomTypeC[j]])*q4)+
	  (a2[atomTypeC[j]])*exp(-(b2[atomTypeC[j]])*q4)+
	  (a3[atomTypeC[j]])*exp(-(b3[atomTypeC[j]])*q4)+
	  (a4[atomTypeC[j]])*exp(-(b4[atomTypeC[j]])*q4)+
	  (a5[atomTypeC[j]])*exp(-(b5[atomTypeC[j]])*q4)+(c[atomTypeC[j]]);
      }
      for (j=0; j < *natomTypesS; j++) {
	  ffMS[i][j] = (a1[atomTypeS[j]])*exp(-(b1[atomTypeS[j]])*q4)+
	    (a2[atomTypeS[j]])*exp(-(b2[atomTypeS[j]])*q4)+
	    (a3[atomTypeS[j]])*exp(-(b3[atomTypeS[j]])*q4)+
	    (a4[atomTypeS[j]])*exp(-(b4[atomTypeS[j]])*q4)+
	    (a5[atomTypeS[j]])*exp(-(b5[atomTypeS[j]])*q4)+(c[atomTypeS[j]]);
      }
    }
  }
  else { // neutron scattering
    Rprintf("Calculating neutron scattering\n");
    q4 = Q[0]/(4*M_PI);
    for (i=0; i < *len; i++) {
      for (j=0; j < *natomTypesC; j++) {
	ffMC[i][j] = scatterLengths[atomTypeC[j]];
      }
      for (j=0; j < *natomTypesS; j++) {
	ffMS[i][j] = scatterLengths[atomTypeS[j]];
      }
    }
  }

  for (i=0; i < *len; i++) {
    ffAvC[i] = ffAvS[i] = ffAvCS[i] = fAC[i] = fAS[i] = fACS[i] = 0;
    for (j=0; j < *nrowC; j++) {
      for (k=0; k < *nrowC; k++) {

	q4 = Q[i]/(4*M_PI);

	ff1 = ffMC[i][atomTypeC[j]];
	ff2 = ffMC[i][atomTypeC[k]];

	ffAvC[i] += ff2*ff1;

	if(j==0)
	  fAC[i] += ff1;
      }
    }

    ffAvC[i] = ffAvC[i]/(pow((double)*nrowC,2)) ;
    fAC[i] = fAC[i]/(double)*nrowC;

    for (j=0; j < *nrowS; j++) {
      for (k=0; k < *nrowS; k++) {

	q4 = Q[i]/(4*M_PI);

	ff1 = ffMS[i][atomTypeS[j]];
	ff2 = ffMS[i][atomTypeS[k]];

	ffAvS[i] += ff2*ff1;

	if(j==0)
	  fAS[i] += ff1;
      }
    }

    ffAvS[i] = ffAvS[i]/(pow((double)*nrowS,2));
    fAS[i] = fAS[i]/(double)*nrowS;

    ffAvCS[i] = (double)(*nrowC)/(double)(*nrowC+*nrowS)*ffAvC[i] + (double)(*nrowS)/(double)(*nrowC+*nrowS)*ffAvS[i];
    fACS[i] = (double)(*nrowC)/(double)(*nrowC+*nrowS)*fAC[i] + (double)(*nrowS)/(double)(*nrowC+*nrowS)*fAS[i];

  }
  for (i=0; i < *len; i++) {

    ntmpC = ntmpS = ntmpCS = 0;

    for (j=0; j < *nrowC; j++) {
      for (k=0; k < *nrowC; k++) {

	ff1 = ffMC[i][atomTypeC[j]];
	ff2 = ffMC[i][atomTypeC[k]];

	ff = ff2*ff1;


	ntmpC += (1+res[i])*((double)*nrowC/ffAvC[i]) + pow((ff1 - fAC[i]),2);
      }
    }
    for (j=0; j < *nrowS; j++) {
      for (k=0; k < *nrowS; k++) {

	ff1 = ffMS[i][atomTypeS[j]];
	ff2 = ffMS[i][atomTypeS[k]];

	ff = ff2*ff1;

	ntmpS += (1+res[i])*((double)*nrowS/ffAvS[i]) + pow((ff1 - fAS[i]),2);
      }
    }
    for (j=0; j < *nrowC; j++) {
      for (k=0; k < *nrowS; k++) {

	ff1 = ffMC[i][atomTypeC[j]];
	ff2 = ffMS[i][atomTypeS[k]];

	ff = ff2*ff1;

	ntmpCS += (1+res[i])*( (((double)*nrowS+(double)*nrowC)/2) /ffAvCS[i]) + pow(( (ff1+ff2)/2 - fACS[i]),2);
      }
    }
    if(ntmpC < 0 || ntmpS < 0 || ntmpCS < 0)
      Rprintf("%d  %f  %f   %f  %f  %f   %f\n",i,ntmpC, ntmpCS, ntmpS,ffAvCS[i],ffAvC[i],ffAvS[i]);

    res[i] = ntmpC*(double)(*nrowC)/(double)(*nrowC+*nrowS) + ntmpCS + ntmpS*(double)(*nrowS)/(double)(*nrowC+*nrowS);
  }
  Rprintf("S: %d C: %d \n", *nrowS, *nrowC);
}




/* ***************************************************************************
	Broad PDF
	Analytically broaden a PDF using Gaussians	
	see 
	  ''Mullen, K.M. & Levin, I., 2011. Mitigation of errors in pair distribution function analysis of nanoparticles. Journal of Applied Crystallography, 44, pp.788–797.''
	for description
*************************************************************************** */
void broadPDF(double *res, double *pdf, double *r, double *rp, int *len, int *lenp, 
                 double *sigi, double *sigj, double *delta, int *n) 
{
  int i,j,k;
  double dr, x1, x2, sig, mul;

	
  dr = r[1]-r[0];  
  for (k=0; k<*lenp; k++){
    if(*n==0)
      sig = *sigi+*sigj - (*delta/rp[k])* sqrt(*sigi)*sqrt(*sigj);
    else    
      sig = (*sigi+*sigj)*(1-(*delta/pow(rp[k],*n)));
	//      sig = (*sigi+*sigj)*max(0,(1-(*delta/pow(rp[k],*n))));
    
	mul = 1.0 / (sqrt(2*M_PI* sig));
    for (i=1; i<*len; i++ ){
	  x1 = -0.5*pow((rp[k]-r[i]),2)/sig;
	  x2 = -0.5*pow((rp[k]+r[i]),2)/sig;
	  res[i] +=  (exp(x1)-exp(x2)) * pdf[k] * rp[k] * mul; 
	}
	
  } 
  
  for (i=0; i<*len; i++){
    res[i] = res[i]*dr / r[i];
  }
  

}



/* ***************************************************************************
	Exact calculation of Scattering Function S(Q) for the array of identical atoms
		
*************************************************************************** */


/* calcTotalScatt */
void calcTotalScatt(double *res, double *Q, int *len,
		    double *np, int *nrow, double *scale, double *sigma, int *type,
		    double *a1, double *b1, 
			double *a2, double *b2, 
			double *a3,  double *b3, 
			double *a4, double *b4, 
			double *a5, double *b5, 
			double *c,  double *scatterLengths, 
		    int *useN, double *n, double *delta) 

{
  int i, j, k;
  
  double ** nanop = (double **)R_alloc( *nrow, sizeof( double * ));
  for( int i = 0; i < *nrow; i++ )
    nanop[i] = (double *)R_alloc(3, sizeof( double ));
  double ** dist =  (double **)R_alloc( *nrow, sizeof( double * ));
  for( int i = 0; i < *nrow; i++ )
    dist[i] = (double *)R_alloc(*nrow, sizeof( double ));

  double * ff1 =  (double *)R_alloc( *len, sizeof( double ));
  double * ff2 =  (double *)R_alloc( *len, sizeof( double ));
  double ffAv = 0;

  double ntmp, q4, ff;

/* very strange bug. 'useN' doesn't seemed to be passed correctly */  
   if ( (*n < 1e-6) & (*delta < 1e-6) ) 
	*useN = 0;
   else
    *useN = 1;
  
  
  
  j = 0;
  for (i=0; i < *nrow; i++) { 
    nanop[i][0] = np[j];
    nanop[i][1] = np[j+1];
    nanop[i][2] = np[j+2];
    j+=3;
  }
  
  for (j=0; j < *nrow; j++) {
    for (k=0; k < *nrow; k++) { 
      dist[j][k] = sqrt(pow(nanop[j][0]-nanop[k][0],2)+ 
			pow(nanop[j][1]-nanop[k][1],2)+ 
			pow(nanop[j][2]-nanop[k][2],2));
      
    }
  }
  
  
  if(*type) { // x-ray scattering 
    for(i=0; i<*len; i++){ 
	  q4 = pow(Q[i]/(4*M_PI),2);
      ff1[i] = a1[0]*exp(-b1[0]*q4) + a2[0]*exp(-b2[0]*q4) +
	         a3[0]*exp(-b3[0]*q4) + a4[0]*exp(-b4[0]*q4) + 
		     a5[0]*exp(-b5[0]*q4)+c[0]; 
      ff2[i] =  a1[1]*exp(-b1[1]*q4) + a2[1]*exp(-b2[1]*q4) +
	         a3[1]*exp(-b3[1]*q4) + a4[1]*exp(-b4[1]*q4) + 
		     a5[1]*exp(-b5[1]*q4)+c[1]; 
	}
  }
  else { // neutron scattering 
    for(i=0; i<*len; i++){ 
      ff1[i] = scatterLengths[0];
	  ff2[i] = scatterLengths[1];
	}
  }
 
  for (i=0; i < *len; i++) {
    ntmp = 0;
    ffAv = 0;
    for (j=0; j < *nrow; j++) {
      for (k=0; k < *nrow; k++) {
	    q4 = Q[i]/(4*M_PI);
	    ff = ff1[i]*ff1[i];
		if(dist[j][k]!=0){ 
		  if(*sigma==0) // base case, no broadening
	        ntmp = ntmp + (ff *(sin(Q[i]*dist[j][k]))/(Q[i]*dist[j][k]));
		  else if(*useN==0) { //  broadening not using delta, n
			ntmp = ntmp + (ff * (exp(-.5*(*sigma+*sigma)*pow(Q[i],2)))*(sin(Q[i]*dist[j][k]))/(Q[i]*dist[j][k])); 
		  }
		  else {// use delta, n 
		  	q4 = (*sigma+*sigma)*(1.0-*delta/pow(dist[j][k],*n));
//		  	q4 = (*sigma+*sigma)*max(0,(1.0-*delta/pow(dist[j][k],*n)));
	        ntmp = ntmp + (ff * (exp(-.5*q4*pow(Q[i],2)))*(sin(Q[i]*dist[j][k]))/(Q[i]*dist[j][k]));
		  }	  
		}
	  }
	}
    res[i] =  scale[i] * ntmp;


  }

}


/* ***************************************************************************
	Fast computaion of Scattering function for the array of identical atoms
		
	_! with arbitrary accuracy but not exact !_

for the parameters meaning as well as accuracy estimation please refer to
Cervellino et. al. 
J of Comp Chem / Vol. 27, No. 9
*************************************************************************** */
void fastCalcTotalScatt(double *res, double *Q, int *len,
            double *np, int *nrow, double *scale, double *sigma, int *type,
            double *a1, double *b1, 
			double *a2, double *b2, 
			double *a3, double *b3,
			double *a4, double *b4,
		    double *a5, double *b5,
			double *c,  double *scatterLengths,
		    int *useN, double *n, double *delta,
			double *dr, double *diam, double *del, double *eps) 
{
  int i, j, k;

/*  clock_t start, end;
  double cpu_time_used;
*/
  double ** nanop = (double **)R_alloc( *nrow, sizeof( double * ));
  for( int i = 0; i < *nrow; i++ )
    nanop[i] = (double *)R_alloc(3, sizeof( double ));

  double q4, mul;

  double* weights;
  double omega, x1, x2, dist2;
  int M_d, p, S;
  double rmin = 0;
  int pseudoN;

  if(fabs(*dr)>=1e-8)
    pseudoN = floor((*diam*1.1-rmin)/ *dr);
  else
    pseudoN = *nrow * (*nrow-1)/2;
	
  double * d_l =  (double *)R_alloc( pseudoN, sizeof( double ));
  int * mu_l =  (int *)R_alloc( pseudoN, sizeof( int ));
  double * ff1 =  (double *)R_alloc( *len, sizeof( double ));
  double * ff2 =  (double *)R_alloc( *len, sizeof( double ));

/* eps: defines what to consider as similar distances for distance distribution, the lower the more precise
   M_d: number of all possible interatomic distances
   d_l: all possible different interatomic distances
   dist: calculated interatomic distances (just all of them)
   nrow: number of atoms under current consideration
   nrow2: remaining number of atoms in the particle (sorry for that notaion :) ) Zero for uniform particle. Number
          of atoms in the core if shell is under considerations, and vice versa. Is used for the scaling factor.  
   ff1: scattering factor for the atoms under consideration
   ff2: scattering factor for the other atoms   
*/  
  j = 0;
  for (i=0; i < *nrow; i++) {
    nanop[i][0] = np[j];
    nanop[i][1] = np[j+1];
    nanop[i][2] = np[j+2];
    j+=3;
  }

  if(*type) { // x-ray scattering 
    for(i=0; i<*len; i++){ 
	  q4 = pow(Q[i]/(4*M_PI),2);
      ff1[i] = a1[0]*exp(-b1[0]*q4) + a2[0]*exp(-b2[0]*q4) +
	         a3[0]*exp(-b3[0]*q4) + a4[0]*exp(-b4[0]*q4) + 
		     a5[0]*exp(-b5[0]*q4)+c[0]; 
      ff2[i] =  a1[1]*exp(-b1[1]*q4) + a2[1]*exp(-b2[1]*q4) +
	         a3[1]*exp(-b3[1]*q4) + a4[1]*exp(-b4[1]*q4) + 
		     a5[1]*exp(-b5[1]*q4)+c[1]; 
	}
  }
  else { // neutron scattering 
    for(i=0; i<*len; i++){ 
      ff1[i] = scatterLengths[0];
	  ff2[i] = scatterLengths[1];
	}
  }
     
/*   using honest  ->  */  
  if(fabs(*dr) < 1e-8){
    double * dist =  (double *)R_alloc(pseudoN, sizeof( double ));
  
    for (j=0; j < *nrow; j++) {
      for (k=0; k < j; k++) {
	    p = j*(j-1)/2 + k;
        dist[p] = sqrt(pow(nanop[j][0]-nanop[k][0],2)+
		  	  pow(nanop[j][1]-nanop[k][1],2)+
			  pow(nanop[j][2]-nanop[k][2],2));
      }
    }
  
    q_sort(dist, 0, pseudoN-1);

    for (j=1; j < pseudoN; j++)
	  mu_l[j]=0;

    p=0;
    d_l[0]=dist[0];
    mu_l[0]=1;
    for (j=0; j < pseudoN-1; j++) {
      if (dist[j+1] < dist[j] + *eps) {
	    mu_l[p]++;
  	  }
	  else{
        p++;
	    d_l[p] = dist[j+1];
	    mu_l[p] = 1;
	  }
    }
    M_d = p+1;
  }
/*   <- using honest   */  
/*   using histograms  ->  */
  else {
    for (j=0; j < pseudoN; j++){
	  mu_l[j]=0;
	  d_l[j]=rmin + (j+0.5)* *dr;
    }
	
    for (j=0; j < *nrow; j++) {
      for (k=0; k < j; k++) {
        dist2 = sqrt(pow(nanop[j][0]-nanop[k][0],2)+
		  	  pow(nanop[j][1]-nanop[k][1],2)+
			  pow(nanop[j][2]-nanop[k][2],2));
	    p = round( (dist2-rmin)/ *dr );	
	    mu_l[p]++;
      }
    }

    p = 0;
    for (j=0; j < pseudoN; j++){
      if(mu_l[j] > 0.5){
	    mu_l[p] = mu_l[j];
	    d_l[p] = d_l[j];
	    p ++; 
	  }
    }
    M_d = p;
  }
  
 /* <- using histograms */
    
  
  omega = 1.5* *del;
  S = floor((d_l[M_d-1]+omega*8.4904)/ *del);
  weights = (double*) R_alloc(S,sizeof(double));

  if (*useN==0){
    for (k=1; k <= S; k++) {
      weights[k-1] = 0;
	  for (j=0; j < M_d; j++){
	    x1 = k* *del-d_l[j];
	    x2 = k* *del+d_l[j];
	    weights[k-1] += k*pow(*del,2)*mu_l[j]/d_l[j]* ( gauss(x1, omega) - gauss(x2, omega) );
   	  }
    }
  }
  else {
    for (k=1; k <= S; k++) {
      weights[k-1] = 0;
	  for (j=0; j < M_d; j++){
	    x1 = k* *del-d_l[j];
	    x2 = k* *del+d_l[j];
	    q4 = (*sigma+*sigma)*(1.0-*delta/pow(d_l[j],*n));
//	    q4 = (*sigma+*sigma)*max(0,(1.0-*delta/pow(d_l[j],*n)));
	    weights[k-1] += k*pow(*del,2)*mu_l[j]/d_l[j]* ( gauss2(x1, omega, q4) - gauss2(x2, omega, q4) );
   	  }
    }
  }

  calcSum(res, weights, Q, S, *len, *del);
	
  for (i=0; i<*len; i++) {
    res[i] =2*ff1[i]*ff1[i]*res[i]*Fgauss(Q[i]/(2*M_PI), omega)* scale[i];
	if(*useN==0)
	    res[i] *= exp(-.5*(*sigma+*sigma)*pow(Q[i],2));

  }
}


/* fastCalcTotalScattAv */
/* ***************************************************************************
		Fast computaion of Scattering function for the array of identical atoms
        optimized by guessing interatomic distances from lattice structure;
		accuracy: arbitrary but certain error is introduced if the atoms were moved from their positions to adjust structure
                (i.e. lattice is not ideal)
				
		=== Given the array of particles from simPart, calculates average signal ===

for the parameters meaning as well as accuracy estimation please refer to
Cervellino et. al. 
J of Comp Chem / Vol. 27, No. 9
*/

		
											
							
void fastCalcTotalScattAv(int *layer_start, int *layer_end, int *layers_num,
            double *res, double *Q, int *len,
            double *np, int *nrow, double *scale, double *sigma, int *type,
            double *a1, double *b1, 
			double *a2, double *b2, 
			double *a3, double *b3,
			double *a4, double *b4,
		    double *a5, double *b5,
			double *c,  double *scatterLengths,
		    int *useN, double *n, double *delta,
			double *dr, double *diam, double *del, double *eps
			) {

  int i, j, k, l;
//  clock_t start, end;
//  double cpu_time_used, t1,t2,t3,t4;
  
  double ** nanop = (double **)R_alloc( *nrow, sizeof( double * ));
  for( int i = 0; i < *nrow; i++ )
    nanop[i] = (double *)R_alloc(3, sizeof( double ));

  double q4, mul;

  double* weights;
  double omega, x1, x2, dist2;
  int M_d, p, len_max, S;
  double rmin = 0;
  int pseudoN;
  
  if(fabs(*dr)>1e-8)
    pseudoN = floor((1.1*diam[*layers_num-1]-rmin)/ *dr);
  else
    pseudoN = *nrow * (*nrow-1)/2;
	
//  int * mu_l =  (int *)R_alloc( d_l_lens[*layers_num-1], sizeof( int ));
  double * dist =  (double *)R_alloc( *nrow * (*nrow-1)/2, sizeof( double ));
  double * res2 =  (double *)R_alloc( *len, sizeof( double ));
  double * ff1 =  (double *)R_alloc( *len, sizeof( double ));
  double * ff2 =  (double *)R_alloc( *len, sizeof( double ));
  double * d_l =  (double *)R_alloc( pseudoN, sizeof( double ));
  int * mu_l =  (int *)R_alloc( pseudoN, sizeof( int ));
  int * d_l_lens = (int *)R_alloc(*layers_num, sizeof(int));

 /* 
   eps: defines what to consider as similar distances for distance distribution, the lower the more precise
   del: defines pseudogrid, the lower the better
   M_d: number of all possible interatomic distances
   d_l: array; all possible different interatomic distances 
   layer_start: array; the arrays of atoms for the array of particles is given in a single array np,
                layer_start and layer_end defines how to cut each particle from that array.	
   layer_end: see layer_start
   layers_num: number of particles in the array   
   dist: array; contains calculated interatomic distances (just all of them)
   nrow: number of atoms under current consideration
   nrow_tot: common number of atoms in the particle. Is used for the scaling factor.  
   ff1: scattering factor for the atoms under consideration
   ff2: scattering factor for the other atoms 
   res2: temporary array for res   
*/   
  

  j = 0;
  for (i=0; i < *nrow; i++) {
    nanop[i][0] = np[j];
    nanop[i][1] = np[j+1];
    nanop[i][2] = np[j+2];
    j+=3;
  }

  if(*type) { // x-ray scattering 
    for(i=0; i<*len; i++){ 
	  q4 = pow(Q[i]/(4*M_PI),2);
      ff1[i] = a1[0]*exp(-b1[0]*q4) + a2[0]*exp(-b2[0]*q4) +
	         a3[0]*exp(-b3[0]*q4) + a4[0]*exp(-b4[0]*q4) + 
		     a5[0]*exp(-b5[0]*q4)+c[0]; 
      ff2[i] =  a1[1]*exp(-b1[1]*q4) + a2[1]*exp(-b2[1]*q4) +
	         a3[1]*exp(-b3[1]*q4) + a4[1]*exp(-b4[1]*q4) + 
		     a5[1]*exp(-b5[1]*q4)+c[1]; 
	}
  }
  else { // neutron scattering 
    for(i=0; i<*len; i++){ 
      ff1[i] = scatterLengths[0];
	  ff2[i] = scatterLengths[1];
	}
  }


  if(fabs(*dr) < 1e-8){  
    for (j=0; j < *nrow; j++) {
      for (k=0; k < j; k++) {
	    p = j*(j-1)/2 + k;
        dist[p] = sqrt(pow(nanop[j][0]-nanop[k][0],2)+
		  	  pow(nanop[j][1]-nanop[k][1],2)+
			  pow(nanop[j][2]-nanop[k][2],2));
      }
    }
  
    q_sort(dist, 0, pseudoN-1);

    p=0;
    d_l[0]=dist[0];
    for (j=0; j < pseudoN-1; j++) {
      if (dist[j+1] >= dist[j] + *eps) {
        p++;
	    d_l[p] = dist[j+1];
  	  }
    }
    M_d = p+1;
  }             
  /*   <- using honest   */  
  /*   using histograms  ->  */
  else {        
    for (j=0; j < pseudoN; j++){
	  mu_l[j]=0;
	  d_l[j]=rmin + (j+0.5)* *dr;
    }
    for (j=0; j < *nrow; j++) {
      for (k=0; k < j; k++) {
        dist2 = sqrt(pow(nanop[j][0]-nanop[k][0],2)+
		  	  pow(nanop[j][1]-nanop[k][1],2)+
			  pow(nanop[j][2]-nanop[k][2],2));
	    p = round( (dist2-rmin)/ *dr );	
	    mu_l[p]++;
      }
    }
    p = 0;
    for (j=0; j < pseudoN; j++){
      if(mu_l[j] > 0.5){
	    mu_l[p] = mu_l[j];
	    d_l[p] = d_l[j];
	    p ++; 
	  }
    }
    M_d = p;	
  }

  d_l[M_d] = d_l[M_d-1]*2;
  S = 0;
  for (j=0; j < *layers_num-1; j++) {
	for (k=S; k < M_d; k++){
	  if( (d_l[k] <= diam[j]) && (d_l[k+1] >= diam[j]) ){
	    d_l_lens[j]=k+1;
		S=k+1;
	    break;
	  }
	}
  } 
  d_l_lens[*layers_num-1]=M_d;  

  for (j=0; j < d_l_lens[*layers_num-1]; j++) {
	mu_l[j]=0;
  }  
  for (i=0; i<*len; i++) {
    res[i] =0.0;
	res2[i] = 0.0;
  }
  
/* *****************************************************
            Main cycle over all layers 
******************************************************** */
  for (l=0; l<*layers_num; l++){
    
	p=0;  
    for (j=layer_start[l]-1; j < layer_end[l]; j++) {
      for (k=layer_start[l]-1; k < j; k++) {
        dist[p] = sqrt(pow(nanop[j][0]-nanop[k][0],2)+
		    pow(nanop[j][1]-nanop[k][1],2)+
		    pow(nanop[j][2]-nanop[k][2],2));
	    p = p+1;
      }	  
    }

    omega = 1.5* *del;
    S = floor((d_l[d_l_lens[l]-1]+omega*8.4904)/ *del);
    weights = (double*) R_alloc(S,sizeof(double));

	for (i=0; i < d_l_lens[*layers_num-1]; i++){
	  mu_l[i]=0;
    }  
    for (j=0; j < p; j++) {
      for (i=0; i < d_l_lens[l]; i++){
	    if ( fabs(dist[j] - d_l[i]) < *eps) { // dist - interatomic distances; d_l - all posible interatomic distances;  
          mu_l[i]++;                         // the number of occurrences dist into d_l[i]
	      break;
	    }
	  }
    }
	if (*useN==0){
      for (k=1; k <= S; k++) {
        weights[k-1] = 0;
	    for (j=0; j < d_l_lens[l]; j++){
	      x1 = k* *del-d_l[j];
	      x2 = k* *del+d_l[j];
	      weights[k-1] += k*pow( *del,2)*mu_l[j]/d_l[j]* ( gauss(x1, omega) - gauss(x2, omega) );
   	    }
      }
    }
    else {
      for (k=1; k <= S; k++) {
        weights[k-1] = 0;
	    for (j=0; j < d_l_lens[l]; j++){
	      x1 = k* *del-d_l[j];
	      x2 = k* *del+d_l[j];
	      q4 = (*sigma+*sigma)*(1.0-*delta/pow(d_l[j],*n));
//	      q4 = (*sigma+*sigma)*max(0,(1.0-*delta/pow(d_l[j],*n)));
	      weights[k-1] += k*pow( *del,2)*mu_l[j]/d_l[j]* ( gauss2(x1, omega, q4) - gauss2(x2, omega, q4) );
     	  }
      }
    }
	
	calcSum(res2, weights, Q, S, *len,  *del);
    for (i=0; i<*len; i++) 
      res[i] = res[i] + res2[i] * scale[l* *len+i];

  }
/* *************************************************
                  End of CYCLE over l (layers)
   ************************************************* */
  for (i=0; i<*len; i++) {
    res[i] =2*ff1[i]*ff1[i]*res[i]*Fgauss(Q[i]/(2*M_PI), omega) / *layers_num;
	if(*useN==0)
	  res[i] *= exp(-.5*(*sigma+*sigma)*pow(Q[i],2));
  }
	
}

/* ***************************************************************************

   Exact computaion of Scattering function for the 2 different arrays of atoms
   
*************************************************************************** */
void calcTotalScatt_CS(double *res, double *Q, int *len, double *scale,
		       double *np_mu, double *np_nu, int *nrow_mu, int *nrow_nu, int *type,
		       double *sigma_mu, double *sigma_nu,
			   double *a1, double *b1, 
		       double *a2, double *b2, 
			   double *a3, double *b3, 
			   double *a4, double *b4, 
			   double *a5, double *b5, 
			   double *c,  double *scatterLengths, 
			   int *useN, double *n, double *delta)
{

  int i, j, k;
  
  double ** nanop_mu = (double **)R_alloc( *nrow_mu, sizeof( double * ));
  for( int i = 0; i < *nrow_mu; i++ )
    nanop_mu[i] = (double *)R_alloc(3, sizeof( double ));
  double ** nanop_nu = (double **)R_alloc( *nrow_nu, sizeof( double * ));
  for( int i = 0; i < *nrow_nu; i++ )
    nanop_nu[i] = (double *)R_alloc(3, sizeof( double ));

  double ** distCS =  (double **)R_alloc( *nrow_mu, sizeof( double * ));
  for( int i = 0; i < *nrow_mu; i++ )
    distCS[i] = (double *)R_alloc(*nrow_nu, sizeof( double ));
  
  double * ff1 =  (double *)R_alloc( *len, sizeof( double ));
  double * ff2 =  (double *)R_alloc( *len, sizeof( double ));
  
  
  double ff;

  double ntmpCS, q4;

    j = 0;
    for (i=0; i < *nrow_mu; i++) {
  	  nanop_mu[i][0] = np_mu[j];
	  nanop_mu[i][1] = np_mu[j+1];
	  nanop_mu[i][2] = np_mu[j+2];
	  j+=3;
    }
    j = 0;
    for (i=0; i < *nrow_nu; i++) {
	  nanop_nu[i][0] = np_nu[j];
	  nanop_nu[i][1] = np_nu[j+1];
	  nanop_nu[i][2] = np_nu[j+2];
	  j+=3;
    }
  


  for (j=0; j < *nrow_mu; j++) {
    for (k=0; k < *nrow_nu; k++) { 
      distCS[j][k] = sqrt(pow(nanop_mu[j][0]-nanop_nu[k][0],2)+ 
			  pow(nanop_mu[j][1]-nanop_nu[k][1],2)+ 
			  pow(nanop_mu[j][2]-nanop_nu[k][2],2));
    }
  }
  
  
  if(*type) { // x-ray scattering 
    for(i=0; i<*len; i++){ 
	  q4 = pow(Q[i]/(4*M_PI),2);
      ff1[i] = a1[0]*exp(-b1[0]*q4) + a2[0]*exp(-b2[0]*q4) +
	         a3[0]*exp(-b3[0]*q4) + a4[0]*exp(-b4[0]*q4) + 
		     a5[0]*exp(-b5[0]*q4)+c[0]; 
      ff2[i] = a1[1]*exp(-b1[1]*q4) + a2[1]*exp(-b2[1]*q4) +
	       a3[1]*exp(-b3[1]*q4) + a4[1]*exp(-b4[1]*q4) + 
		   a5[1]*exp(-b5[1]*q4)+c[1];   
	}
  }
  else { // neutron scattering 
    for(i=0; i<*len; i++){ 
      ff1[i] = scatterLengths[0]; 
	  ff2[i] = scatterLengths[1];
	}
  }
  
  
  
  
  for (i=0; i < *len; i++) {
    ntmpCS = 0;    
    for (j=0; j < *nrow_mu; j++) {
      for (k=0; k < *nrow_nu; k++) { 
		ff = ff2[i]*ff1[i];
		if(distCS[j][k]!=0) {
          if(*sigma_mu==0&&*sigma_nu==0) // base case 
	        ntmpCS = ntmpCS + (ff *(sin(Q[i]*distCS[j][k]))/(Q[i]*distCS[j][k]));
	      else if(*useN==0) { //  not using delta, n
	        ntmpCS = ntmpCS + (ff * (exp(-.5*(*sigma_mu+*sigma_nu)*pow(Q[i],2)))*(sin(Q[i]*distCS[j][k]))/(Q[i]*distCS[j][k])); 
	      }
	      else {// use delta, n 
	        q4 = (*sigma_mu+*sigma_nu)*(1.0-*delta/pow(distCS[j][k],*n));
//	        q4 = (*sigma_mu+*sigma_nu)*max(0,(1.0-*delta/pow(distCS[j][k],*n)));
	        ntmpCS = ntmpCS + (ff * (exp(-.5*q4*pow(Q[i],2)))*(sin(Q[i]*distCS[j][k]))/(Q[i]*distCS[j][k]));
	      }		
		}
      }
    }
    res[i] =  2*ntmpCS* scale[i]; //(*nrow_nu + *nrow_mu) / pow(*nrow_mu*ff1[i] + *nrow_nu*ff2[i], 2);
	
  }
}

/* ***************************************************************************
		Fast computaion of Scattering function for the 2 different arrays of atoms
	    _! with arbitrary accuracy but not exact !_

for the parameters meaning as well as accuracy estimation please refer to
Cervellino et. al. 
J of Comp Chem / Vol. 27, No. 9
*/

void fastCalcTotalScatt_CS(double *res, double *Q, int *len, double *scale,
            double *np_mu, double *np_nu, int *nrow_mu, int *nrow_nu, int *type,
		    double *sigma_mu, double *sigma_nu,
			double *a1, double *b1, 
			double *a2, double *b2, 
			double *a3, double *b3, 
			double *a4, double *b4, 
		    double *a5, double *b5, 
			double *c,  double *scatterLengths, 
			int *useN, double *n, double *delta,
    	    double *dr, double *rmax, double *del, double *eps) 
{

  int i, j, k;

  double sigma, q4;
  double x1, x2;

  double* weights;
  double  omega;
  int M_d, p, pp, S;
  double rmin = 0;
  int pseudoN;
  		
  if(fabs(*dr)>1e-8)
    pseudoN = floor((*rmax-rmin)/ *dr);
  else
    pseudoN = *nrow_nu* *nrow_mu;

	
/*  clock_t start, end;
  double cpu_time_used;
*/

  double rel, mul;

  double * d_l =  (double *)R_alloc( pseudoN, sizeof( double ));
  int * mu_l =  (int *)R_alloc( pseudoN, sizeof( int ));
  double  dist2;
  double * ff1 =  (double *)R_alloc( *len, sizeof( double ));
  double * ff2 =  (double *)R_alloc( *len, sizeof( double ));

  double ** nanop_mu = (double **)R_alloc( *nrow_mu, sizeof( double * ));
  for( int i = 0; i < *nrow_mu; i++ )
    nanop_mu[i] = (double *)R_alloc(3, sizeof( double ));
  double ** nanop_nu = (double **)R_alloc( *nrow_nu, sizeof( double * ));
  for( int i = 0; i < *nrow_nu; i++ )
    nanop_nu[i] = (double *)R_alloc(3, sizeof( double ));

	
  if(*type) { // x-ray scattering 
    for(i=0; i<*len; i++){ 
	  q4 = pow(Q[i]/(4*M_PI),2);
      ff1[i] = a1[0]*exp(-b1[0]*q4) + a2[0]*exp(-b2[0]*q4) +
	         a3[0]*exp(-b3[0]*q4) + a4[0]*exp(-b4[0]*q4) + 
		     a5[0]*exp(-b5[0]*q4)+c[0]; 
      ff2[i] = a1[1]*exp(-b1[1]*q4) + a2[1]*exp(-b2[1]*q4) +
	       a3[1]*exp(-b3[1]*q4) + a4[1]*exp(-b4[1]*q4) + 
		   a5[1]*exp(-b5[1]*q4)+c[1];   
	}
  }
  else { // neutron scattering 
    for(i=0; i<*len; i++){ 
      ff1[i] = scatterLengths[0]; 
	  ff2[i] = scatterLengths[1];
	}
  }
  

  j = 0;
  for (i=0; i < *nrow_mu; i++) {
    nanop_mu[i][0] = np_mu[j];
    nanop_mu[i][1] = np_mu[j+1];
    nanop_mu[i][2] = np_mu[j+2];
	j+=3;
  }
  j = 0;
  for (i=0; i < *nrow_nu; i++) {
	nanop_nu[i][0] = np_nu[j];
	nanop_nu[i][1] = np_nu[j+1];
	nanop_nu[i][2] = np_nu[j+2];
	j+=3;
  }

  
/*   using honest  ->  */  
  if(fabs(*dr) < 1e-8){	
    double * dist =  (double *)R_alloc( pseudoN, sizeof( double ));

    for (j=0; j < *nrow_nu; j++) {
	  for (k=0; k < *nrow_mu; k++) {
	    pp = j* *nrow_mu+k;
		p = pp/(*nrow_nu) + *nrow_mu*(pp % *nrow_nu);
	    dist[p] = sqrt(pow(nanop_nu[j][0]-nanop_mu[k][0],2)+
		  	  pow(nanop_nu[j][1]-nanop_mu[k][1],2)+
			  pow(nanop_nu[j][2]-nanop_mu[k][2],2));
	  }
    }

    q_sort(dist, 0, pseudoN-1);
	
    for (j=0; j < pseudoN; j++) {
	  mu_l[j]=0;
    }

    p=0;
    d_l[0]=dist[0];
    mu_l[0]=1;
    for (j=0; j < pseudoN-1; j++) {
      if (dist[j+1] < dist[j] + *eps) {
	    mu_l[p]++;
	  }
	  else{
        p++;
	    d_l[p] = dist[j+1];
	    mu_l[p] = 1;
  	  }
    }
    M_d = p+1;
  }
/*   <- using honest   */  
/*   using histograms  ->  */ 
  else {
    for (j=0; j < pseudoN; j++){
	  mu_l[j]=0;
	  d_l[j]=rmin + (j+0.5)* *dr;
    }
    for (j=0; j < *nrow_nu; j++) {
	  for (k=0; k < *nrow_mu; k++) {
	    dist2 = sqrt(pow(nanop_nu[j][0]-nanop_mu[k][0],2)+
		  	  pow(nanop_nu[j][1]-nanop_mu[k][1],2)+
			  pow(nanop_nu[j][2]-nanop_mu[k][2],2));
	    p = round( (dist2-rmin)/ *dr );	
	    mu_l[p]++;
      }
    }

    p = 0;
    for (j=0; j < pseudoN; j++){
      if(mu_l[j] > 0.5){
	    mu_l[p] = mu_l[j];
	    d_l[p] = d_l[j];
	    p ++; 
	  }
    }
    M_d = p;
  }
/* <- using histograms */  
  
  
  // NOW: d_l --  various distances
  //      mu_l -- multiplicity of the l-th distance d_l
  //      M_d -- number of distinct nonzero interatomic distances


  omega = 1.5* *del;
  S = floor((d_l[M_d-1]+omega*8.4904)/ *del);
  double * divi =  (double *)R_alloc( M_d, sizeof( double ));
  for (j=0; j < M_d; j++){
	divi[j] = mu_l[j]/d_l[j];
  }
  
  rel = 1/ (2.0*omega*omega);
  weights = (double*) R_alloc(S,sizeof(double));
  
  if (*useN==0){
    for (k=1; k <= S; k++) {
      weights[k-1] = 0;
	  mul = k*pow(*del,2)/(omega*sqrt(2.0*M_PI));
	  for (j=0; j < M_d; j++){
	    x1 = k* *del-d_l[j];
	    x2 = k* *del+d_l[j];
	    weights[k-1] += mul* divi[j] * ( exp(-pow(x1,2)*rel) - exp(-pow(x2,2)*rel) );
	  }
    }
  }
  else {
    for (k=1; k <= S; k++) {
      weights[k-1] = 0;
	  mul = k*pow(*del,2);
	  for (j=0; j < M_d; j++){
	    x1 = k* *del-d_l[j];
	    x2 = k* *del+d_l[j];
	    q4 = (*sigma_mu + *sigma_nu)*(1.0-*delta/pow(d_l[j],*n));
//	    q4 = (*sigma_mu + *sigma_nu)*max(0,(1.0-*delta/pow(d_l[j],*n)));
	    weights[k-1] += mul* divi[j] *( gauss2(x1, omega, q4) - gauss2(x2, omega, q4) );
   	  }
    }
  }
  
  sigma = (*sigma_mu + *sigma_nu);
  calcSum(res, weights, Q, S, *len, *del);

  for (i=0; i<*len; i++) {
 //   mul = (*nrow_nu + *nrow_mu) / pow(*nrow_mu*ff1[i] + *nrow_nu*ff2[i], 2);
    res[i] = 2*ff1[i]*ff2[i]*res[i]*Fgauss(Q[i]/(2*M_PI), omega)* scale[i];
	if(*useN==0)
	   res[i] *= exp(-.5*sigma*pow(Q[i],2));

  }

}



/* ***************************************************************************
		very fast computaion of Scattering function for the 2 different arrays of atoms
		calculates S(Q) with arbitrary accuracy but not exact
		!no assumtion about interatomics distances is made! (in contrast to the fastCalcTotalScattAv) so no additional error is introduced
		So fast is only because we calculate all possible interatomic distances for all particles in the array once and then use it for all layers
		
notes about parameters and accuracy

due to 
Cervellino et. al. 
J of Comp Chem / Vol. 27, No. 9

*/
			
void fastCalcTotalScattAv_CS(int *layer_start, int *layer_end, int *layerS_start, int *layerS_end, int *layers_num, 
            double *res, double *Q, int *len, double *scale,
            double *np_mu, double *np_nu, int *nrow_mu, int *nrow_nu, int *type,
		    double *sigma_mu, double *sigma_nu,
			double *a1, double *b1, 
			double *a2, double *b2, 
			double *a3, double *b3, 
			double *a4, double *b4, 
		    double *a5, double *b5, 
			double *c,  double *scatterLengths, 
			int *useN, double *n, double *delta, double *del){

  int i, j, k, l;

  double sigma, q4;
  double x1, x2;

  double* weights;
  double  omega;
  int M_d, p, p1, p2, len_max, S;
//  const double eps=1e-3;    // defines what to count as similar distances for distance distribution
//  if(Q[*len-1] < 1.0)
//    *del *=7;
  clock_t start, end;
  double cpu_time_used, t1,t2,t3,t4, t5;


  double rel, mul;

  double * dist =  (double *)R_alloc( *nrow_nu* *nrow_mu, sizeof( double ));
  double * d_l =  (double *)R_alloc( *nrow_nu* *nrow_mu, sizeof( double ));
  double * res2 =  (double *)R_alloc( *len, sizeof( double ));
  double * ff1 =  (double *)R_alloc( *len, sizeof( double ));
  double * ff2 =  (double *)R_alloc( *len, sizeof( double ));


  double ** nanop_mu = (double **)R_alloc( *nrow_mu, sizeof( double * ));
  for( int i = 0; i < *nrow_mu; i++ )
    nanop_mu[i] = (double *)R_alloc(3, sizeof( double ));
  double ** nanop_nu = (double **)R_alloc( *nrow_nu, sizeof( double * ));
  for( int i = 0; i < *nrow_nu; i++ )
    nanop_nu[i] = (double *)R_alloc(3, sizeof( double ));
	
	
  if(*type) { // x-ray scattering 
    for(i=0; i<*len; i++){ 
	  q4 = pow(Q[i]/(4*M_PI),2);
      ff1[i] = a1[0]*exp(-b1[0]*q4) + a2[0]*exp(-b2[0]*q4) +
	         a3[0]*exp(-b3[0]*q4) + a4[0]*exp(-b4[0]*q4) + 
		     a5[0]*exp(-b5[0]*q4)+c[0]; 
      ff2[i] = a1[1]*exp(-b1[1]*q4) + a2[1]*exp(-b2[1]*q4) +
	       a3[1]*exp(-b3[1]*q4) + a4[1]*exp(-b4[1]*q4) + 
		   a5[1]*exp(-b5[1]*q4)+c[1];   
	}
  }
  else { // neutron scattering 
    for(i=0; i<*len; i++){ 
      ff1[i] = scatterLengths[0]; 
	  ff2[i] = scatterLengths[1];
	}
  }
  

  j = 0;
  for (i=0; i < *nrow_mu; i++) {
    nanop_mu[i][0] = np_mu[j];
	nanop_mu[i][1] = np_mu[j+1];
	nanop_mu[i][2] = np_mu[j+2];
	j+=3;
  }
  j = 0;
  for (i=0; i < *nrow_nu; i++) {
    nanop_nu[i][0] = np_nu[j];
	nanop_nu[i][1] = np_nu[j+1];
	nanop_nu[i][2] = np_nu[j+2];
	j+=3;
  }

  for (i=0; i<*len; i++) {
    res[i] =0.0;
  }
  
  
/* Calculation of interatomic distances
   dist - array, will be sorted later
*/
  len_max = *nrow_nu* *nrow_mu;	
  p = 0;	
  

  for (j=layerS_start[0]-1; j < layerS_end[0]; j++) {
    for (k=layer_start[0]-1; k < layer_end[0]; k++) {
      dist[p] = sqrt(pow(nanop_nu[j][0]-nanop_mu[k][0],2)+
	 	  pow(nanop_nu[j][1]-nanop_mu[k][1],2)+
	   	  pow(nanop_nu[j][2]-nanop_mu[k][2],2));
	  p = p+1;
    } 
  }
  

  for (i=1; i < *layers_num; i++){
    for (j=layerS_start[i]-1; j < layerS_end[i]; j++) {
      for (k=layer_end[i-1]; k < layer_end[i] ; k++) {
        dist[p] = sqrt(pow(nanop_nu[j][0]-nanop_mu[k][0],2)+
	   	  pow(nanop_nu[j][1]-nanop_mu[k][1],2)+  /* new_layer_to_core x new_shell */
	   	  pow(nanop_nu[j][2]-nanop_mu[k][2],2));
	    p = p+1;
      } 
    }
    for (j=layerS_end[i-1]; j < layerS_end[i]; j++) {
      for (k=layer_start[0]-1; k < layer_end[i-1] ; k++) {
        dist[p] = sqrt(pow(nanop_nu[j][0]-nanop_mu[k][0],2)+
	   	  pow(nanop_nu[j][1]-nanop_mu[k][1],2)+  /* old_core x new_layer_to_shell */
	   	  pow(nanop_nu[j][2]-nanop_mu[k][2],2));
	    p = p+1;
      } 
    }
  }
  M_d = p;

  q_sort(dist, 0, M_d-1);

/* d_l - all possible interatomic distances, sorted
         (that doesn't matter, actually) */  
  p=0;
  d_l[0]=dist[0];
  for (j=0; j < M_d-1; j++) {
    if (dist[j+1] > dist[j] + 1e-7) {
      p++;
	  d_l[p] = dist[j+1];
	}
  }
  M_d = p+1;


  omega = 1.5* *del;
  S = floor((d_l[M_d-1]+omega*8.4904)/ *del);
  
  weights = (double*) R_alloc(S,sizeof(double));
  int * mu_l =  (int *)R_alloc( M_d, sizeof( int ));

  
  double ** F = (double **)R_alloc( S, sizeof( double * ));
  for( int i = 0; i < S; i++ )
    F[i] = (double *)R_alloc(M_d, sizeof( double ));	
	
  rel = 1/ (2.0*omega*omega);
  
  if (*useN==0){
    for (k=1; k <= S; k++) {
	  mul = k*pow(*del,2)/(omega*sqrt(2.0*M_PI));
	  for (j=0; j < M_d; j++){
	    x1 = k* *del-d_l[j];
	    x2 = k* *del+d_l[j];
	    F[k-1][j] = mul/d_l[j]*( exp(-pow(x1,2)*rel) - exp(-pow(x2,2)*rel) );
	  }
    }
  }
  else {
    for (k=1; k <= S; k++) {
      weights[k-1] = 0;
	  mul = k*pow(*del,2);
	  for (j=0; j < M_d; j++){
	    x1 = k* *del-d_l[j];
	    x2 = k* *del+d_l[j];
	    q4 = (*sigma_mu + *sigma_nu)*(1.0-*delta/pow(d_l[j],*n));
//	    q4 = (*sigma_mu + *sigma_nu)*max(0,(1.0-*delta/pow(d_l[j],*n)));
		F[k-1][j] = mul/d_l[j]*( gauss2(x1, omega, q4) - gauss2(x2, omega, q4) );
   	  }
    }
  }
  
  
/* *****************************************************
            Main cycle over all layers 
******************************************************** */  
  for (l=0; l<*layers_num; l++){
  
    p = 0;
    for (j=layerS_start[l]-1; j < layerS_end[l]; j++) {
      for (k=layer_start[l]-1; k < layer_end[l]; k++) {
        dist[p] = sqrt(pow(nanop_nu[j][0]-nanop_mu[k][0],2)+
	   	    pow(nanop_nu[j][1]-nanop_mu[k][1],2)+
	   	    pow(nanop_nu[j][2]-nanop_mu[k][2],2));
	    p = p+1;
      } 
    }

	q_sort(dist, 0, p-1);
  
    for (j=0; j < M_d; j++){
	  mu_l[j]=0;
    }
    j = 0;
	i = 0;
	
    while((j < p) && (i < M_d)){	    
	  if ( fabs(dist[j] - d_l[i]) < 1e-7){  
	    mu_l[i]++;
		j++;
	  }
	  else
		i++;
	}
	
	double * divi =  (double *)R_alloc( M_d, sizeof( double ));
    for (j=0; j < M_d; j++){
	  divi[j] = mu_l[j]/d_l[j];
    }
	
    for (k=1; k <= S; k++) {
      weights[k-1] = 0;
	  for (j=0; j < M_d; j++){
	    weights[k-1] += mu_l[j]*F[k-1][j];
	  }
    }

	
    sigma = (*sigma_mu + *sigma_nu);

	calcSum(res2, weights, Q, S, *len, *del);
	   
//    p1 = layer_end[l] - layer_start[l] + 1;
//    p2 = layerS_end[l] - layerS_start[l] + 1;

    for (i=0; i<*len; i++) {
      res[i] = res[i] + res2[i]* scale[l* *len + i]; // * (p1+p2) / pow(p1*ff1[i] + p2*ff2[i], 2);
    }  
	
  }
  
/* *****************************************************
            End of cycle
******************************************************** */ 


  for (i=0; i<*len; i++) {
    res[i] =2*ff1[i]*ff2[i]*res[i]*Fgauss(Q[i]/(2*M_PI), omega)/ (*layers_num);
    if(*useN==0)
	  res[i] *= exp(-.5*sigma*pow(Q[i],2));
  }
  
  
}

/* *********************************************
   Simple Gauss function */
double gauss(double x, double omega){
  double g;
  g = 1.0/(omega*sqrt(2.0*M_PI))*exp(-x*x/(2.0*omega*omega));
  return g;
}


/* *********************************************
   Simple Gauss function 2 */
double gauss2(double x, double omega, double sigma){
  double g;
  g = 1.0/(sqrt(sigma+omega*omega)* sqrt(2.0*M_PI))*exp(-x*x/(2.0*(sigma+omega*omega)));
  return g;
}


/* *********************************************
   Convoluted Gauss function */
double Fgauss(double q, double omega){
  double c;
  c = exp(2.0*M_PI*M_PI*omega*omega*q*q);
  return c;
}

/* *********************************************
   max function */
double max(double a, double b){
  double c;
  if(a>b)
    c = a;
  return c;
}


/* *********************************************
    Calculates weights for 
   *********************************************  */
void calcSum(double res[], double weights[], double Q[], int S, int len, double del){
int n, k;
double sn,cn, ukm2, ukm1,uk;

 for (n=0; n<len; n++) {
    sn = sin(Q[n]*del)/(Q[n]*del);
	cn = cos(Q[n]*del);
	ukm2 = 1;
	ukm1 = 2*cn;

	res[n]= weights[0] + weights[1]*cn;
	for (k=2; k<S; k++){
	  uk =  2*cn * ukm1 - ukm2;
	  res[n] += weights[k]/(k+1) * uk;
	  ukm2 = ukm1;
	  ukm1 = uk;
	}

  }
}

/* *********************************************
                 QUICK SORT 
   *********************************************  */
void q_sort(double numbers[], int left, int right){

	int l_hold, r_hold;
	double pivot;

	if (right-left==1) {  // if only 2 elements
		if (numbers[left]>numbers[right]) {
			pivot = numbers[left];
			numbers[left]=numbers[right];
			numbers[right] = pivot;
			return;
		}
	}
	if (right-left==2) { // if only 3 elements
		if(numbers[left]>numbers[left+1]) {
			// 3 possible  BAC, CAB, CBA
			if (numbers[left]>numbers[right]) { // CAB, CBA
				if (numbers[left+1]>numbers[right]) { // CBA
					pivot = numbers[left];
					numbers[left] = numbers[right];
					numbers[right] = pivot;
				}
				else { // CAB
					pivot = numbers[left];
					numbers[left] = numbers[left+1];
					numbers[left+1] = numbers[right];
					numbers[right] = pivot;
				}
			}
			else { // BAC
				pivot = numbers[left];
				numbers[left]=numbers[left+1];
				numbers[left+1] = pivot;
			}
		}
		else { // ABC, ACB, BCA
			if (numbers[left]>numbers[right]) { // BCA
				pivot = numbers[left];
				numbers[left] = numbers[right];
				numbers[right] = numbers[left+1];
				numbers[left+1] = pivot;
			}
			else { // ABC, ACB
				if (numbers[left+1]>numbers[right]) { // ACB
					pivot = numbers[left+1];
					numbers[left+1]=numbers[right];
					numbers[right] = pivot;
				}
			}
		}
	}

	l_hold = left;
	r_hold = right;
	pivot = numbers[left];

	while (left < right)
	{
		while ((numbers[right] >= pivot) && (left < right)) {
			right--;
		}
		if (left != right)
		{
			numbers[left] = numbers[right];
			left++;
		}
		while ((numbers[left] <= pivot) && (left < right)) {
			left++;
		}
		if (left != right)
		{
			numbers[right] = numbers[left];
			right--;
		}
	}
	numbers[left] = pivot;

	if (l_hold < left)
		q_sort(numbers, l_hold, left-1);
	if (r_hold > left)
		q_sort(numbers, left+1, r_hold);
}
