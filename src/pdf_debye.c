#include <R.h>
#include <Rmath.h>
#include <Rdefines.h> 


/* calcPDF */

void calcPDF(double *res, double *r, int *len, double *np, int *nrow,
	     double *calpha, double *dr, double *minR, double *p) {

  int i, j, k;
  double ** nanop = Calloc( *nrow * sizeof( double * ) , double*);
  for( int i = 0; i < *nrow; i++ )
    nanop[ i ] = Calloc(3 * sizeof( double ) , double);
 
  double dist, ntmp;

  j = 0;
  for (i=0; i < *nrow; i++) { 
    nanop[i][0] = np[j];
    nanop[i][1] = np[j+1];
    nanop[i][2] = np[j+2];
    j+=3;
  }
  for (j=0; j < *nrow; j++) {
    for (k=0; k < *nrow; k++) { 
	
      dist = sqrt(pow(nanop[j][0]-nanop[k][0],2)+ 
		  pow(nanop[j][1]-nanop[k][1],2)+ 
		  pow(nanop[j][2]-nanop[k][2],2));
      
      i = round( (dist-*minR)/ *dr);
      if(i >= 0 && i < *len) 
	res[i]++;
    } 
  }
  for (i=0; i < *len; i++) {
    res[i] = (res[i] /  (double)*nrow) / 
      (4.0*M_PI*pow(r[i],2)* *dr * *p * *calpha);
  }
}

/* calcPDF_CS */

void calcPDF_CS(double *res, double *r, int *len, double *npC, int
		*nrowC, double *npS, int *nrowS, double *calpha,
		double *dr, double *minR, double *p) {

  int i, j, k;
  double nanopC[*nrowC][3];
  double nanopS[*nrowS][3];
  double dist, ntmp;

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
  
  for (j=0; j < *nrowS; j++) {
    for (k=0; k < *nrowC; k++) { 
      dist = sqrt(pow(nanopS[j][0]-nanopC[k][0],2)+ 
			pow(nanopS[j][1]-nanopC[k][1],2)+ 
			pow(nanopS[j][2]-nanopC[k][2],2));
      i = round( (dist-*minR)/ *dr);
      if(i >= 0 && i < *len) 
	res[i]++;
    } 
  }
  // mult by two since do not loop by matrix twice
  for (i=0; i < *len; i++) {
    res[i] = 2 * (res[i] / ( (double)*nrowC+(double)*nrowS)) / 
      (4.0*M_PI*pow(r[i],2)* *dr * *p * *calpha);
  }
}

/* calcTotalScatt */

void calcTotalScatt(double *res, double *Q, int *len, double *minQ,
		    double *dQ, double *np, int *nrow, double *a1,
		    double *b1, double *a2, double *b2, double *a3,
		    double *b3, double *a4, double *b4, double *c) {

  int i, j, k;
  
  double ** nanop = Calloc( *nrow * sizeof( double * ) , double*);
  for( int i = 0; i < *nrow; i++ )
    nanop[ i ] = Calloc(3 * sizeof( double ) , double);
  double ** dist = Calloc( *nrow * sizeof( double * ) , double*);
  for( int i = 0; i < *nrow; i++ )
    dist[ i ] = Calloc(*nrow * sizeof( double ) , double);

  double *ff = Calloc( *len * sizeof( double) , double);
  
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
      (*a4)*exp(-(*b4)*q4)+(*c);
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
	if(Q[i]*dist[j][k]!=0)
	  res[i] = res[i] + ff[i] * 
	    (sin(Q[i]*dist[j][k]))/(Q[i]*dist[j][k]);
	else 
	  res[i] = res[i] + ff[i];
      }
    }
  }
}


/* calcQDepPDF */

void calcQDepPDF(double *res, 
		    double *Q, double *r, 
		    int *len, double *np, int *nrow,
		    double *a1, double *b1, 
		    double *a2, double *b2,
		    double *a3, double *b3, 
		    double *a4, double *b4,
		    double *c) {
  

  int i, j, k;
  
  double ** nanop = Calloc( *nrow * sizeof( double * ) , double*);
  for( int i = 0; i < *nrow; i++ )
    nanop[ i ] = Calloc(3 * sizeof( double ) , double);
  double ** dist = Calloc( *nrow * sizeof( double * ) , double*);
  for( int i = 0; i < *nrow; i++ )
    dist[ i ] = Calloc(*nrow * sizeof( double ) , double);

  double *ff = Calloc( *len * sizeof( double) , double);
  
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
      (*a4)*exp(-(*b4)*q4)+(*c);
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
    ntmp = 0;
    for (j=0; j < *nrow; j++) {
      for (k=0; k < *nrow; k++) { 
	if(dist[j][k]!=0) 
	  ntmp = ntmp + (ff[i] *(sin(Q[i]*dist[j][k]))/(Q[i]*dist[j][k]));
      }
    }
    res[i] = 1/((double)*nrow*ff[i]) * ntmp * sin(Q[i]* *r) * Q[i];
  }
 
}
/* calcQDepPDFAux */

void calcQDepPDFAux(double *res, 
		    double *Q, 
		    int *len, double *np, int *nrow,
		    double *a1, double *b1, 
		    double *a2, double *b2,
		    double *a3, double *b3, 
		    double *a4, double *b4,
		    double *c) {
  
  int i, j, k;
  
  double ** nanop = Calloc( *nrow * sizeof( double * ) , double*);
  for( int i = 0; i < *nrow; i++ )
    nanop[ i ] = Calloc(3 * sizeof( double ) , double);
  double ** dist = Calloc( *nrow * sizeof( double * ) , double*);
  for( int i = 0; i < *nrow; i++ )
    dist[ i ] = Calloc(*nrow * sizeof( double ) , double);

  double *ff = Calloc( *len * sizeof( double) , double);
  
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
      (*a4)*exp(-(*b4)*q4)+(*c);
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
    ntmp = 0;
    for (j=0; j < *nrow; j++) {
      for (k=0; k < *nrow; k++) { 
	if(dist[j][k]!=0) 
	  ntmp = ntmp + (ff[i] *(sin(Q[i]*dist[j][k]))/(Q[i]*dist[j][k]));
      }
    }
    res[i] = 1/((double)*nrow*ff[i]) * ntmp;
  }
}

/* calcRedTotalScatt */
void calcRedTotalScatt(double *res, double *Q, int *len, double *minQ,
		       double *dQ, double *np, int *nrow, double *a1,
		       double *b1, double *a2, double *b2, double *a3,
		       double *b3, double *a4, double *b4, double *c)
		       {
  
  int i, j, k;

  
  double ** nanop = Calloc( *nrow * sizeof( double * ) , double*);
  for( int i = 0; i < *nrow; i++ )
    nanop[ i ] = Calloc(3 * sizeof( double ) , double);
  double ** dist = Calloc( *nrow * sizeof( double * ) , double*);
  for( int i = 0; i < *nrow; i++ )
    dist[ i ] = Calloc(*nrow * sizeof( double ) , double);

  double *ff = Calloc( *len * sizeof( double) , double);
  
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
      (*a4)*exp(-(*b4)*q4)+(*c);
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

void simPart(double *res, int *lenres, double *base, int *lenbase, 
	     double *shiftx, double *shifty, double *shiftz, 
	     int *a, int *b, int *c) {

  double sh;
  int i, j, k, h, cnt;
  cnt = 0;
 
  for (i=-*a; i < *a+1; i++) { 
     for (j=-*b; j < *b+1; j++) { 
       for (k=-*c; k < *c+1; k++) { 

	 for (h=0; h < (*lenbase-2); h+=3) { 
	   res[cnt] = base[h]+i**shiftx;
	   res[cnt+1] = base[h+1]+j**shifty;
	   res[cnt+2] = base[h+2]+k**shiftz;
	   
	   cnt+=3;
	 }
       }
     }
  }
}
/* calcTotalScatt */

void calcTotalScatterSimple(double *res, double *Q, int *len, double *minQ,
		    double *dQ, double *np, int *nrow) {

  int i, j, k;
   
  double ** nanop = Calloc( *nrow * sizeof( double * ) , double*);
  for( int i = 0; i < *nrow; i++ )
    nanop[ i ] = Calloc(3 * sizeof( double ) , double);
  double ** dist = Calloc( *nrow * sizeof( double * ) , double*);
  for( int i = 0; i < *nrow; i++ )
    dist[ i ] = Calloc(*nrow * sizeof( double ) , double);
 
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
	if(Q[i]*dist[j][k]!=0)
	  res[i] = res[i] * 
	    (sin(Q[i]*dist[j][k]))/(Q[i]*dist[j][k]);

      }
    }
  }
}
