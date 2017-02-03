#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "panelIA.h"
#include "Laguerre.h"
#include "Jacobi.h"
#include "Hermite.h"

#define SMOOTH 1

//implementation: make quadjac
//		./quadjac -p=%whatever desired p value is

void paramcircle(double theta, double *x);
void paramsquare(double theta, double *x, int i);
void paramtriangle(double theta, double *x, int i);
void paramellipse(double theta, double *x);
double f2(double *x, double t);
void getfquad(double h, int i, int p, panel *pX, double xLeg[], double wLeg[], double solnquad[]);
double l2err(double h, int i, int p, panel *pX, double rhs, double xLeg[], double wLeg[]);
double getfnode(double h, int i, int p, double *x);
/* LU factorization */
extern void  dgetrf_(int *nRow, int *nCol, double *A, int *lda, int *pvt, int *inf); 
/* LU solve */
extern void  dgetrs_(char *trans, int *n, int *nRhs, double *A, int *lda, int *pvt, double *B, int *ldb, int *inf); 
void MtVadd( int n, int m, double *A, double *x, double *y );
void MtVsub( int n, int m, double *A, double *x, double *y );
double g(double x[], double t, double ny[]);
void getdubdn(double h, int i, panel *p, double dubdnsoln[]);
double max(double myArray[], int size);

double x0[2];

void main(int nargs, char *argv[]) {

  double h, alph, theta, tmax=1.0;
  double *xLeg, *wLeg;
  double **vtx;
  panel *pnls;
  double xdiff, ydiff, norm, t0, t1, ny[2], fnodal;
  double solnquad[2];
  double *fquad;
  double *fnode;
  double **ftild, **f;
  double *B, *C;
  double **V, **K;
  double solution[8]={0,0,0,0,0,0,0,0};
  double **rhs;
  double *innerrhs;
  double solnnodal;
  double **dubdn;
  double *z;
  double solndubdn[2];

  double *y;
  double val1, val2, yy;
  double maxerr;
  double l2error, s;
  int *pvt;
  int p, d, dd, i, j, k, info, nVtx, nPnls, nRhs = 1, test=0, dmax, pnlsX, pnlsY, Xidxv0, Xidxv1, Yidxv0, Yidxv1, idxv0, idxv1;
  char trans='N';
  dmax = 10; //max number of time steps
  nVtx = 5;//number of vertices
  dd = 0;
//build gauss legendre nodes
  p = 5; 
  x0[0] = 0.5;
//  x0[0] = 0.0;
  x0[1] = 0.0;

  /* parse the command line */
  for ( i=1; i<nargs; i++ )
    if ( argv[i][0] == '-' )
      switch ( argv[i][1] ) {
     case 't': test = atoi( argv[i]+3 );
        break;
      case 'p': p = atoi( argv[i]+3 ); //number of nodes
        break;
     case 'N': nVtx = atoi( argv[i]+3 ); //max # of vertices
        break;
     case 'M': dmax = atoi( argv[i]+3 ); //max # of time steps
        break;
     case 'T': tmax = atoi( argv[i]+3 ); //time step size numerator
        break;
     case 'd': dd = atoi( argv[i]+3 ); //which matrix to print out
        break;
      }
  nPnls = nVtx; //number of panels, -1 if panels dont wrap around
  h = tmax/dmax; //time step size
  printf("N = %d  M = %d  T = %f  p = %d \n",nVtx, dmax, tmax, p);

  xLeg = (double*)malloc( p*sizeof(double) );
  wLeg = (double*)malloc( p*sizeof(double) );
  Jacobi( p, 0.0, 0.0, xLeg, wLeg);

//build vertices
  vtx = (double **)calloc(nVtx, sizeof(double *));
  for(k = 0; k<nVtx; k++){
    vtx[k] = (double *)calloc(2, sizeof(double));
  }


#if SMOOTH
//circular geometry
  for(k = 0; k < nVtx; k++){
    theta = 2*M_PI*k/nVtx;
    paramcircle(theta, vtx[k]);
//    printf("vtx[%d][0] = %f  vtx[%d][1] = %f\n",k,vtx[k][0],k,vtx[k][1]);
  }
printf("circular geometry\n");
#else
//square geometry note: nVtx must be divisible by 4
  int vtxpside;
  vtxpside = nVtx/4;
  for(i = 1; i < 5; i ++){
    for(k = 0; k < vtxpside; k ++){
      theta = k/(double)(vtxpside);
      paramsquare(theta, vtx[k+(i-1)*vtxpside], i);
//      printf("vtx[%d][0] = %f  vtx[%d][1] = %f\n",k+(i-1)*vtxpside,vtx[k+(i-1)*vtxpside][0],k+(i-1)*vtxpside,vtx[k+(i-1)*vtxpside][1]);
    }
  }
printf("square geometry\n");
#endif

/*
//triangular geometry
  int vtxpside;
  vtxpside = nVtx/3;
  for(i = 1; i < 4; i ++){
    for(k = 0; k < vtxpside; k ++){
      theta = k/(double)(vtxpside);
      paramtriangle(theta, vtx[k+(i-1)*vtxpside], i);
      printf("vtx[%d][0] = %f  vtx[%d][1] = %f\n",k+(i-1)*vtxpside,vtx[k+(i-1)*vtxpside][0],k+(i-1)*vtxpside,vtx[k+(i-1)*vtxpside][1]);
    }
  }
printf("triangular geometry\n");


//elliptic geometry
  for(k = 0; k < nVtx; k++){
    theta = 2*M_PI*k/nVtx;
    paramellipse(theta, vtx[k]);
    printf("vtx[%d][0] = %f  vtx[%d][1] = %f\n",k,vtx[k][0],k,vtx[k][1]);
  }
*/
/*
for(i = 0; i < nVtx; i ++){
  printf("vtx[%d][0] = %f  vtx[%d][1] = %f\n",i,vtx[i][0],i,vtx[i][1]);
}
*/
//build panels
  pnls = (panel *)malloc(nPnls*sizeof(panel));
  for(k=0;k < nPnls; k++){
    pnls[k].v1 = vtx[k];
    pnls[k].idxv1 = k;  
    pnls[k].v0 = vtx[(k+1)%nPnls];
    pnls[k].idxv0 = (k+1)%nPnls; 
//printf("pnls[%d].v0[0]=%.3f  .v0[1]=%.3f  .v1[0]=%.3f .v1[1]=%.3f\n",k,pnls[k].v0[0], pnls[k].v0[1],pnls[k].v1[0],pnls[k].v1[1]);
    xdiff = vtx[k][0]-vtx[(k+1)%nPnls][0];
    ydiff = vtx[k][1]-vtx[(k+1)%nPnls][1];
    norm = sqrt(xdiff*xdiff + ydiff*ydiff);
    pnls[k].len = norm;
    t0 = xdiff/norm;
    t1 = ydiff/norm; 
    pnls[k].ny = (double *)malloc(2*sizeof(double));
    pnls[k].ny[0] = -t1; 
    pnls[k].ny[1] = t0; 
//    printf("%d nrm = [%lf, %lf]\n", k, pnls[k].ny[0], pnls[k].ny[1]);
  }

//build ftilde
  ftild = (double **)calloc(dmax+1, sizeof(double*));
  for(d = 0; d <= dmax; d++){ //loop over time steps
    fquad = (double *)calloc(nVtx, sizeof(double));
    for(k = 0; k < nPnls; k++){ //loop over panels (space)
      idxv0 = pnls[k].idxv0;  
      idxv1 = pnls[k].idxv1;  
      getfquad(h, d, p, &pnls[k], xLeg, wLeg, solnquad);
      fquad[k] += solnquad[0]+solnquad[1];
    }
    for(i = 0; i < nVtx; i++){
      //fquad[i] *= 0.5;
      //printf("fquad = %f\n",fquad[i]);
    }
    ftild[d] = fquad;
/*
      for(k = 0; k < nVtx; k++){
        printf("ftild[%d][%d] = %f\n",d,k,ftild[d][k]);
      }
*/
  }

//build f
  f = (double **)calloc(dmax+1, sizeof(double*));
  for(d = 0; d <= dmax; d++){
    fnode = (double *)calloc(nVtx, sizeof(double)); 
    for(j = 0; j < nVtx; j++){
      solnnodal = getfnode(h, d, p, vtx[j]);
      fnode[j] = solnnodal;
//      printf("fnode = %f\n",fnode[j]);
    }
    f[d] = fnode;
/*
    for(k = 0; k < nVtx; k++){
      printf("f[%d][%d] = %f\n",d,k,f[d][k]);
    }
*/
  }


//build the matrices V, K
  V = (double **)calloc(dmax+1, sizeof(double*));//dmax number of timesteps 
  K = (double **)calloc(dmax+1, sizeof(double*));

//generate B,C

  for(d = 0; d <= dmax; d++){
    B = (double *)calloc(nPnls*nPnls, sizeof(double));
    C = (double *)calloc(nVtx*nPnls, sizeof(double));
    for(pnlsX = 0; pnlsX < nPnls; pnlsX++){
      for(pnlsY = 0; pnlsY < nPnls; pnlsY++){
        Xidxv0 = pnls[pnlsX].idxv0; 
        Xidxv1 = pnls[pnlsX].idxv1;
        Yidxv0 = pnls[pnlsY].idxv0;
        Yidxv1 = pnls[pnlsY].idxv1;
        panelIA(d, p, h, &pnls[pnlsX], &pnls[pnlsY], xLeg, wLeg, solution);
        B[pnlsX + nPnls*pnlsY] = solution[0]+solution[1]+solution[2]+solution[3];
        C[pnlsX + nPnls*Yidxv0] += solution[4]+solution[6];
        C[pnlsX + nPnls*Yidxv1] += solution[5]+solution[7];
//        printf("single layer: %f %f %f %f\ndouble layer: %f %f %f %f\n",solution[0],solution[1],solution[2],solution[3],solution[4],solution[5],solution[6],solution[7]);
      }
    }
    V[d] = B;
    K[d] = C;
  }


//print out matrices
/*
for(dd = 0; dd <= dmax; dd++){
  printf("----matrix V[%d]------\n",dd);
  for(i = 0; i < nVtx; i++){
    for(j = 0; j < nVtx; j++){
      printf("%f ",V[dd][i*nVtx+j]);
    } printf("\n");
  }
  printf("----matrix K[%d]------\n",dd);
  for(i = 0; i < nVtx; i++){
    for(j = 0; j < nVtx; j++){
      printf("%.11f ",K[dd][i*nVtx+j]);
    } printf("\n");
  }
}
*/
//full rhs
  rhs = (double **)calloc(dmax+1, sizeof(double*));
  for(d = 0; d <= dmax; d++){ //loop over time steps
    innerrhs = (double *)calloc(nVtx, sizeof(double));
    rhs[d] = innerrhs;
  }

//allocate space for pivots
  pvt  = (int*)calloc(nVtx, sizeof(int));

//LU factorize A[0]
  dgetrf_(&nVtx, &nVtx, V[0], &nVtx, pvt, &info); 


printf("rhs\n");
//solve the system
  for(d = 0; d <= dmax; d++){
    //history
    for(k = 0; k < nVtx; k++){
      rhs[d][k] = -.5*ftild[d][k]; 
//      printf("solnrhsf[%d][%d] = %f \n",d,k,rhs[d][k]);
    } 
    for ( j=0; j < d; j++ ){ 
      MtVsub( nPnls, nVtx, V[d-j], rhs[j], rhs[d]);
/*
      for(k = 0; k < nVtx; k++){
        printf("solnrhs[%d][%d] = %f\n",d,k,rhs[j][k]);
      } 
*/
    }

/*
    for(k = 0; k < nVtx; k++){
      printf("solnrhsV[%d][%d] = %f \n",d,k,rhs[d][k]);
    } 
*/

    for ( j = 0; j < d; j++){
      MtVadd( nPnls, nVtx, K[d-j], f[j], rhs[d]);
    }
    MtVadd( nPnls, nVtx, K[0], f[d], rhs[d]);

/*
    for(k = 0; k < nVtx; k++){
      printf("solnrhsK[%d][%d] = %f \n",d,k,rhs[d][k]);
    } 
*/
    dgetrs_( &trans, &nVtx, &nRhs, V[0], &nVtx, pvt, rhs[d], &nVtx, &info); 

//print out solution

    printf("rhs[%d]\n",d);
    for(k = 0; k < nVtx; k++){
      printf("%.16f \n",rhs[d][k]); 
    }

  }

//error calculations

//max error
printf("dubdn\n");
  dubdn = (double **)calloc(dmax+1, sizeof(double*));
  for(d = 0; d <= dmax; d++){
    z = (double *)calloc(nVtx, sizeof(double));
    for(j = 0; j < nVtx; j++){
      getdubdn(h, d, &pnls[j], solndubdn);
      z[j] += 0.5*(solndubdn[0]+solndubdn[1]);
    }
    dubdn[d] = z;

    printf("dubdn[%d]\n",d);
    for(j = 0; j < nVtx; j++){
       printf("%.10f \n",dubdn[d][j]);
    }

  }

  double err[dmax+1];
  for(d = 0; d <= dmax; d++){
    y = (double *)calloc(nVtx, sizeof(double));
    for(j = 0; j < nVtx; j++){
      val1 = rhs[d][j];
      val2 = dubdn[d][j];
      y[j] = fabs(rhs[d][j]-dubdn[d][j]);
      printf("dubdn[%d][%d]: %.10f rhs=%.10f err=%lg\n",d,j,dubdn[d][j], rhs[d][j], y[j]);
    }
    err[d] = max(y,nVtx);
  }

  maxerr = max(err,dmax+1);
  printf("maxerr = %.10f\n",maxerr);

  //l2 error
  l2error = 0.0;
  for(d = 0; d <= dmax; d++){
    for(i = 0; i < nPnls; i++){
      s = l2err(h, d, p, &pnls[i], rhs[d][i], xLeg, wLeg);
      l2error += s;
    }
  }
  printf("l2 error = %.10f\n",sqrt(l2error));


//free allocated variables
  free(y);
  free(z);
  free(dubdn);
  free(B);
  free(C);
  free(V);
  free(K);
  free(ftild);
  free(f);
  free(rhs);
  free(fquad);
  free(fnode);
  free(pnls); 
  //free(vtx);

}

/*
[-0.1934349516453270,-0.4847693538752217,-0.1842260020001012 ,-0.1488078977100917,-0.0892623873391871,-0.0699365335331945,-0.0523916634114661,-0.0424350769565085,-0.0346014473903442,-0.0289864087608800,-0.0245735330322142]
[-0.1735143359,-0.5404320864,-0.3789424698,-0.2572772279,-0.1824097469,-0.1350877036,-0.1037253068,-0.0820087942,-0.0664012505,-0.0548292807,-0.0460220261]
*/

///////////////////////////
//subroutines
///////////////////////////

void paramcircle(double theta, double *x){
  x[0] = cos(theta);
  x[1] = sin(theta);
}

void paramsquare(double theta, double *x, int i){
  if(i==1){ //side 1
    x[0] = 1*(1-theta) + 1*theta;
    x[1] = -1*(1-theta) + 1*theta;
  }else if(i==2){ //side 2
    x[0] = 1*(1-theta) + -1*theta;
    x[1] = 1*(1-theta) + 1*theta;
  }else if(i==3){ //side 3
    x[0] = -1*(1-theta) + -1*theta;
    x[1] = 1*(1-theta) + -1*theta;
  }else if(i==4){ //side 4
    x[0] = -1*(1-theta) + 1*theta;
    x[1] = -1*(1-theta) + -1*theta;
  }
}

void paramtriangle(double theta, double *x, int i){
  if(i==1){ //side 1
    x[0] = 1*(1-theta) + 0*theta;
    x[1] = -1*(1-theta) + 1*theta;
  }else if(i==2){ //side 2
    x[0] = 0*(1-theta) + -1*theta;
    x[1] = 1*(1-theta) + -1*theta;
  }else if(i==3){ //side 3
    x[0] = -1*(1-theta) + 1*theta;
    x[1] = -1*(1-theta) + -1*theta;
  }
}

void paramellipse(double theta, double *x){
  x[0] = cos(theta);
  x[1] = 10*sin(theta);
}

// right hand side
double f( double *x, double t) {
  return t;
} /* fcn */	

double f2(double *x, double t) {
  double r1, r2, r, soln;
  r1 = x[0]-x0[0];
  r2 = x[1]-x0[1];
  r = r1*r1 + r2*r2;
  soln = exp(-r/(4*t))/(4*M_PI*t);
  return soln;
  }
  
void getfquad(double h, int i, int p, panel *pX, double xLeg[], double wLeg[], double solnquad[]){
  int i1, i2;
  double fquad, x[2], tquad, length;
  solnquad[0] = solnquad[1] = 0.0;
  length  = (pX->v1[0] - pX->v0[0])*(pX->v1[0] - pX->v0[0]);
  length += (pX->v1[1] - pX->v0[1])*(pX->v1[1] - pX->v0[1]);
  length = sqrt(length);

  for(i1 = 0; i1 < p; i1++){ //time loop
    tquad = h*(i + xLeg[i1]);   
    for(i2 = 0; i2 < p; i2++){ //space loop
      x[0] =  xLeg[i2]*pX->v0[0] + (1-xLeg[i2])*pX->v1[0];
      x[1] =  xLeg[i2]*pX->v0[1] + (1-xLeg[i2])*pX->v1[1];  
      //fcn = f(x,t);   
      fquad = f2(x,tquad); 
      fquad *= wLeg[i1]*wLeg[i2];    
      solnquad[0] += fquad*xLeg[i2];
      solnquad[1] += fquad*(1 - xLeg[i2]);
    }
  }
  solnquad[0] *= h*length;
  solnquad[1] *= h*length;
//  printf("solnquad[0] = %f solutionquad[1] = %f \n",solnquad[0],solnquad[1]);
}

double l2err(double h, int i, int p, panel *pX, double rhs, double xLeg[], double wLeg[]){
  int i1, i2;
  double t, integ2, integ, x[2], a[2], len;
  integ = 0.0;

  a[0] = pX->v1[0] - pX->v0[0];
  a[1] = pX->v1[1] - pX->v0[1];
  len = sqrt( a[0]*a[0] + a[1]*a[1] );

  for(i1 = 0; i1 < p; i1++){//time loop
    t = h*(i + xLeg[i1]);
    for(i2 = 0; i2 < p; i2++){//space loop
      x[0] =  xLeg[i2]*pX->v0[0] + (1-xLeg[i2])*pX->v1[0];
      x[1] =  xLeg[i2]*pX->v0[1] + (1-xLeg[i2])*pX->v1[1];  
      integ2 = g(x, t, pX->ny) - rhs;
      integ += wLeg[i1]*wLeg[i2]*integ2*integ2;
    }
  }
  return h*len*integ;
}


double getfnode(double h, int i, int p, double *x){
  int i1, i2;
  double tnodal, solnnodal;
  tnodal = h*(i + 0.5);
  solnnodal = f2(x,tnodal); 
//  printf("solnnodal = %f\n",solnnodal);
  return solnnodal;
}

/* Johannes Tausch
 * Computes the matrix-vector product
 *    y += A*x
 * where A is an n-by-n matrix and adds the result to the vector y
 */
// modified to take in rectangular matrices
void MtVadd( int n, int m, double *A, double *x, double *y){
  int i, j;
//rows n, columns m
  for ( i=0; i<m; i++ ) {
    for ( j=0; j<n; j++ ) {
      y[i] += A[i+n*j]*x[j];
    }
  }
} /* MtVadd */

 /* Johannes Tausch
 * Computes the matrix-vector product
 *    y -= A*x
 * where A is an n-by-n matrix and subtracts the result from the vector y
 */
// modified to take in rectangular matrices
void MtVsub( int n, int m, double *A, double *x, double *y ){
  int i, j;
//rows n, columns m
  for ( i=0; i<m; i++ ) {
    for ( j=0; j<n; j++ ) {
      y[i] -= A[i+n*j]*x[j];
    }
  }
} /* MtVsub */

double g(double x[], double t, double ny[]) {
  double r0, r1, r, rny, soln;
  r0 = x[0]-x0[0];
  r1 = x[1]-x0[1];
  r = r0*r0 + r1*r1;
  rny = r0*ny[0] + r1*ny[1];
  soln = -exp(-r/(4*t))*rny/(4*M_PI*2*t*t);
  return soln;
  }
/*
  fcnVal[0] = exp(-rr/(4*del))/(del*4*M_PI);
  fcnVal[1] = drr*fcnVal[0]/(2*del);
*/

void getdubdn(double h, int i, panel *p, double dubdnsoln[]){
  double t;
  t = h*(i + 0.5);
  dubdnsoln[0] = 0.0; dubdnsoln[1] = 0.0;
  dubdnsoln[0] = g(p->v0,t, p->ny); 
  dubdnsoln[1] = g(p->v1,t, p->ny);
}

/*http://stackoverflow.com/questions/1690428/finding-max-number-in-an-array-c-programming
*/
double max(double myArray[], int size) {
    double maxValue = myArray[0];
    int i;

    for (i = 1; i < size; ++i) {
        if ( myArray[i] > maxValue ) {
            maxValue = myArray[i];
        }
    }
    return maxValue;
}



