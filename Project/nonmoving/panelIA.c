//compile object file:
//gcc -c panelIA.c

#include <math.h>
#include <stdio.h>
#include "panelIA.h"

#define TOL 1e-12

int nrcommonvtx(panel *p1, panel *p2);
void intgrnonsing(int d, int p, double h, panel *p1, panel *p2, double a[], double b[], double xLeg[], double wLeg[], double solution[]);
void intgrd1nc1(int d, int p, double h, panel *p1, panel *p2, double a[], double b[], double xLeg[], double wLeg[], double solution[]);
void intgrd1nc2(int d, int p, double h, panel *p1, panel *p2, double a[], double b[], double xLeg[], double wLeg[], double solution[]);
void intgrd0nc1(int d, int p, double h, panel *p1, panel *p2, double a[], double b[], double xLeg[], double wLeg[], double solution[]);
void intgrd0nc2(int d, int p, double h, panel *p1, panel *p2, double a[], double b[], double xLeg[], double wLeg[], double solution[]);
void kernel(double r[], double del, double *ny, double fcnVal[]);

void panelIA(int d, int p, double h, panel *p1, panel *p2, double xLeg[], double wLeg[], double solution[]) {
  double *v0, *v1,  a[2], *w0, *w1, b[2];
  int nc;
  solution[0] = solution[1] = solution[2] = solution[3] = solution[4] = solution[5] = solution[6] = solution[7] = 0.0;
  nc = nrcommonvtx(p1,p2);
  v0 = p1->v0;
  v1 = p1->v1;
  w0 = p2->v0;
  w1 = p2->v1;
  if(d>1){
    a[0] = v1[0]-v0[0];
    a[1] = v1[1]-v0[1];
    b[0] = w1[0]-w0[0];
    b[1] = w1[1]-w0[1];
    intgrnonsing(d, p, h, p1, p2, a, b, xLeg, wLeg, solution);
  }else if(d==1){
    if(nc==2){//p and q are the same
      a[0] = v0[0]-v1[0];
      a[1] = v0[1]-v1[1];
      b[0] = w0[0]-w1[0];
      b[1] = w0[1]-w1[1];
      intgrd1nc2(d, p, h, p1, p2, a, b, xLeg, wLeg, solution);
    } else if(nc==1) {
      a[0] = v1[0]-v0[0];
      a[1] = v1[1]-v0[1];
      b[0] = w0[0]-w1[0];
      b[1] = w0[1]-w1[1];
      intgrd1nc1(d, p, h, p1, p2, a, b, xLeg, wLeg, solution);
    } else if(nc==-1) {
      a[0] = v0[0]-v1[0];
      a[1] = v0[1]-v1[1];
      b[0] = w1[0]-w0[0];
      b[1] = w1[1]-w0[1];
      intgrd1nc1(d, p, h, p1, p2, a, b, xLeg, wLeg, solution);
    } else if(nc==0){
      a[0] = v0[0]-v1[0];
      a[1] = v0[1]-v1[1];
      b[0] = w0[0]-w1[0];
      b[1] = w0[1]-w1[1];   
      intgrnonsing(d, p, h, p1, p2, a, b, xLeg, wLeg, solution);
    }
  }else if(d==0){
    if(nc==2){
      a[0] = v0[0]-v1[0];
      a[1] = v0[1]-v1[1];
      b[0] = w0[0]-w1[0];
      b[1] = w0[1]-w1[1];   
      intgrd0nc2(d, p, h, p1, p2, a, b, xLeg, wLeg, solution);      
    } else if(nc==1) {
      a[0] = v1[0]-v0[0];
      a[1] = v1[1]-v0[1];
      b[0] = w0[0]-w1[0];
      b[1] = w0[1]-w1[1];
      intgrd0nc1(d, p, h, p1, p2, a, b, xLeg, wLeg, solution);
    } else if(nc==-1) {
      a[0] = v0[0]-v1[0];
      a[1] = v0[1]-v1[1];
      b[0] = w1[0]-w0[0];
      b[1] = w1[1]-w0[1];
      intgrd0nc1(d, p, h, p1, p2, a, b, xLeg, wLeg, solution);
    } else if(nc==0){
      a[0] = v0[0]-v1[0];
      a[1] = v0[1]-v1[1];
      b[0] = w0[0]-w1[0];
      b[1] = w0[1]-w1[1];   
      intgrnonsing(d, p, h, p1, p2, a, b, xLeg, wLeg, solution);

    }    
  }
}

int nrcommonvtx(panel *p1, panel *p2) {
  if (p1->v0==p2->v0){

    return 2;//panels are the same
    
  }else if (p1->v0==p2->v1) {

    return 1;//panels share 1 vertex, pX on left, pY on right
    
  }else if (p1->v1==p2->v0) {

    return -1;//panels share 1 vertex, pX on right, pY on left

  }else {

    return 0;//panels are not neighbors

  }
}


void intgrnonsing(int d, int p, double h, panel *p1, panel *p2, double a[], double b[], double xLeg[], double wLeg[], double solution[]) {
  double *xhat, *yhat, *that, *tauhat, *v0, *v1, *w0, *w1, *ny, r[2], del, jac, fcnVal[2], dfcn, fcn, norma, normb;
  int i1, i2, i3, i4;
  v0 = p1->v0;
  v1 = p1->v1;
  w0 = p2->v0;
  w1 = p2->v1;
  ny = p2->ny;
  xhat = xLeg;
  yhat = xLeg;
  tauhat = xLeg;
  that = xLeg; 
  norma = p1->len;
  normb = p2->len;
  for(i1 = 0; i1 < p; i1++) { //xhat loop
    for(i2 = 0; i2 < p; i2++) { //yhat loop
      for(i3 = 0; i3 < p; i3++) { //that loop
	for(i4 = 0; i4 < p; i4++) { //tauhat loop
	  del = h*(d + that[i3] - tauhat[i4]);
	  r[0] = v1[0]*(1 - xhat[i1]) + v0[0]*xhat[i1] - 
	    w1[0]*(1 - yhat[i2]) - w0[0]*yhat[i2];
	  r[1] = v1[1]*(1 - xhat[i1]) + v0[1]*xhat[i1] -
	    w1[1]*(1 - yhat[i2]) - w0[1]*yhat[i2];
          jac = h*h*norma*normb;
	  //heat kernel
          if(del < TOL){
            fcn = 0.0;
            dfcn = 0.0;
          }
          else{
            kernel(r, del, ny, fcnVal);
	    fcn = fcnVal[0]*jac;
            dfcn = fcnVal[1]*jac;
          }

	  //weights
	  fcn *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
     
	  //test functions
	  solution[0] += fcn*yhat[i2]*xhat[i1];
	  solution[1] += fcn*(1 - yhat[i2])*xhat[i1];
	  solution[2] += fcn*yhat[i2]*(1 - xhat[i1]);
	  solution[3] += fcn*(1 - yhat[i2])*(1 - xhat[i1]);

	  solution[4] += dfcn*yhat[i2]*xhat[i1];
	  solution[5] += dfcn*(1 - yhat[i2])*xhat[i1];
	  solution[6] += dfcn*yhat[i2]*(1 - xhat[i1]);
	  solution[7] += dfcn*(1 - yhat[i2])*(1 - xhat[i1]);

	}
      }
    }
  }
}

void intgrd1nc1(int d, int p, double h, panel *p1, panel *p2, double a[], double b[], double xLeg[], double wLeg[], double solution[]) {
  double *eta, *xi1, *xi2, *xi3, *omega, *that, *tauhat, *ny, *v0, *v1, *w0, *w1, norma, normb, rho, etasq, xi1sq, xi2sq, del1, del2, del3, fcn1, fcn2, fcn3, fcn4, fcn, eta3, r1[2], r2[2], r3[2], rr1, rr2, rr3, r[2], fcnVal[2], dfcn1, dfcn2, dfcn3, dfcn4, rr, xhat1, yhat1, xhat2, yhat2, xhat3, yhat3, jac;
  int i1, i2, i3, i4;
  v0 = p1->v0;
  v1 = p1->v1;
  w0 = p2->v0;
  w1 = p2->v1;
  ny = p2->ny;
  eta = xLeg;
  xi1 = xLeg;
  xi2 = xLeg;
  xi3 = xLeg;
  omega = xLeg; 
  norma = p1->len;
  normb = p2->len;
  rho = norma*norma/(4*h);
  for(i1 = 0; i1 < p; i1++) { //xi1 loop
    for(i2 = 0; i2 < p; i2++) { //xi2 loop
      for(i3 = 0; i3 < p; i3++) { //xi3 loop
	for(i4 = 0; i4 < p; i4++) { //eta loop
	  eta3 = eta[i4]*eta[i4]*eta[i4];
	  del1 = h*eta[i4]*eta[i4]*(xi2[i2]*xi2[i2] + xi3[i3]*xi3[i3]);
	  del2 = h*eta[i4]*eta[i4]*(1 + xi3[i3]*xi3[i3]);
	  r1[0] = eta[i4]*a[0] - eta[i4]*xi1[i1]*b[0];
	  r1[1] = eta[i4]*a[1] - eta[i4]*xi1[i1]*b[1];
          jac = norma*normb*4*h*h*eta[i4]*xi2[i2]*eta[i4]*xi3[i3]*eta3;
          kernel(r1, del1, ny, fcnVal);
          fcn1 = fcnVal[0]*jac;
          dfcn1 = fcnVal[1]*jac;

	  r2[0] = eta[i4]*xi1[i1]*a[0] - eta[i4]*b[0];
	  r2[1] = eta[i4]*xi1[i1]*a[1] - eta[i4]*b[1];
          jac = norma*normb*4*h*h*eta[i4]*xi2[i2]*eta[i4]*xi3[i3]*eta3;
          kernel(r2, del1, ny, fcnVal);
          fcn2 = fcnVal[0]*jac;
          dfcn2 = fcnVal[1]*jac;

	  r3[0] = eta[i4]*xi1[i1]*a[0] - eta[i4]*xi2[i2]*b[0];
	  r3[1] = eta[i4]*xi1[i1]*a[1] - eta[i4]*xi2[i2]*b[1];
          jac = norma*normb*4*h*h*eta[i4]*eta[i4]*xi3[i3]*eta3;
          kernel(r3, del2, ny, fcnVal);
	  fcn3 = fcnVal[0]*jac;
          dfcn3 = fcnVal[1]*jac;

          jac = norma*normb*4*h*h*eta[i4]*xi3[i3]*eta[i4]*eta3;
          kernel(r3, del2, ny, fcnVal);
	  fcn4 = fcnVal[0]*jac;
          dfcn4 = fcnVal[1]*jac;

	  //weights, jacobian 
	  fcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  fcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  fcn3 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  fcn4 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn3 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn4 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];

	  //test functions
          xhat1 = eta[i4];
          yhat1 = eta[i4]*xi1[i1];
          xhat2 = eta[i4]*xi1[i1];
          yhat2 = eta[i4];
          xhat3 = eta[i4]*xi1[i1];
          yhat3 = eta[i4]*xi2[i2];

	  solution[0] += fcn1*yhat1*xhat1 + fcn2*yhat2*xhat2 + fcn3*yhat3*xhat3 + fcn4*yhat3*xhat3;
	  solution[1] += fcn1*(1 - yhat1)*xhat1 + fcn2*(1 - yhat2)*xhat2 + fcn3*(1 - yhat3)*xhat3 + fcn4*(1 - yhat3)*xhat3;
	  solution[2] += fcn1*yhat1*(1 - xhat1) + fcn2*yhat2*(1 - xhat2) + fcn3*yhat3*(1 - xhat3) + fcn4*yhat3*(1 - xhat3);
	  solution[3] += fcn1*(1 - yhat1)*(1 - xhat1) + fcn2*(1 - yhat2)*(1 - xhat2) + fcn3*(1 - yhat3)*(1 - xhat3)+ fcn4*(1 - yhat3)*(1 - xhat3); 

	  solution[4] += dfcn1*yhat1*xhat1 + dfcn2*yhat2*xhat2 + dfcn3*yhat3*xhat3 + dfcn4*yhat3*xhat3;
	  solution[5] += dfcn1*(1 - yhat1)*xhat1 + dfcn2*(1 - yhat2)*xhat2 + dfcn3*(1 - yhat3)*xhat3 + dfcn4*(1 - yhat3)*xhat3;
	  solution[6] += dfcn1*yhat1*(1 - xhat1) + dfcn2*yhat2*(1 - xhat2) + dfcn3*yhat3*(1 - xhat3) + dfcn4*yhat3*(1 - xhat3);
	  solution[7] += dfcn1*(1 - yhat1)*(1 - xhat1) + dfcn2*(1 - yhat2)*(1 - xhat2) + dfcn3*(1 - yhat3)*(1 - xhat3)+ dfcn4*(1 - yhat3)*(1 - xhat3); 	
	}
      }
    }
  }
}

void intgrd1nc2(int d, int p, double h,panel *p1, panel *p2, double a[], double b[], double xLeg[], double wLeg[], double solution[]) {
  double *eta, *xi1, *xi2, *omega, *ny, *v0, *v1, *w0, *w1, norma, normb, rho, etasq, xi1sq, xi2sq, del, fcn1, fcn2, fcn3, xhat1a, xhat1b, xhat2a, xhat2b, xhat3a, xhat3b, xhat9b, yhat1a, yhat1b, yhat2a, yhat2b, yhat3a, yhat3b, fcnVal[2], r[2], jac, dfcn1, dfcn1b, dfcn2, dfcn2b, dfcn3, dfcn3b, fcn1b, fcn2b, fcn3b;
  int i1, i2, i3, i4;
  v0 = p1->v0;
  v1 = p1->v1;
  w0 = p2->v0;
  w1 = p2->v1;
  ny = p2->ny;
  eta = xLeg;
  xi1 = xLeg;
  xi2 = xLeg;
  omega = xLeg; 
  norma = p1->len;
  normb = p2->len;
  rho = norma*norma/(4*h);
  for(i1 = 0; i1 < p; i1++) { //xi1 loop
    for(i2 = 0; i2 < p; i2++) { //xi2 loop
      for(i3 = 0; i3 < p; i3++) { //eta loop
	for(i4 = 0; i4 < p; i4++) { //omega loop
          etasq = eta[i3]*eta[i3];
          xi1sq = xi1[i1]*xi1[i1];
          xi2sq = xi2[i2]*xi2[i2];

          r[0] = a[0]*eta[i3]*xi2[i2];
          r[1] = a[1]*eta[i3]*xi2[i2];
          del = h*etasq*(1 + xi1sq);
          jac = norma*norma*h*h*4*eta[i3]*eta[i3]*xi1[i1]*etasq*(1-eta[i3]*xi2[i2]);
          kernel(r, del, ny, fcnVal);
          fcn1 = fcnVal[0]*jac;
          dfcn1 = fcnVal[1]*jac;

          r[0] = -a[0]*eta[i3]*xi2[i2];
          r[1] = -a[1]*eta[i3]*xi2[i2];
          kernel(r, del, ny, fcnVal);
          fcn1b = fcnVal[0]*jac;
          dfcn1b = fcnVal[1]*jac;

          r[0] = a[0]*eta[i3]*xi2[i2];
          r[1] = a[1]*eta[i3]*xi2[i2];
          del = h*etasq*(xi1sq + 1);
          jac = norma*norma*h*h*4*eta[i3]*xi1[i1]*eta[i3]*etasq*(1-eta[i3]*xi2[i2]);          
          kernel(r, del, ny, fcnVal);
          fcn2 = fcnVal[0]*jac;
          dfcn2 = fcnVal[1]*jac;

          r[0] = -a[0]*eta[i3]*xi2[i2];
          r[1] = -a[1]*eta[i3]*xi2[i2];
          kernel(r, del, ny, fcnVal);
          fcn2b = fcnVal[0]*jac;
          dfcn2b = fcnVal[1]*jac;

          r[0] = a[0]*eta[i3];
          r[1] = a[1]*eta[i3];
          del = h*etasq*(xi1sq + xi2sq);
          jac = norma*norma*h*h*4*eta[i4]*xi1[i1]*eta[i3]*xi2[i2]*etasq*(1-eta[i3]);
          kernel(r, del, ny, fcnVal);
          fcn3 = fcnVal[0]*jac;
          dfcn3 = fcnVal[1]*jac;

          r[0] = -a[0]*eta[i3];
          r[1] = -a[1]*eta[i3];
          kernel(r, del, ny, fcnVal);
          fcn3b = fcnVal[0]*jac;
          dfcn3b = fcnVal[1]*jac;

          //weights, jacobian 
	  fcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn3 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  fcn1b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn2b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn3b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn3 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn1b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn2b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn3b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];

	  //test functions
          xhat1a = eta[i3]*xi2[i2] + omega[i4]*(1-eta[i3]*xi2[i2]);
          yhat1a = omega[i4]*(1-eta[i3]*xi2[i2]);
          xhat1b = omega[i4]*(1-eta[i3]*xi2[i2]);
          yhat1b = eta[i3]*xi2[i2] + omega[i4]*(1-eta[i3]*xi2[i2]);
          xhat2a = xhat1a;
          yhat2a = xhat1a;
          xhat2b = xhat1b;
          yhat2b = yhat1b;
          xhat3a = eta[i3] + omega[i4]*(1-eta[i3]);
          yhat3a = omega[i4]*(1-eta[i3]);
          xhat3b = omega[i4]*(1-eta[i3]);
          yhat3b = eta[i3] + omega[i4]*(1-eta[i3]);


	  solution[0] += fcn1*yhat1a*xhat1a + fcn1b*yhat1b*xhat1b + fcn2*yhat2a*xhat2a + fcn2b*yhat2b*xhat2b + fcn3*yhat3a*xhat3a + fcn3b*yhat3b*xhat3b;

	  solution[1] += fcn1*(1 - yhat1a)*xhat1a + fcn1b*(1 - yhat1b)*xhat1b + fcn2*(1 - yhat2a)*xhat2a + fcn2b*(1 - yhat2b)*xhat2b + fcn3*(1 - yhat3a)*xhat3a + fcn3b*(1 - yhat3b)*xhat3b;

	  solution[2] += fcn1*yhat1a*(1 - xhat1a) + fcn1b*yhat1b*(1 - xhat1b) + fcn2*yhat2a*(1 - xhat2a) + fcn2b*yhat2b*(1 - xhat2b) + fcn3*yhat3a*(1 - xhat3a) + fcn3b*yhat3b*(1 - xhat3b);

	  solution[3] += fcn1*(1 - yhat1a)*(1 - xhat1a) + fcn1b*(1 - yhat1b)*(1 - xhat1b) + fcn2*(1 - yhat2a)*(1 - xhat2a) + fcn2b*(1 - yhat2b)*(1 - xhat2b) + fcn3*(1 - yhat3a)*(1 - xhat3a) + fcn3b*(1 - yhat3b)*(1 - xhat3b); 


	  solution[4] += dfcn1*yhat1a*xhat1a + dfcn1b*yhat1b*xhat1b + dfcn2*yhat2a*xhat2a + dfcn2b*yhat2b*xhat2b + dfcn3*yhat3a*xhat3a + dfcn3b*yhat3b*xhat3b;

	  solution[5] += dfcn1*(1 - yhat1a)*xhat1a + dfcn1b*(1 - yhat1b)*xhat1b + dfcn2*(1 - yhat2a)*xhat2a + dfcn2b*(1 - yhat2b)*xhat2b + dfcn3*(1 - yhat3a)*xhat3a + dfcn3b*(1 - yhat3b)*xhat3b;

	  solution[6] += dfcn1*yhat1a*(1 - xhat1a) + dfcn1b*yhat1b*(1 - xhat1b) + dfcn2*yhat2a*(1 - xhat2a) + dfcn2b*yhat2b*(1 - xhat2b) + dfcn3*yhat3a*(1 - xhat3a) + dfcn3b*yhat3b*(1 - xhat3b);

	  solution[7] += dfcn1*(1 - yhat1a)*(1 - xhat1a) + dfcn1b*(1 - yhat1b)*(1 - xhat1b) + dfcn2*(1 - yhat2a)*(1 - xhat2a) + dfcn2b*(1 - yhat2b)*(1 - xhat2b) + dfcn3*(1 - yhat3a)*(1 - xhat3a) + dfcn3b*(1 - yhat3b)*(1 - xhat3b); 

	}
      }
    }
  }
}






void intgrd0nc1(int d, int p, double h, panel *p1, panel *p2, double a[], double b[], double xLeg[], double wLeg[], double solution[]) {
  double *omega, *eta, *xi1, *xi2, *ny, *v0, *v1, *w0, *w1,  etasq, xi2sq, r[2], fcnVal[2], del, jac, fcn1, fcn2, fcn3, xhat1, xhat2, xhat3, yhat1, yhat2, yhat3, norma, normb, dfcn1, dfcn2, dfcn3;
  int i1, i2, i3, i4;
  double delta;
  delta = 0.2947183339744; //1/(sqrt(ln(1/(10^-5)))
  v0 = p1->v0;
  v1 = p1->v1;
  w0 = p2->v0;
  w1 = p2->v1;
  ny = p2->ny;
  omega = xLeg;
  eta = xLeg;
  xi1 = xLeg;
  xi2 = xLeg; 
  norma = p1->len;
  normb = p2->len;
  for(i1 = 0; i1 < p; i1++) { //xi1 loop
    for(i2 = 0; i2 < p; i2++) { //xi2 loop
      for(i3 = 0; i3 < p; i3++) { //eta loop
	for(i4 = 0; i4 < p; i4++) { //sigma loop
          etasq = eta[i3]*eta[i3];
          xi2sq = xi2[i2]*xi2[i2];

          r[0] = eta[i3]*a[0] - eta[i3]*xi1[i1]*b[0];
          r[1] = eta[i3]*a[1] - eta[i3]*xi1[i1]*b[1];
          del = h*etasq*xi2sq;
          jac = norma*normb*h*h*2*eta[i3]*xi2[i2]*etasq*(1-etasq*xi2sq);
          kernel(r, del, ny, fcnVal);
          fcn1 = fcnVal[0]*jac;
          dfcn1 = fcnVal[1]*jac;
//dfcn1 = 0.0;
          r[0] = eta[i3]*xi1[i1]*a[0] - eta[i3]*b[0];
          r[1] = eta[i3]*xi1[i1]*a[1] - eta[i3]*b[1];
          del = h*etasq*xi2sq;
          jac = norma*normb*h*h*2*eta[i3]*xi2[i2]*etasq*(1-etasq*xi2sq);
          kernel(r, del, ny, fcnVal);
          fcn2 = fcnVal[0]*jac;
          dfcn2 = fcnVal[1]*jac;
//dfcn2=0.0;
          r[0] = eta[i3]*xi1[i1]*a[0] - eta[i3]*xi2[i2]*b[0];
          r[1] = eta[i3]*xi1[i1]*a[1] - eta[i3]*xi2[i2]*b[1];
          del = h*etasq;
          jac = norma*normb*h*h*2*eta[i3]*etasq*(1-etasq);
          kernel(r, del, ny, fcnVal);
          fcn3 = fcnVal[0]*jac;
          dfcn3 = fcnVal[1]*jac;
//dfcn3=0.0;
//          fcn2 = fcn3 = 0;
//          fcn1 = fcn3 = 0;
//          fcn1 = fcn2 = 0;
          //weights, jacobian 
	  fcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn3 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn3 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];

	  //test functions
          xhat1 = eta[i3];
          yhat1 = eta[i3]*xi1[i1];
          xhat2 = eta[i3]*xi1[i1];
          yhat2 = eta[i3];
          xhat3 = eta[i3]*xi1[i1];
          yhat3 = eta[i3]*xi2[i2];
	  solution[0] += fcn1*yhat1*xhat1 + fcn2*yhat2*xhat2 +fcn3*yhat3*xhat3;
	  solution[1] += fcn1*(1 - yhat1)*xhat1 + fcn2*(1 - yhat2)*xhat2 + fcn3*(1 - yhat3)*xhat3;
	  solution[2] += fcn1*yhat1*(1 - xhat1) + fcn2*yhat2*(1 - xhat2) +fcn3*yhat3*(1 - xhat3);
	  solution[3] += fcn1*(1 - yhat1)*(1 - xhat1) + fcn2*(1 - yhat2)*(1 - xhat2) + fcn3*(1 - yhat3)*(1 - xhat3); 

	  solution[4] += dfcn1*yhat1*xhat1 + dfcn2*yhat2*xhat2 +dfcn3*yhat3*xhat3;
	  solution[5] += dfcn1*(1 - yhat1)*xhat1 + dfcn2*(1 - yhat2)*xhat2 + dfcn3*(1 - yhat3)*xhat3;
	  solution[6] += dfcn1*yhat1*(1 - xhat1) + dfcn2*yhat2*(1 - xhat2) +dfcn3*yhat3*(1 - xhat3);
	  solution[7] += dfcn1*(1 - yhat1)*(1 - xhat1) + dfcn2*(1 - yhat2)*(1 - xhat2) + dfcn3*(1 - yhat3)*(1 - xhat3); 

	}
      }
    }
  }
}



void intgrd0nc2(int d, int p, double h, panel *p1, panel *p2, double a[], double b[], double xLeg[], double wLeg[], double solution[]) {
  double *omega, *delta, *eta, *xi, *v0, *v1, *w0, *w1, fcnVal[2], r[2], del, *ny, jac, fcn1, fcn1b, fcn2, fcn2b, dfcn1, dfcn1b, dfcn2, dfcn2b, xhat1a, xhat1b, xhat2a, xhat2b, yhat1a, yhat1b, yhat2a, yhat2b, norma, normb, xisq, etasq;
  int i1, i2, i3, i4;
  v0 = p1->v0;
  v1 = p1->v1;
  w0 = p2->v0;
  w1 = p2->v1;
  ny = p2->ny;
  omega = xLeg;
  delta = xLeg;
  eta = xLeg;
  xi = xLeg; 
  norma = p1->len;
  normb = p2->len;
  double rr, drr;
  for(i1 = 0; i1 < p; i1++) { //eta loop
    for(i2 = 0; i2 < p; i2++) { //xi loop
      for(i3 = 0; i3 < p; i3++) { //omega loop
	for(i4 = 0; i4 < p; i4++) { //delta loop
          etasq = eta[i1]*eta[i1];
          xisq = xi[i2]*xi[i2];

          r[0] = a[0]*eta[i1]*xi[i2];
          r[1] = a[1]*eta[i1]*xi[i2];
          del = h*etasq;
          jac = norma*norma*h*h*2*eta[i1]*(1-eta[i1]*xi[i2])*(1-etasq)*eta[i1];
          kernel(r, del, ny, fcnVal);
          fcn1 = fcnVal[0]*jac;
          dfcn1 = fcnVal[1]*jac;

          r[0] = -a[0]*eta[i1]*xi[i2];
          r[1] = -a[1]*eta[i1]*xi[i2];
          kernel(r, del, ny, fcnVal);
          fcn1b = fcnVal[0]*jac;
          dfcn1b = fcnVal[1]*jac;

          r[0] = a[0]*eta[i1];
          r[1] = a[1]*eta[i1];
          del = h*etasq*xisq;
          jac = norma*norma*h*h*2*eta[i1]*xi[i2]*(1-eta[i1])*(1-etasq*xisq)*eta[i1];
          kernel(r, del, ny, fcnVal);
          fcn2 = fcnVal[0]*jac;
          dfcn2 = fcnVal[1]*jac;

          r[0] = -a[0]*eta[i1];
          r[1] = -a[1]*eta[i1];
          kernel(r, del, ny, fcnVal);
          fcn2b = fcnVal[0]*jac;
          dfcn2b = fcnVal[1]*jac;

          //weights, jacobian 
	  fcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn1b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  fcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          fcn2b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn1 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn1b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
	  dfcn2 *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];
          dfcn2b *= wLeg[i1]*wLeg[i2]*wLeg[i3]*wLeg[i4];

	  //test functions
          xhat1a = eta[i1]*xi[i2] + omega[i3]*(1-eta[i1]*xi[i2]);
          yhat1a = omega[i3]*(1-eta[i1]*xi[i2]);
          xhat1b = omega[i3]*(1-eta[i1]*xi[i2]);
          yhat1b = eta[i1]*xi[i2] + omega[i3]*(1-eta[i1]*xi[i2]);

          xhat2a = eta[i1] + omega[i3]*(1-eta[i1]);
          yhat2a = omega[i3]*(1-eta[i1]);
          xhat2b = omega[i3]*(1-eta[i1]);
          yhat2b = eta[i1] + omega[i3]*(1-eta[i1]);


	  solution[0] += fcn1*yhat1a*xhat1a + fcn1b*yhat1b*xhat1b + fcn2*yhat2a*xhat2a + fcn2b*yhat2b*xhat2b;
	  solution[1] += fcn1*(1 - yhat1a)*xhat1a + fcn1b*(1 - yhat1b)*xhat1b + fcn2*(1 - yhat2a)*xhat2a + fcn2b*(1 - yhat2b)*xhat2b;
	  solution[2] += fcn1*yhat1a*(1 - xhat1a) + fcn1b*yhat1b*(1 - xhat1b) + fcn2*yhat2a*(1 - xhat2a) + fcn2b*yhat2b*(1 - xhat2b);
	  solution[3] += fcn1*(1 - yhat1a)*(1 - xhat1a) + fcn1b*(1 - yhat1b)*(1 - xhat1b) + fcn2*(1 - yhat2a)*(1 - xhat2a) + fcn2b*(1 - yhat2b)*(1 - xhat2b); 

	  solution[4] += dfcn1*yhat1a*xhat1a + dfcn1b*yhat1b*xhat1b + dfcn2*yhat2a*xhat2a + dfcn2b*yhat2b*xhat2b;
	  solution[5] += dfcn1*(1 - yhat1a)*xhat1a + dfcn1b*(1 - yhat1b)*xhat1b + dfcn2*(1 - yhat2a)*xhat2a + dfcn2b*(1 - yhat2b)*xhat2b;
	  solution[6] += dfcn1*yhat1a*(1 - xhat1a) + dfcn1b*yhat1b*(1 - xhat1b) + dfcn2*yhat2a*(1 - xhat2a) + dfcn2b*yhat2b*(1 - xhat2b);
	  solution[7] += dfcn1*(1 - yhat1a)*(1 - xhat1a) + dfcn1b*(1 - yhat1b)*(1 - xhat1b) + dfcn2*(1 - yhat2a)*(1 - xhat2a) + dfcn2b*(1 - yhat2b)*(1 - xhat2b); 

	}
      }
    }
  }
}

void kernel(double r[], double del, double ny[], double fcnVal[]){
  double rr, drr;
  rr = r[0]*r[0] + r[1]*r[1];
  drr = r[0]*ny[0] + r[1]*ny[1];
/* 
  printf("ny[0] = %f ny[1] = %f\n",ny[0],ny[1]);
  printf("r[0] = %f r[1] = %f\n",r[0],r[1]);
  printf("drr = %f \n",drr);
*/
  fcnVal[0] = exp(-rr/(4*del))/(del*4*M_PI);
  fcnVal[1] = drr*fcnVal[0]/(2*del);
}
