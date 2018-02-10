#include "MeanField.h"
#include <math.h>

double** MeanField::inference(int* converged) {
  // define a threshold for convergence
  double dBel = 0.0;
  double epsilon = pow(10.,-200);

  double** prev_beliefs = new double*[ia_mrf->N];
  for (int i=0; i<ia_mrf->N; i++) {
    prev_beliefs[i] = new double[ia_mrf->V[i]];
    for (int xi=0; xi<ia_mrf->V[i]; xi++) {
      prev_beliefs[i][xi] = ia_beliefs[i][xi];
    }
  }
  int t;
  for (t=1; t<=mf_maxIter; t++) {
    for (int i=0; i<ia_mrf->N; i++) {
      double sum_qi = 0.0;
      // define qi(xi) = Psi(xi,i)*exp{ sum_j_Ni[sum_xj[ qj(xj) * ln(Psi(xi,xj)) ]] }
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	double sum_j_xj = 0.0;
	for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	  int j = ia_mrf->adjMat[i][n];
	  for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	    sum_j_xj += (ia_beliefs[j][xj] * log(ia_mrf->pairPotential(i,n,xi,xj) + epsilon));
	  }
	}

	ia_beliefs[i][xi] = ia_mrf->localMat[i][xi] * exp(sum_j_xj);
	sum_qi += ia_beliefs[i][xi];
      }
      // normalize sum_xi(qi(xi)) = 1
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	ia_beliefs[i][xi] /= sum_qi;
      }
    }
    
    // check convergence + assign current beliefs into previous beliefs
    dBel = 0.0;
    for (int i=0; i<ia_mrf->N; i++) {

      double norm_dBel_i = 0.0;
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	norm_dBel_i += pow((ia_beliefs[i][xi] - prev_beliefs[i][xi]), 2.0);
	prev_beliefs[i][xi] = ia_beliefs[i][xi];
      }
      norm_dBel_i = pow(norm_dBel_i, 0.5);

      dBel += norm_dBel_i;
    }
    if (dBel<mf_th) {
      break;
    }
  }
  // free memory
  for (int i=0; i<ia_mrf->N; i++) {
    delete[] prev_beliefs[i];
  }
  delete[] prev_beliefs;
  prev_beliefs = 0;

  if (dBel<=mf_th) {
    (*converged) = t;
  }
  else {
    (*converged) = -1;
  }
  
  return ia_beliefs;
}
