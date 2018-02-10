#include "LogMeanField.h"
#include "MathFunctions.h"
#include <math.h>
#include <limits>

double** LogMeanField::inference(int* converged) {
  // define a threshold for convergence
  double dBel = mf_th + 1.0;

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
      double min_qi = HUGE_VAL;
      // define ln-qi(xi) = Ei(xi,i) + sum_j_Ni[sum_xj[ exp(ln-qj(xj)) * Eij(xi,xj) ]]
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	double sum_j_xj = 0.0;
	for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	  int j = ia_mrf->adjMat[i][n];
	  for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	    sum_j_xj += (exp(-ia_beliefs[j][xj]/ia_mrf->getTemperature()) * ia_mrf->pairPotential(i,n,xi,xj));
	  }
	}

	ia_beliefs[i][xi] = ia_mrf->localMat[i][xi] + sum_j_xj;
	if (ia_beliefs[i][xi] < min_qi) {
	  min_qi = ia_beliefs[i][xi];
	}
      }
      // normalize sum_xi(ln-qi(xi)) = 0
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	ia_beliefs[i][xi] -= min_qi;
      }
    }
    
    // check convergence + assign current beliefs into previous beliefs
    dBel = -HUGE_VAL;
    for (int i=0; i<ia_mrf->N; i++) {
      double norm_dBel_i = -HUGE_VAL;
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	norm_dBel_i = AddLog(norm_dBel_i, 2*AbsSubLog(-ia_beliefs[i][xi], -prev_beliefs[i][xi])); // + log |e^thisVal - e^otherVal|^2
	prev_beliefs[i][xi] = ia_beliefs[i][xi];
      }
      norm_dBel_i = 0.5*norm_dBel_i;

      dBel = AddLog(dBel, norm_dBel_i);
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
  
  if (!lmf_logBels) {
    for (int i=0; i<ia_mrf->N; i++) {
      double sum_beliefs_i = 0.0;
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	ia_beliefs[i][xi] = exp(- ia_beliefs[i][xi] / ia_mrf->getTemperature());
	sum_beliefs_i += ia_beliefs[i][xi];
      }
      if (sum_beliefs_i > 0.0) {
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  ia_beliefs[i][xi] /= sum_beliefs_i;
	}
      }
    }
  }
  return ia_beliefs;
}
