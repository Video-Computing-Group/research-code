#include "LogPairsGBP.h"
#include "MathFunctions.h"
#include <math.h>
#include <iostream>
#include "mex.h"

using namespace std;

double**** LogPairsGBP::calcPairBeliefs() {

  double**** new_pairBeliefs = new double***[ia_mrf->N];

  for (int i=0; i<ia_mrf->N; i++) {
    new_pairBeliefs[i] = new double**[ia_mrf->neighbNum(i)];
    for (int n=0; n<ia_mrf->neighbNum(i); n++) {
      new_pairBeliefs[i][n] = 0;
      int j = ia_mrf->adjMat[i][n];
      if (i<j) {
	double min_beliefs_ij = numeric_limits<double>::infinity();
	new_pairBeliefs[i][n] = new double*[ia_mrf->V[i]];
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  new_pairBeliefs[i][n][xi] = new double[ia_mrf->V[j]];
	  for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	    new_pairBeliefs[i][n][xi][xj] = ia_mrf->pairEnergy(i,n,xi,xj);
	    for (int ni=0; ni<ia_mrf->neighbNum(i); ni++) {
	      int k = ia_mrf->adjMat[i][ni];
	      int nk = 0;
	      while (ia_mrf->adjMat[k][nk] != i) {
		nk++;
	      }
	      if (k!=j) {
		new_pairBeliefs[i][n][xi][xj] += pgbp_messages[k][nk][xi];
	      }
	      else {
		new_pairBeliefs[i][n][xi][xj] += (pgbp_messages[k][nk][xi]*pgbp_beta[i]);
	      }
	    }
	    for (int nj=0; nj<ia_mrf->neighbNum(j); nj++) {
	      int k = ia_mrf->adjMat[j][nj];
	      if (k!=i) {
		int nk = 0;
		while (ia_mrf->adjMat[k][nk] != j) {
		  nk++;
		}
		new_pairBeliefs[i][n][xi][xj] += pgbp_messages[k][nk][xj];
	      }
	      else {
		new_pairBeliefs[i][n][xi][xj] += (pgbp_messages[i][n][xj]*pgbp_beta[j]);
	      }
	    }

	    if (new_pairBeliefs[i][n][xi][xj] < min_beliefs_ij) {
	      min_beliefs_ij = new_pairBeliefs[i][n][xi][xj];
	    }
	  }
	}

	// normalize the ij-beliefs
	if (lpgbp_logBels) {
	  for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	    for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	      new_pairBeliefs[i][n][xi][xj] -= min_beliefs_ij;
	    }
	  }
	}
	else {
	  double sum_beliefs_ij = 0.0;
	  for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	    for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	      new_pairBeliefs[i][n][xi][xj] -= min_beliefs_ij;
	      new_pairBeliefs[i][n][xi][xj] = exp(- new_pairBeliefs[i][n][xi][xj] / ia_mrf->getTemperature());
	      sum_beliefs_ij += new_pairBeliefs[i][n][xi][xj];
	    }
	  }
	  if (sum_beliefs_ij > 0.0) {
	    for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	      for (int xj=0; xj<ia_mrf->V[j]; xj++) {
		new_pairBeliefs[i][n][xi][xj] /= sum_beliefs_ij;
	      }
	    }
	  }
	}
      }
    }
  }
  freePairBeliefs();
  pgbp_pairBeliefs = new_pairBeliefs;
  new_pairBeliefs = 0;

  return pgbp_pairBeliefs;
}


double** LogPairsGBP::inference(int* converged) {
  double dBel = pgbp_th+1.0;
  int nIter = 0;
  
  double*** new_messages = 0;
  if (pgbp_strategy == PARALLEL) {
    new_messages = new double**[ia_mrf->N];
    for (int i=0; i<ia_mrf->N; i++) {
      new_messages[i] = new double*[ia_mrf->neighbNum(i)];
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	int j = ia_mrf->adjMat[i][n];
	new_messages[i][n] = new double[ia_mrf->V[j]];
	for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	  new_messages[i][n][xj] = pgbp_messages[i][n][xj];
	}
      }
    }
  }
  while (dBel>pgbp_th && nIter<pgbp_maxIter) {
    nIter++;

    for (int i=0; i<ia_mrf->N; i++) {

      // init the incoming messages to 1
      double* incoming = new double[ia_mrf->V[i]];
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	incoming[xi] = 0.0;
      }
      // get incoming messages
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	int j = ia_mrf->adjMat[i][n];
	int nj = 0;
	while (ia_mrf->adjMat[j][nj] != i) {
	  nj++;
	}
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  incoming[xi] += pgbp_messages[j][nj][xi];
	}
      }
      
      // calculate outgoing messages
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	int j = ia_mrf->adjMat[i][n];
	int nj = 0;
	while (ia_mrf->adjMat[j][nj] != i) {
	  nj++;
	}

	double norm_msg = HUGE_VAL;
	double* outgoing = new double[ia_mrf->V[j]];
	
	for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	  
	  switch (pgbp_sumOrMax) {
	    case SUM:
	      outgoing[xj] = -HUGE_VAL;
	      break;
	    case MAX:
	      outgoing[xj] = numeric_limits<double>::infinity();
	      break;
	    default:
	      break;
	  }
	    
	  for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	    double outM = ia_mrf->pairEnergy(i,n,xi,xj) + incoming[xi] +
	      (pgbp_messages[j][nj][xi] * (pgbp_beta[i] - 1));

	    switch (pgbp_sumOrMax) {
	      case SUM:
    		outgoing[xj] = AddLogFactor(outgoing[xj],outM,-ia_mrf->getTemperature());
		break;
	      case MAX:
		if (outM < outgoing[xj]) {
		  outgoing[xj] = outM;
		}
		break;
	      default:
		break;
	    }
	  }
	  outgoing[xj] *= (1 / (1 - pgbp_beta[j]));
	  if (outgoing[xj] < norm_msg) {
	    norm_msg = outgoing[xj];
	  }
	}
	for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	  if (norm_msg != -HUGE_VAL) {
	    outgoing[xj] -= norm_msg;
	  }
	    
	  switch (pgbp_strategy) {

	    case SEQUENTIAL:
	      pgbp_messages[i][n][xj] = outgoing[xj];
	      break;

	    case PARALLEL:
	      new_messages[i][n][xj] = outgoing[xj];
	      break;

	    default:
	      break;
	  }
	}
	delete[] outgoing;
	outgoing = 0;
      }
      delete[] incoming;
      incoming = 0;
    }

    if (pgbp_strategy == PARALLEL) {
      for (int i=0; i<ia_mrf->N; i++) {
	for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	  int j = ia_mrf->adjMat[i][n];
	  for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	    pgbp_messages[i][n][xj] = new_messages[i][n][xj];
	  }
	}
      }
    }

    // update beliefs and check for convergence
    
     dBel = -HUGE_VAL;
    
    double** new_beliefs = new double*[ia_mrf->N];
    
    for (int i=0; i<ia_mrf->N; i++) {
      new_beliefs[i] = new double[ia_mrf->V[i]];
      double min_beliefs_i = HUGE_VAL;
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	new_beliefs[i][xi] = 0.0;
	for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	  int j = ia_mrf->adjMat[i][n];
	  int nj = 0;
	  while (ia_mrf->adjMat[j][nj] != i) {
	    nj++;
	  }
	  new_beliefs[i][xi] += pgbp_messages[j][nj][xi];
	}
	if (new_beliefs[i][xi] < min_beliefs_i) {
	  min_beliefs_i = new_beliefs[i][xi];
	}
      }
      double norm_dBel_i = -HUGE_VAL;
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
 	new_beliefs[i][xi] -= min_beliefs_i;
	norm_dBel_i = AddLog(norm_dBel_i, 2*AbsSubLog(-new_beliefs[i][xi], -ia_beliefs[i][xi])); // + log |e^thisVal - e^otherVal|^2
      }
      norm_dBel_i = 0.5*norm_dBel_i;

      dBel = AddLog(dBel, norm_dBel_i);
    }
    freeBeliefs();
    ia_beliefs = new_beliefs;
    new_beliefs = 0;

  }

  if (pgbp_strategy == PARALLEL) {
    for (int i=0; i<ia_mrf->N; i++) {
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	delete[] new_messages[i][n];
      }
      delete[] new_messages[i];
    }
    delete[] new_messages;
    new_messages = 0;
  }

  if (dBel<=pgbp_th) {
    (*converged) = nIter;
    mexPrintf("c-LogPairsGBP: converged in %d iterations\n",nIter);
    //    cout << "c-GBP: converged in " << nIter << " iterations " << endl;
  }
  else {
    (*converged) = -1;
    mexPrintf("c-LogPairsGBP: did not converge\n");
    //    cout << "c-GBP: did not converge" << endl;
  }

  if (!lpgbp_logBels) {
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
