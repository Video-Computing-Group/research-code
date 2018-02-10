#include "PairsGBP.h"
#include <math.h>
#include <iostream>
#include "mex.h"

using namespace std;
PairsGBP::~PairsGBP() {
  freeMessages();
  freePairBeliefs();
  freeBeta();
}

void PairsGBP::initMessages(double*** initMsg) {

  freeMessages();
  
  // init the messages matrix
  if (initMsg != 0) {
    pgbp_messages = initMsg;
  }
  else {
    // init the messages matrix
    pgbp_messages = new double**[ia_mrf->N];
    for (int i=0; i<ia_mrf->N; i++) {
      pgbp_messages[i] = new double*[ia_mrf->neighbNum(i)];
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	int j = ia_mrf->adjMat[i][n];
	pgbp_messages[i][n] = new double[ia_mrf->V[j]];
	for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	  pgbp_messages[i][n][xj] = 1.0 / ia_mrf->V[j];
	}
      }
    }

  }
}

void PairsGBP::freeMessages() {
  if (pgbp_messages != 0) {
    // free the messages matrix
    for (int i=0; i<ia_mrf->N; i++) {
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	delete[] pgbp_messages[i][n];
      }
      delete[] pgbp_messages[i];
    }
    delete[] pgbp_messages;
    pgbp_messages = 0;
  }
}

void PairsGBP::initPairBeliefs() {

  freePairBeliefs();
  // init the pair beliefs defultive to p(Xi=xi, Xj=xj) = Psi(xi,i)*Psi(xj,j)*Psi(xi,i,xj,j)
  pgbp_pairBeliefs = new double***[ia_mrf->N];
  for (int i=0; i<ia_mrf->N; i++) {
    pgbp_pairBeliefs[i] = new double**[ia_mrf->neighbNum(i)];
    for (int n=0; n<ia_mrf->neighbNum(i); n++) {
      pgbp_pairBeliefs[i][n] = 0;
      int j = ia_mrf->adjMat[i][n];
      if (i<j) {
	pgbp_pairBeliefs[i][n] = new double*[ia_mrf->V[i]];
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  pgbp_pairBeliefs[i][n][xi] = new double[ia_mrf->V[j]];
	  for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	    pgbp_pairBeliefs[i][n][xi][xj] = ia_mrf->pairPotential(i,n,xi,xj);
	  }
	}
      }
    }
  }
}

void PairsGBP::freePairBeliefs() {
  if (pgbp_pairBeliefs != 0) {
    for (int i=0; i<ia_mrf->N; i++) {
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	if (pgbp_pairBeliefs[i][n] != 0) {
	  for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	    delete[] pgbp_pairBeliefs[i][n][xi];
	  }
	  delete[] pgbp_pairBeliefs[i][n];
	}
      }
      delete[] pgbp_pairBeliefs[i];
    }
    delete[] pgbp_pairBeliefs;
    pgbp_pairBeliefs = 0;
  }
}

void PairsGBP::initBeta(double* doubleCount) {
  freeBeta();
  pgbp_beta = new double[ia_mrf->N];
  if (doubleCount != 0) {
    for (int i=0; i<ia_mrf->N; i++) {
      pgbp_beta[i] = 1.0 - ia_mrf->neighbNum(i) - doubleCount[i];
    }
  }
  else {
    for (int i=0; i<ia_mrf->N; i++) {
      pgbp_beta[i] = 0.0;
    }
  }
}

void PairsGBP::freeBeta() {
  if (pgbp_beta != 0) {
    delete[] pgbp_beta;
    pgbp_beta = 0;
  }
}

double**** PairsGBP::calcPairBeliefs() {

  double**** new_pairBeliefs = new double***[ia_mrf->N];

  for (int i=0; i<ia_mrf->N; i++) {
    new_pairBeliefs[i] = new double**[ia_mrf->neighbNum(i)];
    for (int n=0; n<ia_mrf->neighbNum(i); n++) {
      new_pairBeliefs[i][n] = 0;
      int j = ia_mrf->adjMat[i][n];
      if (i<j) {
	double sum_beliefs_ij = 0.0;
	new_pairBeliefs[i][n] = new double*[ia_mrf->V[i]];
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  new_pairBeliefs[i][n][xi] = new double[ia_mrf->V[j]];
	  for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	    new_pairBeliefs[i][n][xi][xj] = ia_mrf->pairPotential(i,n,xi,xj);
	    for (int ni=0; ni<ia_mrf->neighbNum(i); ni++) {
	      int k = ia_mrf->adjMat[i][ni];
	      int nk = 0;
	      while (ia_mrf->adjMat[k][nk] != i) {
		nk++;
	      }
	      if (k!=j) {
		new_pairBeliefs[i][n][xi][xj] *= pgbp_messages[k][nk][xi];
	      }
	      else {
		new_pairBeliefs[i][n][xi][xj] *= pow(pgbp_messages[k][nk][xi],(double)(pgbp_beta[i]));
	      }
	    }
	    for (int nj=0; nj<ia_mrf->neighbNum(j); nj++) {
	      int k = ia_mrf->adjMat[j][nj];
	      if (k!=i) {
		int nk = 0;
		while (ia_mrf->adjMat[k][nk] != j) {
		  nk++;
		}
		new_pairBeliefs[i][n][xi][xj] *= pgbp_messages[k][nk][xj];
	      }
	      else {
		new_pairBeliefs[i][n][xi][xj] *= pow(pgbp_messages[i][n][xj],(double)(pgbp_beta[j]));
	      }
	    }

	    sum_beliefs_ij += new_pairBeliefs[i][n][xi][xj];
	  }
	}

	// normalize the ij-beliefs
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
  freePairBeliefs();
  pgbp_pairBeliefs = new_pairBeliefs;
  new_pairBeliefs = 0;

  return pgbp_pairBeliefs;
}


double** PairsGBP::inference(int* converged) {

  double epsilon = pow(10.,-200);

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
//       double* factor = new double[ia_mrf->neighbNum(i)];
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	incoming[xi] = 1.0;
      }
      // get incoming messages
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	int j = ia_mrf->adjMat[i][n];
	int nj = 0;
	while (ia_mrf->adjMat[j][nj] != i) {
	  nj++;
	}
// 	factor[n] = 0.0;
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  incoming[xi] *= pgbp_messages[j][nj][xi];
// 	  factor[n] += incoming[xi];
	}
// 	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
// 	  incoming[xi] /= factor[n];
// 	}
      }
      
      // calculate outgoing messages
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	int j = ia_mrf->adjMat[i][n];
	int nj = 0;
	while (ia_mrf->adjMat[j][nj] != i) {
	  nj++;
	}

	double sum_outgoing_to_j = 0.0;
	double* outgoing = new double[ia_mrf->V[j]];
	
	for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	  
	  switch (pgbp_sumOrMax) {
	    case SUM:
	      outgoing[xj] = 0.0;
	      break;
	    case MAX:
	      outgoing[xj] = -1.0;
	      break;
	    default:
	      break;
	  }
	    
	  for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	    double outM = ia_mrf->pairPotential(i,n,xi,xj) * 
	      (incoming[xi] / pgbp_messages[j][nj][xi]) * pow(pgbp_messages[j][nj][xi], pgbp_beta[i]);
// 	    double outM = ia_mrf->pairPotential(i,n,xi,xj) * ia_mrf->localMat[i][xi] *
// 	      incoming[xi] * factor[n] / pgbp_messages[j][nj][xi];

	    switch (pgbp_sumOrMax) {
	      case SUM:
		outgoing[xj] += outM;
		break;
	      case MAX:
		if (outM > outgoing[xj]) {
		  outgoing[xj] = outM;
		}
		break;
	      default:
		break;
	    }
	  }
// 	  outgoing[xj] = pow(outgoing[xj], 1.0 - pgbp_beta[i]);
 	  outgoing[xj] = pow(outgoing[xj], 1.0/(1.0 - pgbp_beta[j]));

	  if (outgoing[xj] < epsilon)
	    outgoing[xj] = epsilon;
	  sum_outgoing_to_j += outgoing[xj];
	}
	for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	  if (sum_outgoing_to_j > 0.0) {
	    outgoing[xj] /= sum_outgoing_to_j;
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
//       delete[] factor;
//       factor = 0;
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
    
    dBel = 0.0;
    
    double** new_beliefs = new double*[ia_mrf->N];
    
    for (int i=0; i<ia_mrf->N; i++) {
      new_beliefs[i] = new double[ia_mrf->V[i]];
      double sum_beliefs_i = 0.0;
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	new_beliefs[i][xi] = 1.0;
	for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	  int j = ia_mrf->adjMat[i][n];
	  int nj = 0;
	  while (ia_mrf->adjMat[j][nj] != i) {
	    nj++;
	  }
	  new_beliefs[i][xi] *= pgbp_messages[j][nj][xi];
	}
	sum_beliefs_i += new_beliefs[i][xi];
      }

      double norm_dBel_i = 0.0;
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	if (sum_beliefs_i > 0.0) {
	  new_beliefs[i][xi] /= sum_beliefs_i;
	}
	norm_dBel_i += pow((new_beliefs[i][xi] - ia_beliefs[i][xi]), 2.0);
      }
      norm_dBel_i = pow(norm_dBel_i, 0.5);

      dBel += norm_dBel_i;
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
    mexPrintf("c-PairsGBP: converged in %d iterations\n",nIter);
    //    cout << "c-GBP: converged in " << nIter << " iterations " << endl;
  }
  else {
    (*converged) = -1;
    mexPrintf("c-PairsGBP: did not converge\n");
    //    cout << "c-GBP: did not converge" << endl;
  }
  
  return ia_beliefs;
  
}
