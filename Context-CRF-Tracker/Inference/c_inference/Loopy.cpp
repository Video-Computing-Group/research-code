#include "Loopy.h"
#include <math.h>
#include <iostream>
#include "mex.h"

using namespace std;
Loopy::~Loopy() {
  freeMessages();
  freePairBeliefs();
}

void Loopy::initMessages(double*** initMsg) {

  freeMessages();
  
  // init the messages matrix
  if (initMsg != 0) {
    l_messages = initMsg;
  }
  else {
    // init the messages matrix
    l_messages = new double**[ia_mrf->N];
    for (int i=0; i<ia_mrf->N; i++) {
      l_messages[i] = new double*[ia_mrf->neighbNum(i)];
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	int j = ia_mrf->adjMat[i][n];
	l_messages[i][n] = new double[ia_mrf->V[j]];
	for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	  l_messages[i][n][xj] = 1.0 / ia_mrf->V[j];
	}
      }
    }

  }
}

void Loopy::freeMessages() {
  if (l_messages != 0) {
    // free the messages matrix
    for (int i=0; i<ia_mrf->N; i++) {
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	delete[] l_messages[i][n];
      }
      delete[] l_messages[i];
    }
    delete[] l_messages;
    l_messages = 0;
  }
}

void Loopy::initPairBeliefs() {

  freePairBeliefs();
  // init the pair beliefs defultive to p(Xi=xi, Xj=xj) = Psi(xi,i)*Psi(xj,j)*Psi(xi,i,xj,j)
  l_pairBeliefs = new double***[ia_mrf->N];
  for (int i=0; i<ia_mrf->N; i++) {
    l_pairBeliefs[i] = new double**[ia_mrf->neighbNum(i)];
    for (int n=0; n<ia_mrf->neighbNum(i); n++) {
      l_pairBeliefs[i][n] = 0;
      int j = ia_mrf->adjMat[i][n];
      if (i<j) {
	l_pairBeliefs[i][n] = new double*[ia_mrf->V[i]];
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  l_pairBeliefs[i][n][xi] = new double[ia_mrf->V[j]];
	  for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	    l_pairBeliefs[i][n][xi][xj] = (ia_mrf->localMat[i][xi] *
					   ia_mrf->localMat[j][xj] *
					   ia_mrf->pairPotential(i,n,xi,xj));
	  }
	}
      }
    }
  }
}

void Loopy::freePairBeliefs() {
  if (l_pairBeliefs != 0) {
    for (int i=0; i<ia_mrf->N; i++) {
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	if (l_pairBeliefs[i][n] != 0) {
	  for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	    delete[] l_pairBeliefs[i][n][xi];
	  }
	  delete[] l_pairBeliefs[i][n];
	}
      }
      delete[] l_pairBeliefs[i];
    }
    delete[] l_pairBeliefs;
    l_pairBeliefs = 0;
  }
}

void Loopy::initTRWRho(double** trwRho) {
  freeTRWRho();
  if (trwRho != 0) {
    l_trwRho = new double*[ia_mrf->N];
    for (int i=0; i<ia_mrf->N; i++) {
      l_trwRho[i] = new double[ia_mrf->neighbNum(i)];
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	l_trwRho[i][n] = trwRho[i][n];
      }
    }
  }
}

void Loopy::freeTRWRho() {
  if (l_trwRho != 0) {
    for (int i=0; i<ia_mrf->N; i++) {
      delete[] l_trwRho[i];
      l_trwRho[i] = 0;
    }
    delete[] l_trwRho;
    l_trwRho = 0;
  }
}

double**** Loopy::calcPairBeliefs() {

  if (l_trwRho != 0) {
    return calcPairBeliefsTRBP();
  }
  
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
	    new_pairBeliefs[i][n][xi][xj] = (ia_mrf->localMat[i][xi] *
					     ia_mrf->localMat[j][xj] *
					     ia_mrf->pairPotential(i,n,xi,xj));
	    for (int ni=0; ni<ia_mrf->neighbNum(i); ni++) {
	      int k = ia_mrf->adjMat[i][ni];
	      if (k!=j) {
		int nk = 0;
		while (ia_mrf->adjMat[k][nk] != i) {
		  nk++;
		}
		new_pairBeliefs[i][n][xi][xj] *= l_messages[k][nk][xi];
	      }
	    }
	    for (int nj=0; nj<ia_mrf->neighbNum(j); nj++) {
	      int k = ia_mrf->adjMat[j][nj];
	      if (k!=i) {
		int nk = 0;
		while (ia_mrf->adjMat[k][nk] != j) {
		  nk++;
		}
		new_pairBeliefs[i][n][xi][xj] *= l_messages[k][nk][xj];
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
  l_pairBeliefs = new_pairBeliefs;
  new_pairBeliefs = 0;

  return l_pairBeliefs;
}

double**** Loopy::calcPairBeliefsTRBP() {

  double**** new_pairBeliefs = new double***[ia_mrf->N];

  for (int i=0; i<ia_mrf->N; i++) {
    new_pairBeliefs[i] = new double**[ia_mrf->neighbNum(i)];
    for (int n=0; n<ia_mrf->neighbNum(i); n++) {
      new_pairBeliefs[i][n] = 0;
      int j = ia_mrf->adjMat[i][n];
      if (i<j) {
	int nei_ij = 0;
	while (ia_mrf->adjMat[j][nei_ij] != i) {
	  nei_ij++;
	}
	double sum_beliefs_ij = 0.0;
	new_pairBeliefs[i][n] = new double*[ia_mrf->V[i]];
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  new_pairBeliefs[i][n][xi] = new double[ia_mrf->V[j]];
	  for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	    new_pairBeliefs[i][n][xi][xj] = (ia_mrf->localMat[i][xi] *
					     ia_mrf->localMat[j][xj] *
					     ia_mrf->pairPotential(i,n,xi,xj));
// 	    new_pairBeliefs[i][n][xi][xj] = (ia_mrf->localMat[i][xi] *
// 					     ia_mrf->localMat[j][xj] *
// 					     pow(ia_mrf->pairPotential(i,n,xi,xj),(1.0 / l_trwRho[i][n])));

	    for (int ni=0; ni<ia_mrf->neighbNum(i); ni++) {
	      int k = ia_mrf->adjMat[i][ni];
	      int nk = 0;
	      while (ia_mrf->adjMat[k][nk] != i) {
		nk++;
	      }
	      new_pairBeliefs[i][n][xi][xj] *= pow(l_messages[k][nk][xi],l_trwRho[i][ni]);
	    }
	    new_pairBeliefs[i][n][xi][xj] /= l_messages[j][nei_ij][xi];

	    for (int nj=0; nj<ia_mrf->neighbNum(j); nj++) {
	      int k = ia_mrf->adjMat[j][nj];
	      int nk = 0;
	      while (ia_mrf->adjMat[k][nk] != j) {
		nk++;
	      }
	      new_pairBeliefs[i][n][xi][xj] *= pow(l_messages[k][nk][xj],l_trwRho[j][nj]);
	    }
	    new_pairBeliefs[i][n][xi][xj] /= l_messages[i][n][xj];

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
  l_pairBeliefs = new_pairBeliefs;
  new_pairBeliefs = 0;

  return l_pairBeliefs;
}

double** Loopy::inference(int* converged) {

  if (l_trwRho != 0) {
    return inferenceTRBP(converged);
  }
  
  double epsilon = pow(10.,-200);

  double dBel = l_th+1.0;
  int nIter = 0;
  
  double*** new_messages = 0;
  if (l_strategy == PARALLEL) {
    new_messages = new double**[ia_mrf->N];
    for (int i=0; i<ia_mrf->N; i++) {
      new_messages[i] = new double*[ia_mrf->neighbNum(i)];
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	int j = ia_mrf->adjMat[i][n];
	new_messages[i][n] = new double[ia_mrf->V[j]];
	for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	  new_messages[i][n][xj] = l_messages[i][n][xj];
	}
      }
    }
  }
  while (dBel>l_th && nIter<l_maxIter) {
    nIter++;

    for (int i=0; i<ia_mrf->N; i++) {

      // init the incoming messages to 1
      double* incoming = new double[ia_mrf->V[i]];
      double* factor = new double[ia_mrf->neighbNum(i)];
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
	factor[n] = 0.0;
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  incoming[xi] *= l_messages[j][nj][xi];
	  factor[n] += incoming[xi];
	}
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  incoming[xi] /= factor[n];
	}
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
	  
	  switch (l_sumOrMax) {
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
	    double outM = ia_mrf->pairPotential(i,n,xi,xj) * ia_mrf->localMat[i][xi] *
	      incoming[xi] / l_messages[j][nj][xi];
// 	    double outM = ia_mrf->pairPotential(i,n,xi,xj) * ia_mrf->localMat[i][xi] *
// 	      incoming[xi] * factor[n] / l_messages[j][nj][xi];

	    switch (l_sumOrMax) {
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
	    
	  if (outgoing[xj] < epsilon)
	    outgoing[xj] = epsilon;
	  sum_outgoing_to_j += outgoing[xj];
	}
	for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	  if (sum_outgoing_to_j > 0.0) {
	    outgoing[xj] /= sum_outgoing_to_j;
	  }
	    
	  if (!(outgoing[xj]>0.0)) {
 	    nIter = l_maxIter + 1;
 	    break;
	  }
	  
	  switch (l_strategy) {

	    case SEQUENTIAL:
	      l_messages[i][n][xj] = outgoing[xj];
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
      delete[] factor;
      factor = 0;
    }

    if (l_strategy == PARALLEL) {
      for (int i=0; i<ia_mrf->N; i++) {
	for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	  int j = ia_mrf->adjMat[i][n];
	  for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	    l_messages[i][n][xj] = new_messages[i][n][xj];
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
	new_beliefs[i][xi] = ia_mrf->localMat[i][xi];
	for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	  int j = ia_mrf->adjMat[i][n];
	  int nj = 0;
	  while (ia_mrf->adjMat[j][nj] != i) {
	    nj++;
	  }
	  new_beliefs[i][xi] *= l_messages[j][nj][xi];
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

  if (l_strategy == PARALLEL) {
    for (int i=0; i<ia_mrf->N; i++) {
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	delete[] new_messages[i][n];
      }
      delete[] new_messages[i];
    }
    delete[] new_messages;
    new_messages = 0;
  }

  if (nIter > l_maxIter) {
    (*converged) = -1;
    mexPrintf("c-Loopy: messages decreased to zero, iterating stopped\n");
  }
  else {
    if (dBel<=l_th) {
      (*converged) = nIter;
      //mexPrintf("c-Loopy: converged in %d iterations\n",nIter);
    }
    else {
      (*converged) = -1;
      //mexPrintf("c-Loopy: did not converge after %d iterations\n",nIter);
    }
  }
  
  return ia_beliefs;
  
}

double** Loopy::inferenceTRBP(int* converged) {

  double epsilon = pow(10.,-200);

  double dBel = l_th+1.0;
  int nIter = 0;
  
  double*** new_messages = 0;
  if (l_strategy == PARALLEL) {
    new_messages = new double**[ia_mrf->N];
    for (int i=0; i<ia_mrf->N; i++) {
      new_messages[i] = new double*[ia_mrf->neighbNum(i)];
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	int j = ia_mrf->adjMat[i][n];
	new_messages[i][n] = new double[ia_mrf->V[j]];
	for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	  new_messages[i][n][xj] = l_messages[i][n][xj];
	}
      }
    }
  }

  while (dBel>l_th && nIter<l_maxIter) {
    nIter++;

    for (int i=0; i<ia_mrf->N; i++) {

      // init the incoming messages to 1
      double* incoming = new double[ia_mrf->V[i]];
      double* factor = new double[ia_mrf->neighbNum(i)];
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
	factor[n] = 0.0;
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  incoming[xi] *= pow(l_messages[j][nj][xi], l_trwRho[i][n]);
	  factor[n] += incoming[xi];
	}
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  incoming[xi] /= factor[n];
	}
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
	  
	  switch (l_sumOrMax) {
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
	    double outM = ia_mrf->pairPotential(i,n,xi,xj) * ia_mrf->localMat[i][xi] *
	      incoming[xi] / l_messages[j][nj][xi]; // the pair-potentials are raised by 1/rho in the matlab interface

	    switch (l_sumOrMax) {
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

	  if (outgoing[xj] < epsilon)
	    outgoing[xj] = epsilon;
	  sum_outgoing_to_j += outgoing[xj];
	}

	for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	  if (sum_outgoing_to_j > 0.0) {
	    outgoing[xj] /= sum_outgoing_to_j;
	  }
	    
	  if (!(outgoing[xj]>0.0)) {
 	    nIter = l_maxIter + 1;
 	    break;
	  }
	    
	  switch (l_strategy) {

	    case SEQUENTIAL:
	      l_messages[i][n][xj] = outgoing[xj];
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
      delete[] factor;
      factor = 0;
    }

    if (l_strategy == PARALLEL) {
      for (int i=0; i<ia_mrf->N; i++) {
	for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	  int j = ia_mrf->adjMat[i][n];
	  for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	    l_messages[i][n][xj] = new_messages[i][n][xj];
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
	new_beliefs[i][xi] = ia_mrf->localMat[i][xi];
	for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	  int j = ia_mrf->adjMat[i][n];
	  int nj = 0;
	  while (ia_mrf->adjMat[j][nj] != i) {
	    nj++;
	  }
	  new_beliefs[i][xi] *= pow(l_messages[j][nj][xi],l_trwRho[i][n]);
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

  if (l_strategy == PARALLEL) {
    for (int i=0; i<ia_mrf->N; i++) {
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	delete[] new_messages[i][n];
      }
      delete[] new_messages[i];
    }
    delete[] new_messages;
    new_messages = 0;
  }

  if (nIter > l_maxIter) {
    (*converged) = -1;
    mexPrintf("c-Loopy: messages decreased to zero, iterating stopped\n");
  }
  else {
    if (dBel<=l_th) {
      (*converged) = nIter;
      mexPrintf("c-Loopy: converged in %d iterations\n",nIter);
    }
    else {
      (*converged) = -1;
      mexPrintf("c-Loopy: did not converge after %d iterations\n",nIter);
    }
  }
  
  return ia_beliefs;
  
}
