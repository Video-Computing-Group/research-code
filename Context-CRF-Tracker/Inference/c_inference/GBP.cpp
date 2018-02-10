#include "GBP.h"
#include <math.h>
#include <iostream>
#include "mex.h"

GBP::~GBP() {
  gbp_arcsOrder[FATHER].clear();
  gbp_arcsOrder[SON].clear();
  delete[] gbp_arcsOrder;
  gbp_arcsOrder = 0;
  freeMessages();
}

void GBP::initMessages(double*** initMsg) {  

  defineUpdatingOrder();
  freeMessages();

  // init the messages matrix
  if (initMsg != 0) {
    gbp_messages = initMsg;
  }
  else {
    gbp_messages = new double**[ia_mrf->N];
    for (int i=0; i<ia_mrf->N; i++) {
      gbp_messages[i] = new double*[ia_mrf->neighbNum(i)];
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	int j = ia_mrf->adjMat[i][n];
	int numStates = ia_mrf->V[max(i,j)]; // should be the number of son's state
	gbp_messages[i][n] = new double[numStates];
	for (int xs=0; xs<numStates; xs++) {
	  gbp_messages[i][n][xs] = 1.0 / numStates;
	}
      }
    }
  }
}

void GBP::initBeliefs() {
  // init beliefs
  for (int i=0; i<ia_mrf->N; i++) {

    for (int xi=0; xi<ia_mrf->V[i]; xi++) {
      // start with local potential
      ia_beliefs[i][xi] = ia_mrf->localMat[i][xi];
      for (int nj=0; nj<ia_mrf->neighbNum(i); nj++) {
	int j = ia_mrf->adjMat[i][nj];
	int ni = 0;
	while (ia_mrf->adjMat[j][ni] != i) {
	  ni++;
	}
	if (j<i) { // collect message from a father
	  ia_beliefs[i][xi] *= gbp_messages[j][ni][xi];
	}
	else { // collect message from a son
	  int xj = gbp_assignInd[i][nj][xi];
	  ia_beliefs[i][xi] *= gbp_messages[j][ni][xj];
	}
      }
    }

    normalize(ia_beliefs[i],ia_mrf->V[i],0.0);
  
  }
}

void GBP::defineUpdatingOrder() {

  gbp_numArcs = 0;

  gbp_arcsOrder = new vector<int>[2];
  gbp_arcsOrder[FATHER].clear();
  gbp_arcsOrder[SON].clear();

  bool* node_status = new bool[ia_mrf->N];
  for (int i=0; i<ia_mrf->N; i++) {
    node_status[i] = false;
  }

  vector<int>* stack = new vector<int>[2];
  stack[0].push_back(0);
  stack[1].push_back(-1);
  while (!stack[0].empty()) {
    int v = stack[0].back();
    int prev = stack[1].back();
    stack[0].pop_back();
    stack[1].pop_back();
    if (prev != -1) {
      int n=0;
      if (v<prev) {
	while (ia_mrf->adjMat[v][n] != prev) {
	  n++;
	}
	gbp_arcsOrder[FATHER].push_back(v);
      }
      else {
	while (ia_mrf->adjMat[prev][n] != v) {
	  n++;
	}
	gbp_arcsOrder[FATHER].push_back(prev);
      }
      gbp_arcsOrder[SON].push_back(n);
      gbp_numArcs++;
    }
    if (!node_status[v]) {
      for (int nu=0; nu<ia_mrf->neighbNum(v); nu++) {
	int u = ia_mrf->adjMat[v][nu];
	if (!node_status[u]) {
	  stack[0].push_back(u);
	  stack[1].push_back(v);
	}
      }
      node_status[v] = true;
    }
  }

  stack[0].clear();
  stack[1].clear();
  delete[] stack;
  stack = 0;
  delete[] node_status;
  node_status = 0;
  
}


void GBP::freeMessages() {
  if (gbp_messages != 0) {
    // free the messages matrix
    for (int i=0; i<ia_mrf->N; i++) {
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	delete[] gbp_messages[i][n];
      }
      delete[] gbp_messages[i];
    }
    delete[] gbp_messages;
    gbp_messages = 0;
  }
}

void GBP::calcIncomingMessages(double* incoming_i, int i, int j) { // incoming to i without the messages from j
  
   for (int xi=0; xi<ia_mrf->V[i]; xi++) {
     incoming_i[xi] = ia_mrf->localMat[i][xi];
   }

  for (int nk=0; nk<ia_mrf->neighbNum(i); nk++) {
    int k = ia_mrf->adjMat[i][nk];
    if (k!=j) {
      int ni = 0;
      while (ia_mrf->adjMat[k][ni] != i) {
	ni++;
      }
      double sum_incoming = 0.0;
      if (k<i) { // incoming messages from a father
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  incoming_i[xi] *= gbp_messages[k][ni][xi];
	  sum_incoming += incoming_i[xi];
	}
      }
      else { // incoming messages from a son
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  // need to take the index of the son's state
	  // which complies with the father's state
	  int xk = gbp_assignInd[i][nk][xi];
	  incoming_i[xi] *= gbp_messages[k][ni][xk];
	  sum_incoming += incoming_i[xi];
	}
      }
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	incoming_i[xi] /= sum_incoming;
      }
    }
  }
}

void GBP::normalize(double* dataVec, int Vj, double epsilon) {
  double sum_data = 0.0;
  for (int xj=0; xj<Vj; xj++) {
    sum_data += dataVec[xj];
  }
  for (int xj=0; xj<Vj; xj++) {
    dataVec[xj] /= sum_data;
    if (dataVec[xj] < epsilon)
      dataVec[xj] = epsilon;
  }
}

bool GBP::isMaxMarg(double epsilon) const {
  for (int i=0; i<ia_mrf->N; i++) {
    vector<int> max_i_states;
    getMaximalStates(max_i_states, i, epsilon);
    for (int n=0; n<ia_mrf->neighbNum(i); n++) {
      int j = ia_mrf->adjMat[i][n];
      if (i<j) {
	vector<int> max_j_states;
	getMaximalStates(max_j_states, j, epsilon);
	if (max_i_states.size() != max_j_states.size()) {
	  return false;
	}
	for (int xi=0; xi<max_i_states.size(); xi++) {
	  if (gbp_assignInd[i][n][max_i_states[xi]] != max_j_states[xi]) {
	    return false;
	  }
	}
      }
    }
  }
  return true;
}

bool GBP::isSumMarg(double epsilon) const {
  for (int i=0; i<ia_mrf->N; i++) {
    for (int n=0; n<ia_mrf->neighbNum(i); n++) {
      int j = ia_mrf->adjMat[i][n];
      if (i<j) {
	double* sum_i_states = new double[ia_mrf->V[j]];
	for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	  sum_i_states[xj] = 0.0;
	}
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  int xj = gbp_assignInd[i][n][xi];
	  sum_i_states[xj] += ia_beliefs[i][xi];
	}
	for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	  if (fabs(sum_i_states[xj] - ia_beliefs[j][xj]) > epsilon) {
	    return false;
	  }
	}
      }
    }
  }
  return true;
}

double** GBP::inference(int* converged) {

  double dBel = gbp_th+1.0;
  int nIter = 0;

  while (dBel>gbp_th && nIter<gbp_maxIter) {
    nIter++;
    for (int arc=0; arc<gbp_numArcs; arc++) {
      int i = gbp_arcsOrder[FATHER][arc];
      int nj = gbp_arcsOrder[SON][arc];
      int j = ia_mrf->adjMat[i][nj];
      int ni = 0;
      while (ia_mrf->adjMat[j][ni] != i) {
	ni++;
      }

      double* incoming_i = new double[ia_mrf->V[i]];
      calcIncomingMessages(incoming_i,i,j);

      double* incoming_j = new double[ia_mrf->V[j]];
      calcIncomingMessages(incoming_j,j,i);

      // update outgoing messages
      double* outgoing_j_to_i = incoming_j; // from son j to father i
      double* outgoing_i_to_j = new double[ia_mrf->V[j]]; // from father i to son j
      for (int xj=0; xj<ia_mrf->V[j]; xj++) {

	switch (gbp_sumOrMax) {
	  case SUM:
	    outgoing_i_to_j[xj] = 0.0;
	    break;
	  case MAX:
	    outgoing_i_to_j[xj] = -1.0;
	    break;
	  default:
	    break;
	}
	
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  if (xj == gbp_assignInd[i][nj][xi]) {
	    switch (gbp_sumOrMax) {
	      case SUM:
		outgoing_i_to_j[xj] += incoming_i[xi];
		break;
	      case MAX:
		if (incoming_i[xi] > outgoing_i_to_j[xj]) {
		  outgoing_i_to_j[xj] = incoming_i[xi];
		}
		break;
	      default:
		break;
	    }
	  }
	}
      }

      normalize(outgoing_i_to_j, ia_mrf->V[j], pow(10.,-64));
      normalize(outgoing_j_to_i, ia_mrf->V[j], pow(10.,-64));
      
      // handle double-counting of regions
      
      if (gbp_bethe[j] != 1.0) {
	for (int xj=0; xj<ia_mrf->V[j]; xj++) {

	  double epsilon = pow(10.,-5);
	  if (outgoing_j_to_i[xj]<epsilon)
	    outgoing_j_to_i[xj] = epsilon;

	  if (outgoing_i_to_j[xj]<epsilon)
	    outgoing_i_to_j[xj] = epsilon;

	  double mNij_xj = (pow(outgoing_i_to_j[xj],gbp_bethe[j]) *
			    pow(outgoing_j_to_i[xj],gbp_bethe[j]-1));
	  double mNji_xj = (pow(outgoing_i_to_j[xj],gbp_bethe[j]-1) *
			    pow(outgoing_j_to_i[xj],gbp_bethe[j]));

	  outgoing_i_to_j[xj] = mNij_xj;
	  outgoing_j_to_i[xj] = mNji_xj;
	}
	
	normalize(outgoing_i_to_j, ia_mrf->V[j], pow(10.,-64));
	normalize(outgoing_j_to_i, ia_mrf->V[j], pow(10.,-64));
      }

      // update the messages
      for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	double oldi2j = gbp_messages[i][nj][xj];
	gbp_messages[i][nj][xj] = ((1.0 - gbp_alpha) * gbp_messages[i][nj][xj] +
				   gbp_alpha * outgoing_i_to_j[xj]);
	if (!(gbp_messages[i][nj][xj] > 0.0)) {
	  gbp_messages[i][nj][xj] = oldi2j;
	  nIter = gbp_maxIter + 1;
	  break;
	}

	double oldj2i = gbp_messages[j][ni][xj];
	gbp_messages[j][ni][xj] = ((1.0 - gbp_alpha) * gbp_messages[j][ni][xj] +
				   gbp_alpha * outgoing_j_to_i[xj]);
	if (!(gbp_messages[j][ni][xj] > 0.0)) {
	  gbp_messages[j][ni][xj] = oldj2i;
	  nIter = gbp_maxIter + 1;
	  break;
	}
      }
      
      delete[] incoming_i;
      incoming_i = 0;
      delete[] incoming_j;
      incoming_j = 0;
      delete[] outgoing_i_to_j;
      outgoing_i_to_j = 0;
      outgoing_j_to_i = 0;

      if (nIter>gbp_maxIter) {
	break;
      }
    }

    if (nIter>gbp_maxIter) {
      break;
    }
    
    // calculate the beliefs for regions' assignments and check for convergence

    dBel = 0.0;
    
    double** new_beliefs = new double*[ia_mrf->N];

    for (int i=0; i<ia_mrf->N; i++) {

      new_beliefs[i] = new double[ia_mrf->V[i]];
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	// start with local potential
	new_beliefs[i][xi] = ia_mrf->localMat[i][xi];
	for (int nj=0; nj<ia_mrf->neighbNum(i); nj++) {
	  int j = ia_mrf->adjMat[i][nj];
	  int ni = 0;
	  while (ia_mrf->adjMat[j][ni] != i) {
	    ni++;
	  }
	  if (j<i) { // collect message from a father
	    new_beliefs[i][xi] *= gbp_messages[j][ni][xi];
	  }
	  else { // collect message from a son
	    int xj = gbp_assignInd[i][nj][xi];
	    new_beliefs[i][xi] *= gbp_messages[j][ni][xj];
	  }
	}
      }

      normalize(new_beliefs[i],ia_mrf->V[i],0.0);

      // check convergence
      double norm_dBel_i = 0.0;
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	norm_dBel_i += ( (new_beliefs[i][xi] - ia_beliefs[i][xi]) *
			 (new_beliefs[i][xi] - ia_beliefs[i][xi]) );
      }
      norm_dBel_i = pow(norm_dBel_i, 0.5);
      
      dBel += norm_dBel_i;
    }

    // replace beliefs
    freeBeliefs();
    ia_beliefs = new_beliefs;
    new_beliefs = 0;
  }

  if (nIter > gbp_maxIter) {
    (*converged) = -1;
    mexPrintf("c-GBP: messages decreased to zero, iterating stopped\n");
  }
  else {
    if (dBel<=gbp_th) {
      (*converged) = nIter;
      mexPrintf("c-GBP: converged in %d iterations\n",nIter);
    }
    else {
      (*converged) = -1;
      mexPrintf("c-GBP: did not converge after %d iterations\n",nIter);
    }
  }
  return ia_beliefs;

}
