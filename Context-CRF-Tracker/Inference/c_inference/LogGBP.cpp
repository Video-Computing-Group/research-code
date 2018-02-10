#include "LogGBP.h"
#include "MathFunctions.h"
#include <math.h>
#include <iostream>
#include <limits>
#include "mex.h"

void LogGBP::initMessages(double*** initMsg) {  

  // init the messages matrix
  if (initMsg != 0) {
    freeMessages();
    gbp_messages = initMsg;
  }
  else {
    for (int i=0; i<ia_mrf->N; i++) {
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	int j = ia_mrf->adjMat[i][n];
	int numStates = ia_mrf->V[max(i,j)]; // should be the number of son's state
	for (int xs=0; xs<numStates; xs++) {
 	  gbp_messages[i][n][xs] = 0.0;//-ia_mrf->getTemperature()*log(1.0 / numStates);
	}
      }
    }
  }
}

void LogGBP::initBeliefs() {  

  for (int i=0; i<ia_mrf->N; i++) {
    for (int xi=0; xi<ia_mrf->V[i]; xi++) {
      ia_beliefs[i][xi] = 0.0;
    }
  }
}

void LogGBP::calcIncomingMessages(double* incoming_i, int i, int j) { // incoming to i without the messages from j

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
      if (k<i) { // incoming messages from a father
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  incoming_i[xi] += gbp_messages[k][ni][xi];
	}
      }
      else { // incoming messages from a son
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  // need to take the index of the son's state
	  // which complies with the father's state
	  int xk = gbp_assignInd[i][nk][xi];
	  incoming_i[xi] += gbp_messages[k][ni][xk];
	}
      }
    }
  }

  minNormalize(incoming_i, ia_mrf->V[i]);
}

void LogGBP::minNormalize(double* dataVec, int Vj) {
  double min_data = numeric_limits<double>::infinity();
  for (int xj=0; xj<Vj; xj++) {
    if (dataVec[xj] < min_data) {
      min_data = dataVec[xj];
    }
  }
  for (int xj=0; xj<Vj; xj++) {
    dataVec[xj] -= min_data;
  }
}

double** LogGBP::inference(int* converged) {

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
	    outgoing_i_to_j[xj] = -HUGE_VAL;
	    break;
	  case MAX:
	    outgoing_i_to_j[xj] = numeric_limits<double>::infinity();
	    break;
	  default:
	    break;
	}
	
	for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	  if (xj == gbp_assignInd[i][nj][xi]) {
	    switch (gbp_sumOrMax) {
	      case SUM:
  		outgoing_i_to_j[xj] = AddLogFactor(outgoing_i_to_j[xj],incoming_i[xi],-ia_mrf->getTemperature());
// 		outgoing_i_to_j[xj] = -ia_mrf->getTemperature()*AddLog(-outgoing_i_to_j[xj]/ia_mrf->getTemperature(),
// 								       -incoming_i[xi]/ia_mrf->getTemperature());
		break;
	      case MAX:
		if (incoming_i[xi] < outgoing_i_to_j[xj]) {
		  outgoing_i_to_j[xj] = incoming_i[xi];
		}
		break;
	      default:
		break;
	    }
	  }
	}
      }

      minNormalize(outgoing_i_to_j, ia_mrf->V[j]);
      minNormalize(outgoing_j_to_i, ia_mrf->V[j]);
      
      // handle double-counting of regions
      
      if (gbp_bethe[j] != 1.0) {
	for (int xj=0; xj<ia_mrf->V[j]; xj++) {

	  double mNij_xj = (gbp_bethe[j]*outgoing_i_to_j[xj]) + ((gbp_bethe[j]-1)*outgoing_j_to_i[xj]);
	  double mNji_xj = ((gbp_bethe[j]-1)*outgoing_i_to_j[xj]) + (gbp_bethe[j]*outgoing_j_to_i[xj]);

	  outgoing_i_to_j[xj] = mNij_xj;
	  outgoing_j_to_i[xj] = mNji_xj;
	}
	
	minNormalize(outgoing_i_to_j, ia_mrf->V[j]);
	minNormalize(outgoing_j_to_i, ia_mrf->V[j]);
      }
	
      // update the messages
      double T = ia_mrf->getTemperature();
      for (int xj=0; xj<ia_mrf->V[j]; xj++) {
	double oldi2j = gbp_messages[i][nj][xj];
	gbp_messages[i][nj][xj] = AddLogFactor(gbp_messages[i][nj][xj] - T*log(1.0 - gbp_alpha),
					       outgoing_i_to_j[xj] - T*log(gbp_alpha),
					       -T);
	if (gbp_messages[i][nj][xj] != gbp_messages[i][nj][xj]) {
	  gbp_messages[i][nj][xj] = oldi2j;
	  nIter = gbp_maxIter + 1;
	  break;
	}

	double oldj2i = gbp_messages[j][ni][xj];
	gbp_messages[j][ni][xj] = AddLogFactor(gbp_messages[j][ni][xj] - T*log(1.0 - gbp_alpha),
					       outgoing_j_to_i[xj] - T*log(gbp_alpha),
					       -T);
	if (gbp_messages[j][ni][xj] != gbp_messages[j][ni][xj]) {
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

    dBel = -HUGE_VAL;
    
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
	    new_beliefs[i][xi] += gbp_messages[j][ni][xi];
	  }
	  else { // collect message from a son
	    int xj = gbp_assignInd[i][nj][xi];
	    new_beliefs[i][xi] += gbp_messages[j][ni][xj];
	  }
	}
      }

      minNormalize(new_beliefs[i],ia_mrf->V[i]);

      // check convergence
      double norm_dBel_i = -HUGE_VAL;
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	if (ia_beliefs[i][xi] != ia_beliefs[i][xi]) {
	  cout << "iter = " << nIter << " i = " << i << " xi = " << xi << endl;
	}
    	norm_dBel_i = AddLog(norm_dBel_i, 2*AbsSubLog(-new_beliefs[i][xi], -ia_beliefs[i][xi])); // + log |e^thisVal - e^otherVal|^2
	// + log |e^thisVal - e^otherVal|^2
//   	norm_dBel_i = AddLogFactor(norm_dBel_i,2*AbsSubLogFactor(-new_beliefs[i][xi], -ia_beliefs[i][xi],
//  								 ia_mrf->getTemperature()),
//  				   -ia_mrf->getTemperature());
      }
      norm_dBel_i = 0.5*norm_dBel_i;

      dBel = AddLog(dBel, norm_dBel_i);
//       dBel = AddLogFactor(dBel, norm_dBel_i, -ia_mrf->getTemperature());

    }

    // replace beliefs
    freeBeliefs();
    ia_beliefs = new_beliefs;
    new_beliefs = 0;
  }

  if (nIter > gbp_maxIter) {
    (*converged) = -1;
    mexPrintf("c-LogGBP: messages decreased to zero, iterating stopped\n");
  }
  else {
    if (dBel<=gbp_th) {
      (*converged) = nIter;
      mexPrintf("c-LogGBP: converged in %d iterations\n",nIter);
    }
    else {
      cout << "dBel = " << dBel << " gbp_th = " << gbp_th << endl;
      (*converged) = -1;
      mexPrintf("c-LogGBP: did not converge after %d iterations\n",nIter);
    }
  }

  if (!lgbp_logBels) {
    for (int i=0; i<ia_mrf->N; i++) {
      double sum_beliefs_i = 0.0;
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	ia_beliefs[i][xi] = exp(- ia_beliefs[i][xi] / ia_mrf->getTemperature());
	sum_beliefs_i += ia_beliefs[i][xi];
      }
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	if (sum_beliefs_i > 0.0) {
	  ia_beliefs[i][xi] /= sum_beliefs_i;
	}
      }
    }
  }
  return ia_beliefs;

}
