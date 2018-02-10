#include "MonteCarlo.h"
#include <stdlib.h>
#include <iostream>

MonteCarlo::~MonteCarlo() {
  for (int i=0; i<ia_mrf->N; i++) {
    delete[] mc_Pxi[i];
  }
  delete[] mc_Pxi;
  mc_Pxi = 0;
  delete[] mc_Pi;
  mc_Pi = 0;

  freeState();
}

void MonteCarlo::freeState() {
  if (mc_currX != 0) {
    delete[] mc_currX;
    mc_currX = 0;
  }
}

void MonteCarlo::init(int* startX) {
  mc_Pi = new double[ia_mrf->N];
  mc_Pxi = new double*[ia_mrf->N];
  for (int i=0; i<ia_mrf->N; i++) {
    mc_Pi[i] = 1.0 / ia_mrf->N;
    mc_Pxi[i] = new double[ia_mrf->V[i]];
    for (int xi=0; xi<ia_mrf->V[i]; xi++) {
      mc_Pxi[i][xi] = 1.0 / ia_mrf->V[i];
    }
  }

  initState(startX);
}

void MonteCarlo::initState(int* startX) {
  freeState();
  mc_currX = new int[ia_mrf->N];
  for (int i=0; i<ia_mrf->N; i++) {
    mc_currX[i] = startX[i];
  }
}

// will return an integer i from [0,..,n-1], in the given probability:
// P(i) = q[i]
int MonteCarlo::chooseInteger(double* q, int n) {

  // choose from [0..1]
  double randomal = 1.0 * rand()/(RAND_MAX);

  // choose i as the bin in q where randomal has fallen
  double cumsum = 0.0;
  int i=0;
  for (int j=0; j<n; j++) {
    cumsum += q[j];
    if (randomal<cumsum) {
      i = j;
      break;
    }
  }
  return i;
  
}

void MonteCarlo::getCurrentState(int* currX) const {
  for (int i=0; i<ia_mrf->N; i++) {
    currX[i] = mc_currX[i];
  }
}

void MonteCarlo::reachEquilibrium() {
  for (int t=0; t<mc_burning; t++) {
    transition();
  }
}

double** MonteCarlo::inference(int* converged) {

  // count from the samples the number of times the value xi
  // appeared for node i
  for (int i=0; i<ia_mrf->N; i++) {
    for (int xi=0; xi<ia_mrf->V[i]; xi++) {
      ia_beliefs[i][xi] = 0.0;
    }
  }

  reachEquilibrium();
  int total = mc_S*(mc_interval+1);
  // transit between states untill the burning-time, then sample
  // with the defined interval
  for (int t=1; t<=total; t++) {
    bool sample = ((mc_interval+1) * (int)(t / (mc_interval+1)) == t);

    transition();
    if (sample) {
      for (int i=0; i<ia_mrf->N; i++) {
	int xi = mc_currX[i];
	ia_beliefs[i][xi] += 1.0;
      }
    }
  }
  for (int i=0; i<ia_mrf->N; i++) {
    double sum_i_xi = 0.0;
    for (int xi=0; xi<ia_mrf->V[i]; xi++) {
      ia_beliefs[i][xi] /= (double)mc_S;
      sum_i_xi += ia_beliefs[i][xi];
    }
    if (sum_i_xi > 0.0) {
      for (int xi=0; xi<ia_mrf->V[i]; xi++) {
	ia_beliefs[i][xi] /= sum_i_xi;
      }
    }
  }
  (*converged) = total; // actually irrelevant
  return ia_beliefs;
}
