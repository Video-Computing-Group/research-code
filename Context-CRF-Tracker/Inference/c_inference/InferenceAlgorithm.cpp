#include "InferenceAlgorithm.h"

InferenceAlgorithm::InferenceAlgorithm(MRF const* mrf)
{
  ia_mrf = mrf;

  ia_beliefs = 0;
  initBeliefs();
}

InferenceAlgorithm::~InferenceAlgorithm()
{
  freeBeliefs();
}

void InferenceAlgorithm::initBeliefs() {

  freeBeliefs();
  // init the beliefs defultive to p(Xi=xi) = Psi(xi,i)  
  ia_beliefs = new double*[ia_mrf->N];
  for (int i=0; i<ia_mrf->N; i++) {
    ia_beliefs[i] = new double[ia_mrf->V[i]];
    for (int xi=0; xi<ia_mrf->V[i]; xi++) {
      ia_beliefs[i][xi] = ia_mrf->localMat[i][xi];
    }
  }
}

void InferenceAlgorithm::freeBeliefs() {
  if (ia_beliefs != 0) {
    // free the beliefs-matrix and localMat
    for (int i=0; i<ia_mrf->N; i++) {
      delete[] ia_beliefs[i];
    }
    delete[] ia_beliefs;
    ia_beliefs = 0;
  }
}

void InferenceAlgorithm::getMaximalStates(vector<int>& max_states, int i, double epsilon) const {
  int max_xi = 0;
  for (int xi=0; xi<ia_mrf->V[i]; xi++) {
    if (ia_beliefs[i][xi] > ia_beliefs[i][max_xi]) {
      max_xi = xi;
    }
  }
  max_states.clear();
  for (int xi=0; xi<ia_mrf->V[i]; xi++) {
    if (ia_beliefs[i][xi] > (ia_beliefs[i][max_xi] - epsilon)) {
      max_states.push_back(xi);
    }
  }
}

