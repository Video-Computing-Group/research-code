#include "Gibbs.h"
#include <stdlib.h>
#include <iostream>

void Gibbs::transition() {

  int i = chooseInteger(mc_Pi, ia_mrf->N);
  int qi = ia_mrf->V[i];
  for (int xi=0; xi<qi; xi++) {
    mc_Pxi[i][xi] = ia_mrf->localMat[i][xi];
  }
  
  for (int n=0; n<ia_mrf->neighbNum(i); n++) {
    int j = ia_mrf->adjMat[i][n];
    int xj = mc_currX[j];
    for (int xi=0; xi<qi; xi++) {
      mc_Pxi[i][xi] *= ia_mrf->pairPotential(i,n,xi,xj);
    }
  }

  // normalize Pxi
  double sum_Pxi = 0.0;
  for (int xi=0; xi<qi; xi++) {
    sum_Pxi += mc_Pxi[i][xi];
  }
  for (int xi=0; xi<qi; xi++) {
    mc_Pxi[i][xi] /= sum_Pxi;
  }
  
  mc_currX[i] = chooseInteger(mc_Pxi[i], qi);
  
}
