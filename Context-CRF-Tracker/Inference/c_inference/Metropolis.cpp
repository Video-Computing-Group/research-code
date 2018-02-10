#include "Metropolis.h"
#include <math.h>

void Metropolis::transition() {

  int i = chooseInteger(mc_Pi, ia_mrf->N);
  int xi = mc_currX[i];

  int qi = ia_mrf->V[i];
  double equalPr = 1.0 / (qi - 1.0);
  for (int xj=0; xj<qi; xj++) {
    mc_Pxi[i][xj] = equalPr;
  }
  mc_Pxi[i][xi] = 0.0;
  int new_xi = chooseInteger(mc_Pxi[i], qi);

  // calculate the probabilty of the flip:
  // flipPr = exp(-dE/T) where dE is the delta in energy resulted
  // by the flip, and dE>0. if dE<=0, the flip is done in prob. 1
  double mult_psi_ij_new = ia_mrf->localMat[i][new_xi];
  double mult_psi_ij_old = ia_mrf->localMat[i][xi];
  for (int n=0; n<ia_mrf->neighbNum(i); n++) {
    int j = ia_mrf->adjMat[i][n];
    int xj = mc_currX[j];
    mult_psi_ij_new *= ia_mrf->pairPotential(i,n,new_xi,xj);
    mult_psi_ij_old *= ia_mrf->pairPotential(i,n,xi,xj);
  }
  double flipPr = min(1.0, mult_psi_ij_new / mult_psi_ij_old); // this is exp(-dE/T) or 1.0 for dE<=0
  double randomal = (1.0 * rand()/(RAND_MAX));
  if (randomal <= flipPr) {
    mc_currX[i] = new_xi;
  }
}
