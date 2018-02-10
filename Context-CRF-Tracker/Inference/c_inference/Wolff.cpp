#include "Wolff.h"
#include <stdlib.h>
#include <iostream>

void Wolff::transition() {

  int i = chooseInteger(mc_Pi, ia_mrf->N);

  int new_xi = chooseInteger(mc_Pxi[i], ia_mrf->V[i]);

  int cluster_size = 0;
  if (mc_currX[i] != new_xi) {
    cluster_size = cluster(i, i, new_xi);
  }

}

int Wolff::cluster(int i, int i_nei, int new_xi) {
  int xi = mc_currX[i];
  int size = 0;
  if (i==i_nei || accept(i, i_nei)) {
    size++;
    mc_currX[i] = new_xi;
    // check the neighbours
    for (int n=0; n<ia_mrf->neighbNum(i); n++) {
      int j = ia_mrf->adjMat[i][n];
      int xj = mc_currX[j];
      if (xi==xj) {
	size += cluster(j, i, new_xi);
      }
    }
  }
  return size;
}

bool Wolff::accept(int i, int i_nei) {
  int xi = mc_currX[i];
  double randomal = (1.0 * rand()/(RAND_MAX));
  int ni = 0;
  while (ia_mrf->adjMat[i][ni] != i_nei) {
    ++ni;
  }
  double psi_ij = ia_mrf->pairPotential(i,ni,xi,xi);
  double Padd = 1.0 - (1.0 / psi_ij);
  return (randomal <= Padd);
}
