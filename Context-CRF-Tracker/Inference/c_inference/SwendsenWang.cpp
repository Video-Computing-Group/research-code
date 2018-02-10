#include "SwendsenWang.h"
#include <stdlib.h>
#include "mex.h"


SwendsenWang::~SwendsenWang() {
  for (int i=0; i<ia_mrf->N; i++) {
    delete[] sw_bondFrozen[i];
  }

  delete[] sw_bondFrozen;
  sw_bondFrozen = 0;

  delete[] sw_cluster;
  sw_cluster = 0;

  delete[] sw_labelLabel;
  sw_labelLabel = 0;

  delete[] sw_new;
  sw_new = 0;

  delete[] sw_newChosen;
  sw_newChosen = 0;
  
}

void SwendsenWang::transition() {
  // construct a bond lattice with frozen bonds
  freezeOrMeltBonds();
  // use the Hoshen-Kopelman algorithm to identify and label clusters
  labelClusters();
  // re-set cluster nodes randomly
  flipClusterNodes();
}

void SwendsenWang::init() {
  initializeClusterVariables();
}

void SwendsenWang::initializeClusterVariables() {
  sw_bondFrozen = new bool*[ia_mrf->N];
  for (int i=0; i<ia_mrf->N; i++) {
    sw_bondFrozen[i] = new bool[ia_mrf->neighbNum(i)];
  }
  sw_cluster = new int[ia_mrf->N];
  sw_labelLabel = new int[ia_mrf->N];
  sw_new = new int[ia_mrf->N];
  sw_newChosen = new bool[ia_mrf->N];
}

void SwendsenWang::freezeOrMeltBonds() {
  for (int i=0; i<ia_mrf->N; i++) {
    int xi = mc_currX[i];
    for (int n_ij=0; n_ij<ia_mrf->neighbNum(i); n_ij++) {
      int j = ia_mrf->adjMat[i][n_ij];
      if (i<j) {
	int n_ji = 0;
	while (ia_mrf->adjMat[j][n_ji] != i) {
	  ++n_ji;
	}
	bool freeze = false;
	if (xi == mc_currX[j]) {
	  // freeze the bond <i,j> in probability 1 - exp(-Jij/T)
	  double randomal = (1.0 * rand()/(RAND_MAX));
	  double psi_ij = ia_mrf->pairPotential(i,n_ij,xi,xi);
	  double freezePr = 1.0 - (1.0 / psi_ij);
	  freeze = (randomal <= freezePr);
	}
	sw_bondFrozen[i][n_ij] = freeze;
	sw_bondFrozen[j][n_ji] = freeze;
      }
    }
  }

}

int SwendsenWang::properLabel(int label) {
  while (sw_labelLabel[label] != label) {
    label = sw_labelLabel[label];
  }
  return label;
}


// implementation of the Hoshen-Kopelman algorithm

void SwendsenWang::labelClusters() {
  int label = 0;
  for (int i=0; i<ia_mrf->N; i++) {
    // check number of bonds to previously visited sites and
    // find the smallest proper label
    int bonds = 0;
    int minLabel = label;
    for (int n=0; n<ia_mrf->neighbNum(i); n++) {
      int j = ia_mrf->adjMat[i][n];
      if (j<i && sw_bondFrozen[i][n]) {
	++bonds;
	int pLabel = properLabel(sw_cluster[j]);
	if (pLabel<minLabel) {
	  minLabel = pLabel;
	}
      }
    }
    if (bonds==0) { // need to start a new cluster
      sw_cluster[i] = label;
      sw_labelLabel[label] = label;
      ++label;
    }
    else { // re-label bonded nodes with smallest proper label

      // set current site label to smallest proper label
      sw_cluster[i] = minLabel;

      // re-set the proper label links on the previous labels
      for (int n=0; n<ia_mrf->neighbNum(i); n++) {
	int j = ia_mrf->adjMat[i][n];
	if (j<i && sw_bondFrozen[i][n]) {
	  int pLabel = sw_cluster[j];
	  sw_labelLabel[pLabel] = minLabel;
	  // re-set label on connected sites
	  sw_cluster[j] = minLabel;
	}
      }
    }
  }
}

void SwendsenWang::flipClusterNodes() {
  for (int i=0; i<ia_mrf->N; i++) {
    // random new cluster nodes values have not been set
    sw_newChosen[i] = false;
    // replace all labels by their proper values
    sw_cluster[i] = properLabel(sw_cluster[i]);
  }

  int flips = 0;
  for (int i=0; i<ia_mrf->N; i++) {
    int label = sw_cluster[i];
    // choose a random new value for the cluster,
    // if it has not been done yet
    if (!sw_newChosen[label]) {
      sw_new[label] = chooseInteger(mc_Pxi[i], ia_mrf->V[i]);
      sw_newChosen[label] = true;
    }
    if (mc_currX[i] != sw_new[label]) {
      mc_currX[i] = sw_new[label];
      ++flips;
    }
  }

}

