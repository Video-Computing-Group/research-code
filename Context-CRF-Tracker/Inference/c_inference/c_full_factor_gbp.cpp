#include "mex.h"
#include "fillMethods.h"
#include "LogLoopy.h"
#include "GBP.h"
#include "Gibbs.h"
#include "Wolff.h"
#include "SwendsenWang.h"
#include "Metropolis.h"
#include "MeanField.h"
#include "GBPPreProcessor.h"
#include <iostream>

using namespace std;

// *****************************************************************************
// mexFunction
// *****************************************************************************

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  // **************************************************************************
  // reading input arguments
  // **************************************************************************
  
  if (nrhs != 12) {
    mexErrMsgTxt("Incorrect number of inputs.");
  }

  // check number of output arguments
  if (nlhs > 3) {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  // get factors
  RegionLevel* factors = new RegionLevel();
  fillRegions(prhs[0],*factors);

  // get number of nodes and local potentials
  Potentials* local = 0;
  int* V = 0;
  // get N
  int num_nodes = (int)(mxGetScalar(prhs[1]));
  int local_nd = mxGetNumberOfDimensions(prhs[2]);
  const int* local_dim = mxGetDimensions(prhs[2]);
  bool local_given = (local_nd == 2 && mxIsCell(prhs[2]));
  
  V = new int[num_nodes];

  if (local_given) {
    if (num_nodes != local_dim[0]*local_dim[1]) {
      mexErrMsgTxt("number of nodes inconsistent with size of local potentials\n");
    }
    // get V and local potentials
    local = new Potentials[num_nodes];
    for (int i=0; i<num_nodes; i++) {
      mxArray* local_i = mxGetCell(prhs[2],i);
      int i_nd = mxGetNumberOfDimensions(local_i);
      const int* i_dim = mxGetDimensions(local_i);
      if (i_nd != 2 || i_dim[1] != 1) {
	mexErrMsgTxt("each cell {i} in the local cell-array must be a column vector (num_values(i))x1");
      }

      // get V
      V[i] = i_dim[0];
      // get local
      local[i] = new Potential[V[i]];
      double* localPtr = mxGetPr(local_i);
      for (int xi=0; xi<V[i]; xi++) {
	local[i][xi] = localPtr[xi];
      }
    
    }
  }
  else {
    // get V
    if (local_dim[0]*local_dim[1] > 1) {
      if (num_nodes != local_dim[0]*local_dim[1]) {
	mexErrMsgTxt("number of nodes inconsistent with size of cardinalities\n");
      }
      double* VPtr = mxGetPr(prhs[2]);
      for (int i=0; i<num_nodes; i++) {
	V[i] = (int)(VPtr[i]);
      }
    }
    else {
      int card = (int)(mxGetScalar(prhs[2]));
      for (int i=0; i<num_nodes; i++) {
	V[i] = card;
      }
    }
  }

  // create adjacencies matrix with edges between nodes in the same factor
  int num_factors = factors->size();
  vector<Nodes>* adjMat = new vector<Nodes>();
  adjMat->resize(num_nodes);
  for (int r=0; r<num_factors; r++) {
    Region& reg = (*factors)[r];
    for (int i=0; i<reg.size(); i++) {
      for (int j=0; j<reg.size(); j++) {
	if (i!= j) {
	  bool insert = true;
	  for (int k=0; k<(*adjMat)[reg[i]].size(); k++) {
	    if ((*adjMat)[reg[i]][k] == reg[j]) {
	      insert = false;
	      break;
	    }
	  }
	  if (insert) {
	    (*adjMat)[reg[i]].push_back(reg[j]);
	  }
	}
      }
    }
  }

  MRF* mrf = new MRF(*adjMat);
  for (int i=0; i<num_nodes; i++) {
    mrf->V[i] = V[i];
  }
  delete[] V;

  if (local != 0) {
    mrf->initLocalPotentials();

    // fill localMat
    for (int i=0; i<num_nodes; i++) {
      for (int xi=0; xi<mrf->V[i]; xi++) {
	mrf->localMat[i][xi] = local[i][xi];
      }
      delete[] local[i];
    }
    delete[] local;
  }
  // read the factor-potentials

  Potentials* factorPot = new Potentials[num_factors];
  
  int fac_nd = mxGetNumberOfDimensions(prhs[3]);
  const int* fac_dim = mxGetDimensions(prhs[3]);
  if (fac_nd != 2 || fac_dim[0]*fac_dim[1] != num_factors || !mxIsCell(prhs[3])) {
    mexErrMsgTxt("factor-potentials must be a cell array in size 1x(num_factors), each cell {i} is a column vector in length (num_values(i))");
  }
  for (int i=0; i<num_factors; i++) {
    mxArray* facPot_i = mxGetCell(prhs[3],i);
    const int* i_dim = mxGetDimensions(facPot_i);

    int num_states = i_dim[0] * i_dim[1];
    factorPot[i] = new Potential[num_states];

    // get potentials
    fillDouble(facPot_i,factorPot[i],num_states);
  }

  // read the rest of parameters
  
  int maxIter = (int)(mxGetScalar(prhs[4]));

  SumOrMax sumOrMax = (SumOrMax)((int)(mxGetScalar(prhs[5])));
  double gbp_alpha = mxGetScalar(prhs[6]);

  bool trw = ((int)(mxGetScalar(prhs[7]))) > 0;
  double* countingNode = new double[num_nodes];
  fillDouble(prhs[8], countingNode, num_nodes);

  bool full = ((int)(mxGetScalar(prhs[9]))) > 0;
  
  // get initial messages, if given
  double*** initMsg = 0;
  int initM_nd = mxGetNumberOfDimensions(prhs[10]);
  const int* initM_dim = mxGetDimensions(prhs[10]);
  if ((initM_nd == 2) && mxIsCell(prhs[10])) {
    int num_regs = initM_dim[0] * initM_dim[1];
    initMsg = new double**[num_regs];
    for (int i=0; i<num_regs; i++) {

      mxArray* initMsg_i = mxGetCell(prhs[10],i);

      int msg_i_nd = mxGetNumberOfDimensions(initMsg_i);
      const int* msg_i_dim = mxGetDimensions(initMsg_i);
      if ((msg_i_nd != 2) || !mxIsCell(initMsg_i)) {
	mexErrMsgTxt("each cell {i} in initMsg for GBP should be a cell array in length of number of neighbour-regions to region i\n");
      }
      int Ni = msg_i_dim[0] * msg_i_dim[1];
      initMsg[i] = new double*[Ni];
      for (int n=0; n<Ni; n++) {
	mxArray* initMsg_ij = mxGetCell(initMsg_i,n);
	const int* msg_ij_dim = mxGetDimensions(initMsg_ij);
	int numStates = msg_ij_dim[0] * msg_ij_dim[1];
	initMsg[i][n] = new double[numStates];
	fillDouble(initMsg_ij,initMsg[i][n],numStates);
      }
    }
  }

  // get tepmerature
  double temperature = mxGetScalar(prhs[11]);
  mrf->setTemperature(temperature);
  
  // **************************************************************************
  // pre-process
  // **************************************************************************
  GBPPreProcessor* processor = new GBPPreProcessor(*factors, mrf, trw, full, countingNode, factorPot);
  MRF* reg_mrf = processor->getRegionMRF();
  int*** assignInd = processor->getAssignTable();
  double* bethe = processor->getBethe();

  // **************************************************************************
  // create the algorithm
  // **************************************************************************
  GBP* algorithm = new GBP(reg_mrf,assignInd,bethe,sumOrMax,gbp_alpha,maxIter,initMsg);

  // **************************************************************************
  // make inference
  // **************************************************************************
  int converged;
  double** beliefs = algorithm->inference(&converged);
  bool marg;
  if (sumOrMax==SUM) {
    marg = ((GBP*)algorithm)->isSumMarg(.0001);
  }
  else {
    marg = ((GBP*)algorithm)->isMaxMarg(.0001);
  }
  if (!marg) {
    converged = (converged>0 ? -2 : -1);
    mexWarnMsgTxt("resulted beliefs are not marginalizable\n");
  }

  // extract single beliefs
  double** singleBeliefs = new double*[num_nodes];
  for (int i=0; i<num_nodes; i++) {
    singleBeliefs[i] = new double[mrf->V[i]];      
  }
  processor->extractSingle(beliefs,singleBeliefs,sumOrMax);
  // **************************************************************************
  // assign results to output argument
  // **************************************************************************
  if (nlhs > 0) {
    int nodes_dims[2] = {1,num_nodes};
    plhs[0] = mxCreateCellArray(2,nodes_dims); // single beliefs

    for (int i=0; i<num_nodes; i++) {

      // node beliefs
      int val_dims[2] = {mrf->V[i],1};
      mxArray* bel_i = mxCreateNumericArray(2,val_dims,mxDOUBLE_CLASS, mxREAL);
      double* resBelPtr = mxGetPr(bel_i);
      for (int xi=0; xi<mrf->V[i]; xi++) {
	resBelPtr[xi] = singleBeliefs[i][xi];
      }
      mxSetCell(plhs[0],i,bel_i);

    }
    if (nlhs > 1) {
      plhs[1] = mxCreateDoubleScalar(converged); // convergence flag

      if (nlhs > 2) {
	int factors_dims[2] = {1,num_factors};
	plhs[2] = mxCreateCellArray(2,factors_dims); // factor beliefs

	for (int i=0; i<num_factors; i++) {

	  // factor beliefs
	  int val_dims[2] = {reg_mrf->V[i],1};
	  mxArray* bel_i = mxCreateNumericArray(2,val_dims,mxDOUBLE_CLASS, mxREAL);
	  double* resBelPtr = mxGetPr(bel_i);
	  for (int xi=0; xi<reg_mrf->V[i]; xi++) {
	    resBelPtr[xi] = beliefs[i][xi];
	  }
	  mxSetCell(plhs[2],i,bel_i);

	}
      }
    }
  } 

  // **************************************************************************
  // free memory
  // **************************************************************************
  delete algorithm;
  algorithm = 0;
  
  if (singleBeliefs != 0) {
    for (int i=0; i<num_nodes; i++) {
      delete[] singleBeliefs[i];
    }
    delete[] singleBeliefs;    
    singleBeliefs = 0;
  }
  if (processor != 0) {
    delete processor;
    processor = 0;
  }

  delete factors;
  factors = 0;

  delete[] countingNode;
  countingNode = 0;

  delete mrf;
  mrf = 0;
  delete adjMat;
  adjMat = 0;
}


