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
  
  if (nrhs != 10) {
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
  int num_nodes = 0;
  int local_nd = mxGetNumberOfDimensions(prhs[1]);
  const int* local_dim = mxGetDimensions(prhs[1]);
  bool local_given = (local_nd == 2 && mxIsCell(prhs[1]));
  
  // get N
  num_nodes = local_dim[0]*local_dim[1];
  V = new int[num_nodes];
  local = new Potentials[num_nodes];

  if (local_given) {
    // get V and local potentials
    for (int i=0; i<num_nodes; i++) {
      mxArray* local_i = mxGetCell(prhs[1],i);
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
    // get V and initialize potentials to default uniform potentials
    double* VPtr = mxGetPr(prhs[1]);
    for (int i=0; i<num_nodes; i++) {
      V[i] = (int)(VPtr[i]);
      local[i] = new Potential[V[i]];
      for (int xi=0; xi<V[i]; xi++) {
	local[i][xi] = 1.0;
      }
    }
  }
//   fillLocalPot(prhs[1], local, V, num_nodes);


  int num_factors = factors->size();
  int num_all_regs = num_factors + num_nodes;
    
  // create the region-adj-mat  
  vector<Nodes>* adjMat = new vector<Nodes>();
  adjMat->resize(num_all_regs);
  for (int i=0; i<num_factors; i++) {
    Region& reg_i = (*factors)[i];
    nodes_iterator reg_niter = reg_i.begin();
    while (reg_niter != reg_i.end()) {
      int node = (*reg_niter);
      (*adjMat)[i].push_back(node+num_factors);
      (*adjMat)[node+num_factors].push_back(i);
      reg_niter++;
    }
  }
  // create the region-mrf
  MRF* reg_mrf = new MRF(*adjMat);

  // calculate the number of states
  for (int i=0; i<num_nodes; i++) { // for single nodes
    reg_mrf->V[i+num_factors] = V[i];
  }
  for (int i=0; i<num_factors; i++) { // for factors
    reg_mrf->V[i] = 1;
    Region& reg_i = (*factors)[i];
    nodes_iterator reg_niter = reg_i.begin();
    while (reg_niter != reg_i.end()) {
      int node = (*reg_niter);
      reg_mrf->V[i] *= V[node];
      reg_niter++;
    }
  }
  delete[] V;
  V = 0;

  // initialize regions' potentials
  reg_mrf->initLocalPotentials();
  
  // move the local potentials into the mrf
  for (int i=0; i<num_nodes; i++) {
    for (int xi=0; xi<reg_mrf->V[i+num_factors]; xi++) {
      reg_mrf->localMat[i+num_factors][xi] = local[i][xi];
    }
    delete[] local[i];
  }
  delete[] local;
  local = 0;
  
  // read the factor-potentials
  int fac_nd = mxGetNumberOfDimensions(prhs[2]);
  const int* fac_dim = mxGetDimensions(prhs[2]);
  if (fac_nd != 2 || fac_dim[0]*fac_dim[1] != num_factors || !mxIsCell(prhs[2])) {
    mexErrMsgTxt("factor-potentials must be a cell array in size 1x(num_factors), each cell {i} is a column vector in length (num_values(i))");
  }
  for (int i=0; i<num_factors; i++) {
    mxArray* facPot_i = mxGetCell(prhs[2],i);
    int i_nd = mxGetNumberOfDimensions(facPot_i);
    const int* i_dim = mxGetDimensions(facPot_i);
    if (i_nd != 2 || i_dim[0]*i_dim[1] != reg_mrf->V[i]) {
      mexErrMsgTxt("each cell {i} in the factor-potentials cell-array must be a vector in length (num_values(i))");
    }

    // get potentials
    double* facPotPtr = mxGetPr(facPot_i);
    for (int xi=0; xi<reg_mrf->V[i]; xi++) {
      reg_mrf->localMat[i][xi] = facPotPtr[xi];
    }
  }

  // calculate the assignment-indices-table
  int*** assignInd = new int**[num_factors];
  for (int i=0; i<num_factors; i++) {
    Region& reg_i = (*factors)[i];
    int i_numStates = reg_mrf->V[i];
    
    assignInd[i] = new int*[reg_i.size()];

    Assignment assignment;
    assignment.clear();
    assignment.resize(num_nodes);

    for (int n=0; n<reg_i.size(); n++) {
      int node = reg_i[n];
      assignInd[i][n] = new int[i_numStates];
      for (int i_assignInd=0; i_assignInd<i_numStates; i_assignInd++) {
	reg_i.indToAssign(i_assignInd, assignment, &(reg_mrf->V[num_factors]));
	assignInd[i][n][i_assignInd] = assignment[node];
      }
    }
  }

  
  // read the rest of parameters
  
  int maxIter = (int)(mxGetScalar(prhs[3]));

  SumOrMax sumOrMax = (SumOrMax)((int)(mxGetScalar(prhs[4])));
  double gbp_alpha = mxGetScalar(prhs[5]);

  bool trw = ((int)(mxGetScalar(prhs[6]))) > 0;
  double* countingNode = new double[num_nodes];
  fillDouble(prhs[7], countingNode, num_nodes);
  double* bethe = new double[num_all_regs];
  for (int i=0; i<num_factors; i++) {
    bethe[i] = 1.0;
  }
  for (int i=0; i<num_nodes; i++) {
    double q = (1.0 - countingNode[i]) / max(1.0, double(reg_mrf->neighbNum(i+num_factors)));
    bethe[i+num_factors] = 1.0 / (2.0 - q);
  }
  delete[] countingNode;
  countingNode = 0;

  // get initial messages, if given
  double*** initMsg = 0;
  int initM_nd = mxGetNumberOfDimensions(prhs[8]);
  const int* initM_dim = mxGetDimensions(prhs[8]);
  if ((initM_nd == 2) && mxIsCell(prhs[8])) {
    int num_regs = initM_dim[0] * initM_dim[1];
    initMsg = new double**[num_regs];
    for (int i=0; i<num_regs; i++) {

      mxArray* initMsg_i = mxGetCell(prhs[8],i);

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
  double temperature = mxGetScalar(prhs[9]);
  reg_mrf->setTemperature(temperature);
  
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
    converged = -2;
    mexWarnMsgTxt("resulted beliefs are not marginalizable\n");
  }
  
  // **************************************************************************
  // assign results to output argument
  // **************************************************************************
  if (nlhs > 0) {
    int nodes_dims[2] = {1,num_nodes};
    plhs[0] = mxCreateCellArray(2,nodes_dims); // single beliefs

    for (int i=0; i<num_nodes; i++) {

      // node beliefs
      int val_dims[2] = {reg_mrf->V[i+num_factors],1};
      mxArray* bel_i = mxCreateNumericArray(2,val_dims,mxDOUBLE_CLASS, mxREAL);
      double* resBelPtr = mxGetPr(bel_i);
      for (int xi=0; xi<reg_mrf->V[i+num_factors]; xi++) {
	resBelPtr[xi] = beliefs[i+num_factors][xi];
      }
      mxSetCell(plhs[0],i,bel_i);

    }

    if (nlhs > 1) {
      int factors_dims[2] = {1,num_factors};
      plhs[1] = mxCreateCellArray(2,factors_dims); // factor beliefs

      for (int i=0; i<num_factors; i++) {

	// factor beliefs
	int val_dims[2] = {reg_mrf->V[i],1};
	mxArray* bel_i = mxCreateNumericArray(2,val_dims,mxDOUBLE_CLASS, mxREAL);
	double* resBelPtr = mxGetPr(bel_i);
	for (int xi=0; xi<reg_mrf->V[i]; xi++) {
	  resBelPtr[xi] = beliefs[i][xi];
	}
	mxSetCell(plhs[1],i,bel_i);

      }

      if (nlhs > 2) {
	plhs[2] = mxCreateDoubleScalar(converged); // convergence flag
      }
    }
  } 

  // **************************************************************************
  // free memory
  // **************************************************************************
  delete algorithm;
  algorithm = 0;
  delete reg_mrf;
  reg_mrf = 0;
  delete adjMat;
  adjMat = 0;

  for (int i=0; i<num_factors; i++) {
    for (int n=0; n<(*factors)[i].size(); n++) {
      delete[] assignInd[i][n];
    }
    delete[] assignInd[i];
  }
  delete[] assignInd;
  assignInd = 0;
  
  delete factors;
  factors = 0;

  delete[] bethe;
  bethe = 0;
}


