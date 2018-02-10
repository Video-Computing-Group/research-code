#include "mex.h"
#include "fillMethods.h"
#include "PottsMRF.h"
#include "Metropolis.h"
#include "values.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // **************************************************************************
  // reading input arguments
  // **************************************************************************

  // check number of input arguments.
  // arguments should be:
  // adjMat - 1xN cell array, each cell {i} is a row vector with the indices of i's neighbours,
  // lambda - 1xN cell array, each cell {i} is a row vector with the strength of interaction of
  //          i with each of its neighbours
  //          note: Psi{i,j} = exp( [-lambda(i,j), +lambda(i,j); +lambda(i,j), -lambda(i,j)] )
  // local -  Nx1 cell array, each cell {i} is a row vector of length 2 containing the local
  //          potentials
  // temp, burning_time - both row vectors of the same length, representing the cooling schedule:
  //                      burning_time[i] is the number of steps to be taken at temperature temp[i]
  //
  // note: N = number of nodes
  
  if (nrhs != 5) {
    mexErrMsgTxt("Incorrect number of inputs.");
  }

  // check number of output arguments
  if (nlhs > 2) {
    mexErrMsgTxt("Too many output arguments.");
  }

  // get number of nodes and adjMat
  vector<Nodes>* adjMat = new vector<Nodes>();
  fillAdjMat(prhs[0],*adjMat);
  int num_nodes = adjMat->size();

  // define the MRF
  PottsMRF* mrf = new PottsMRF(*adjMat);
  
  // get local potentials
  fillLocalMat(prhs[2],mrf);
  
  // get pairwise potentials
  fillLambdaMat(prhs[1],mrf);

  // get cooling schedule:
  // 1. get temperatures
  int temp_nd = mxGetNumberOfDimensions(prhs[3]);
  const int* temp_dim = mxGetDimensions(prhs[3]);
  int sched_len = temp_dim[1];
  if (temp_nd != 2 || temp_dim[0] != 1) {
    mexErrMsgTxt("temp must be a row vector: 1x(schedule length)\n");
  }
  double* temp = new double[sched_len];
  double* temp_ptr = mxGetPr(prhs[3]);
  for (int i=0; i<sched_len; i++) {
    temp[i] = temp_ptr[i];
  }
  // 2. get burning times
  int burn_nd = mxGetNumberOfDimensions(prhs[4]);
  const int* burn_dim = mxGetDimensions(prhs[4]);
  if (burn_nd != 2 || burn_dim[0] != 1 || burn_dim[1] != sched_len) {
    mexErrMsgTxt("burning_time must be a row vector: 1x(schedule length)\n");
  }
  int* burn = new int[sched_len];
  double* burn_ptr = mxGetPr(prhs[4]);
  for (int i=0; i<sched_len; i++) {
    burn[i] = (int)(burn_ptr[i]);
  }


  // **************************************************************************
  // running metropolis in the schedule defined
  // **************************************************************************

  int* startX = new int[num_nodes];
  int* currX = new int[num_nodes];
  int* minX = new int[num_nodes];
  double currEnergy = 0.0;
  double minEnergy = MAXDOUBLE;
  for (int i=0; i<num_nodes; i++) {
    double randomal = 1.0 * rand()/(RAND_MAX);
    startX[i] =  (randomal > 0.5 ? 1 : 0);
  }

  Metropolis* metro = new Metropolis(mrf,startX,0,0,0);
  delete[] startX;
  startX = 0;

  for (int t=0; t<sched_len; t++) {
    mrf->setTemperature(temp[t]);
    metro->setBurningTime(burn[t]);
    metro->reachEquilibrium();
    metro->getCurrentState(currX);
    currEnergy = mrf->getEnergy(currX);
    if (currEnergy < minEnergy) {
      minEnergy = currEnergy;
      for (int i=0; i<num_nodes; i++) {
	minX[i] = currX[i];
      }
    }    
  }
  if (minEnergy < currEnergy) {
    mexPrintf("minEnergy = %.4f currEnergy = %.4f\n",minEnergy,currEnergy);
  }
  mexPrintf("energy found: %.4f\n",minEnergy);

  // **************************************************************************
  // writing output arguments
  // **************************************************************************

  if (nlhs > 0) {

    int state_dims[2] = {1, num_nodes};
    plhs[0] = mxCreateNumericArray(2, state_dims, mxDOUBLE_CLASS, mxREAL);
    double* state_ptr = mxGetPr(plhs[0]);
    for (int i=0; i<num_nodes; i++) {
      state_ptr[i] = minX[i] + 1.0;
    }
  }

  // **************************************************************************
  // free memory
  // **************************************************************************

  delete[] currX;
  currX = 0;
  delete[] minX;
  minX = 0;
  delete metro;
  metro = 0;
  delete mrf;
  mrf = 0;
  delete adjMat;
  adjMat = 0;
  
}
