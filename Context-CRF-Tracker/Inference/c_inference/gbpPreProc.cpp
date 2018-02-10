#include "fillMethods.h"
#include "GBPPreProcessor.h"

// *****************************************************************************
// mexFunction
// *****************************************************************************

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // variables for GBP
  int*** assignInd = 0;
  double* bethe = 0;
  GBPPreProcessor* processor = 0;
  MRF* reg_mrf = 0;

  // **************************************************************************
  // reading input arguments
  // **************************************************************************

  if (nrhs != 4) {
    mexErrMsgTxt("Incorrect number of input arguments.");
  }

  // check number of output arguments
  if ((nlhs > 0) && (nlhs != 4)) {
    mexErrMsgTxt("Incorrect number of output arguments.");
  }

  // get number of nodes and adjMat
  vector<Nodes>* adjMat = new vector<Nodes>();
  fillAdjMat(prhs[0],*adjMat);
  int num_nodes = adjMat->size();

  // define the MRF
  MRF* mrf = new MRF(*adjMat);
  // get local potentials
  fillLocalMat(prhs[2],mrf);
  
  // get pairwise potentials
  fillPsiMat(prhs[1],mrf);
  
  // get regions for generalized belief propagation
  RegionLevel* regions = new RegionLevel();
  fillRegions(prhs[3],*regions);


  // **************************************************************************
  // making pre-process
  // **************************************************************************

  processor = new GBPPreProcessor(*regions, mrf);
  regions->clear();
  delete regions;
  regions = 0;
      
  reg_mrf = processor->getRegionMRF();
  assignInd = processor->getAssignTable();
  bethe = processor->getBethe();

  
  // **************************************************************************
  // writing output arguments
  // **************************************************************************

  // output: 1. regions' adj-matrix
  //         2. regions' local potentials
  //         3. assignment-indices (for each parent's assignment, the index of
  //                                the relevant son's assignment)
  //         4. bethe
  
  if (nlhs > 0) {

    // assign: 1. regions' adj-matrix & 2. regions' local potentials
    int num_regs = reg_mrf->N;
    int adj_local_n_bethe_dims[2] = {1,num_regs};
    plhs[0] = mxCreateCellArray(2,adj_local_n_bethe_dims);
    plhs[1] = mxCreateCellArray(2,adj_local_n_bethe_dims);

    for (int i=0; i<num_regs; i++) {

      // adj-matrix
      int adj_dims[2] = {reg_mrf->neighbNum(i),1};
      mxArray* adj_i = mxCreateNumericArray(2,adj_dims,mxDOUBLE_CLASS, mxREAL);
      double* adj_i_ptr = mxGetPr(adj_i);
      for (int n=0; n<reg_mrf->neighbNum(i); n++) {
	adj_i_ptr[n] = (double)(reg_mrf->adjMat[i][n] + 1);
      }
      mxSetCell(plhs[0],i,adj_i);

      // local potentials
      int local_dims[2] = {reg_mrf->V[i],1};
      mxArray* local_i = mxCreateNumericArray(2, local_dims, mxDOUBLE_CLASS, mxREAL);
      double* local_i_ptr = mxGetPr(local_i);
      for (int xi=0; xi<reg_mrf->V[i]; xi++) {
	local_i_ptr[xi] = reg_mrf->localMat[i][xi];
      }
      mxSetCell(plhs[1],i,local_i);

    }

    // assign: 3. assignment indices
    int ass_dims[2] = {num_regs,num_regs};
    plhs[2] = mxCreateCellArray(2,ass_dims);

    for (int i=0; i<num_regs; i++) {
      for (int n=0; n<reg_mrf->neighbNum(i); n++) {
	int j = reg_mrf->adjMat[i][n];
	if (i<j) {
	  int ass_val_dims[2] = {reg_mrf->V[i],1};
	  mxArray* ass_ij = mxCreateNumericArray(2, ass_val_dims, mxDOUBLE_CLASS, mxREAL);
	  double* ass_ij_ptr = mxGetPr(ass_ij);
	  for (int xi=0; xi<reg_mrf->V[i]; xi++) {
	    ass_ij_ptr[xi] = (double)(assignInd[i][n][xi] + 1);
	  }
	  mxSetCell(plhs[2],i + j*num_regs,ass_ij);
	}
      }
    }

    // assign: 4. bethe
    plhs[3] = mxCreateNumericArray(2, adj_local_n_bethe_dims, mxDOUBLE_CLASS, mxREAL);
    double* bethe_ptr = mxGetPr(plhs[3]);
    for (int i=0; i<num_regs; i++) {
      bethe_ptr[i] = bethe[i];
    }
      
  }
  
}
