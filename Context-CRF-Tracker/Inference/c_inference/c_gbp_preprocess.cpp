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
  int* extractSingle = 0;
  int** extractPairs = 0;
  RegionLevel* regions = 0;
  vector<RegionLevel>* in_allRegions = 0;
  // variables for user
  Region** all_regions = 0;

  // **************************************************************************
  // reading input arguments
  // **************************************************************************

  if (nrhs != 8) {
    mexErrMsgTxt("Incorrect number of input arguments.");
  }

  // check number of output arguments
  if ((nlhs > 0) && (nlhs != 7)) {
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
  
  // get regions
  bool allLevels = (int)(mxGetScalar(prhs[4])) > 0;
  if (allLevels) {
    in_allRegions = new vector<RegionLevel>();
    in_allRegions->clear();
    fillRegionLevels(prhs[3],*in_allRegions);
  }
  else {
    regions = new RegionLevel();
    fillRegions(prhs[3],*regions);
  }

  bool trw = (int)(mxGetScalar(prhs[5])) > 0;
  double* countingNode = 0;
  if (trw) {
    countingNode = new double[num_nodes];
    fillDouble(prhs[6], countingNode, num_nodes);
  }

  bool full = (int)(mxGetScalar(prhs[7])) > 0;


  // **************************************************************************
  // making pre-process
  // **************************************************************************

  if (allLevels) {
    processor = new GBPPreProcessor(in_allRegions, mrf, trw, full, countingNode);
  }
  else {
    processor = new GBPPreProcessor(*regions, mrf, trw, full, countingNode);
    regions->clear();
    delete regions;
    regions = 0;
  }

  all_regions = processor->getAllRegions();
  reg_mrf = processor->getRegionMRF();
  assignInd = processor->getAssignTable();
  bethe = processor->getBethe();
  extractSingle = processor->getExtractSingle();
  extractPairs = processor->getExtractPairs();

  
  // **************************************************************************
  // writing output arguments
  // **************************************************************************

  // output: 1. regions - a cell array with vectors of the nodes in each region
  //         2. regions' adj-matrix
  //         3. regions' local potentials
  //         4. assignment-indices (for each parent's assignment, the index of
  //                                the relevant son's assignment)
  //         5. bethe
  //         6. extractSingle (for each node, from which region its belief
  //                           should be extracted)
  //         7. extractPairs  (for each pair of neighbours, from which region
  //                           its pairwise beliefs should be extracted)
  
  if (nlhs > 0) {

    int num_regs = reg_mrf->N;
    int regs_dims[2] = {1,num_regs};
    
    // assign: 1. regions 2. regions' adj-matrix & 3. regions' local potentials
    plhs[0] = mxCreateCellArray(2,regs_dims);
    plhs[1] = mxCreateCellArray(2,regs_dims);
    plhs[2] = mxCreateCellArray(2,regs_dims);

    for (int i=0; i<num_regs; i++) {

      // regions
      Region* region = all_regions[i];
      int reg_i_size = (int)(region->size());
      int reg_i_dims[2] = {1,reg_i_size};
      mxArray* reg_i = mxCreateNumericArray(2,reg_i_dims,mxDOUBLE_CLASS, mxREAL);
      double* reg_i_ptr = mxGetPr(reg_i);
      for (int n=0; n<reg_i_size; n++) {
	reg_i_ptr[n] = (double)((*region)[n]);
      }
      mxSetCell(plhs[0],i,reg_i);
      
      // adj-matrix
      int adj_dims[2] = {1,reg_mrf->neighbNum(i)};
      mxArray* adj_i = mxCreateNumericArray(2,adj_dims,mxDOUBLE_CLASS, mxREAL);
      double* adj_i_ptr = mxGetPr(adj_i);
      for (int n=0; n<reg_mrf->neighbNum(i); n++) {
	adj_i_ptr[n] = (double)(reg_mrf->adjMat[i][n] + 1);
      }
      mxSetCell(plhs[1],i,adj_i);

      // local potentials
      int local_dims[2] = {reg_mrf->V[i],1};
      mxArray* local_i = mxCreateNumericArray(2, local_dims, mxDOUBLE_CLASS, mxREAL);
      double* local_i_ptr = mxGetPr(local_i);
      for (int xi=0; xi<reg_mrf->V[i]; xi++) {
	local_i_ptr[xi] = reg_mrf->localMat[i][xi];
      }
      mxSetCell(plhs[2],i,local_i);

    }

    // assign: 4. assignment indices
    int ass_dims[2] = {1,num_regs};
    plhs[3] = mxCreateCellArray(2,ass_dims);

    for (int i=0; i<num_regs; i++) {
      
      int ass_i_dims[2] = {1,reg_mrf->neighbNum(i)};
      mxArray* ass_i = mxCreateCellArray(2,ass_i_dims);
      
      for (int n=0; n<reg_mrf->neighbNum(i); n++) {
	int j = reg_mrf->adjMat[i][n];
	if (i<j) {
	  int ass_val_dims[2] = {1,reg_mrf->V[i]};
	  mxArray* ass_ij = mxCreateNumericArray(2, ass_val_dims, mxDOUBLE_CLASS, mxREAL);
	  double* ass_ij_ptr = mxGetPr(ass_ij);
	  for (int xi=0; xi<reg_mrf->V[i]; xi++) {
	    ass_ij_ptr[xi] = (double)(assignInd[i][n][xi]);
	  }
	  mxSetCell(ass_i, n, ass_ij);
	}
      }
      
      mxSetCell(plhs[3], i, ass_i);
    }

    // assign: 5. bethe
    plhs[4] = mxCreateNumericArray(2, regs_dims, mxDOUBLE_CLASS, mxREAL);
    double* bethe_ptr = mxGetPr(plhs[4]);
    for (int i=0; i<num_regs; i++) {
      bethe_ptr[i] = bethe[i];
    }

    // assign: 6. extractSingle
    int extract_dims[2] = {1,num_nodes};

    plhs[5] = mxCreateNumericArray(2, extract_dims, mxDOUBLE_CLASS, mxREAL);
    double* exSing_ptr = mxGetPr(plhs[5]);
    for (int i=0; i<num_nodes; i++) {
      exSing_ptr[i] = extractSingle[i];
    }
    
    // assign: 7. extractPairs
    plhs[6] = mxCreateCellArray(2, extract_dims);

    for (int i=0; i<num_nodes; i++) {

      int exPair_i_dim[2] = {1,mrf->neighbNum(i)};
      mxArray* exPair_i = mxCreateNumericArray(2, exPair_i_dim, mxDOUBLE_CLASS, mxREAL);
      double* exPair_i_ptr = mxGetPr(exPair_i);
      for (int n=0; n<mrf->neighbNum(i); n++) {
	exPair_i_ptr[n] = extractPairs[i][n];
      }

      mxSetCell(plhs[6], i, exPair_i);
    }
    
  }

  // **************************************************************************
  // free memory
  // **************************************************************************

  delete mrf;
  mrf = 0;
  delete adjMat;
  adjMat = 0;
  if (countingNode != 0) {
    delete[] countingNode;
    countingNode = 0;
  }
}
