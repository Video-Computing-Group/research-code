#include "fillMethods.h"

// *****************************************************************************
// fillAdjMat
// *****************************************************************************

void fillAdjMat(const mxArray* mxAdj, vector<Nodes>& adjMat) {
  // get adj
  int adj_nd = mxGetNumberOfDimensions(mxAdj);
  const int* adj_dim = mxGetDimensions(mxAdj);
  if (adj_nd != 2 || adj_dim[0] != 1 || !mxIsCell(mxAdj)) {
    mexErrMsgTxt("adjMat must be 1x(num_nodes) cell array");
  }
  // get N
  int N = adj_dim[1];
  adjMat.resize(N);
  
  for (int i=0; i<N; i++) {
    mxArray* adjMx_i = mxGetCell(mxAdj,i);
    int adj_i_nd = mxGetNumberOfDimensions(adjMx_i);
    const int* adj_i_dim = mxGetDimensions(adjMx_i);
    if (adj_i_nd != 2 || adj_i_dim[0] != 1) {
      mexErrMsgTxt("each cell {i} in adjMat must be a row vector 1x(num_of_i_neighb)");
    }
    int neighNum_i = adj_i_dim[1];
    adjMat[i].resize(neighNum_i);
    double* adj_i_Ptr = mxGetPr(adjMx_i);
    for (int n=0; n<neighNum_i; n++) {
      adjMat[i][n] = (int)(adj_i_Ptr[n]) - 1;
    }
  }
}

// *****************************************************************************
// fillLocalMat
// *****************************************************************************

void fillLocalMat(const mxArray* mxLocal, MRF* mrf) {
  int const& N = mrf->N;
  
  int local_nd = mxGetNumberOfDimensions(mxLocal);
  const int* local_dim = mxGetDimensions(mxLocal);
  if (local_nd != 2 || local_dim[1] != N || !mxIsCell(mxLocal)) {
    mexErrMsgTxt("local must be a cell array in size 1x(num_nodes), each cell {i} is a column vector in length (num_values(i))\n");
  }

  // get V
  for (int i=0; i<N; i++) {
    mxArray* local_i = mxGetCell(mxLocal,i);
    int i_nd = mxGetNumberOfDimensions(local_i);
    const int* i_dim = mxGetDimensions(local_i);
    if (i_nd != 2 || i_dim[1] != 1) {
      mexErrMsgTxt("each cell {i} in the local cell-array must be a column vector (num_values(i))x1");
    }

    // get V
    mrf->V[i] = i_dim[0];
  }

  mrf->initLocalPotentials();

  // fill localMat
  for (int i=0; i<N; i++) {
    mxArray* local_i = mxGetCell(mxLocal,i);
    double* localPtr = mxGetPr(local_i);
    for (int xi=0; xi<mrf->V[i]; xi++) {
      mrf->localMat[i][xi] = localPtr[xi];
    }
  }
}

// *****************************************************************************
// fillLocalPot
// *****************************************************************************

void fillLocalPot(const mxArray* mxLocal, Potentials* local, int* V, int& N) {

  int local_nd = mxGetNumberOfDimensions(mxLocal);
  const int* local_dim = mxGetDimensions(mxLocal);
  if (local_nd != 2 || !mxIsCell(mxLocal)) {
    mexErrMsgTxt("local must be a cell array in size 1x(num_nodes), each cell {i} is a column vector in length (num_values(i))");
  }
  // get N
  N = local_dim[1];
  V = new int[N];
  local = new Potentials[N];

  // get V and local potentials
  for (int i=0; i<N; i++) {
    mxArray* local_i = mxGetCell(mxLocal,i);
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

// *****************************************************************************
// fillPsiMat
// *****************************************************************************

void fillPsiMat(const mxArray* mxLambda, MRF* mrf) {

  mrf->initPairPotentials();
  
  int const& N = mrf->N;
  int const* V = mrf->V;
  
  int lambda_nd = mxGetNumberOfDimensions(mxLambda);
  const int* lambda_dim = mxGetDimensions(mxLambda);
  if (lambda_nd != 2 ||
      lambda_dim[0] != 1 || lambda_dim[1] != N ||
      !mxIsCell(mxLambda)) {
    mexWarnMsgTxt("lambda must be a cell array in size: 1x(num_nodes)\nuniform potentials are taken as a default.");
    for (int i=0; i<N; i++) {
      for (int n=0; n<mrf->neighbNum(i); n++) {
	int j = mrf->adjMat[i][n];
	if (i<j) {
	  PairPotentials pairPot_ij = new Potentials[V[i]];
	  for (int xi=0; xi<V[i]; xi++) {
	    pairPot_ij[xi] = new Potential[V[j]];
	    for (int xj=0; xj<V[j]; xj++) {
	      pairPot_ij[xi][xj] = 1.0;
	    }
	  }
	  mrf->assignPairPotential(i, n, pairPot_ij);
	  for (int xi=0; xi<V[i]; xi++) {
	    delete[] pairPot_ij[xi];
	  }
	  delete[] pairPot_ij;
	  pairPot_ij = 0;
	}	
      }
    }
  }
  else {
    // fill lambdaMat
    for (int i=0; i<N; i++) {
      mxArray* lambda_i = mxGetCell(mxLambda, i);
      int i_nd = mxGetNumberOfDimensions(lambda_i);
      const int* i_dim = mxGetDimensions(lambda_i);
      if (i_nd != 2 ||
	  i_dim[0] != 1 || i_dim[1] != mrf->neighbNum(i) ||
	  !mxIsCell(lambda_i)) {
	mexErrMsgTxt("each cell {i} in lambda must be a cell array in size: 1x(neighb_num(i))\n");
      }

      for (int n=0; n<mrf->neighbNum(i); n++) {
	int j = mrf->adjMat[i][n];
	if (i<j) {
	  mxArray* lambda_ij = mxGetCell(lambda_i, n);
	  int ij_nd = mxGetNumberOfDimensions(lambda_ij);
	  const int* ij_dim = mxGetDimensions(lambda_ij);
	  if (ij_nd != 2 || ij_dim[0] != V[i] || ij_dim[1] != V[j]) {
	    mexErrMsgTxt("each cell {n} insize cell {i} in lambda, must be 2-dimensional matrix: (num_values(i))x(num_values(j)), where j is the n-th neighbour of i\n");
	  }
      
	  double* lambda_ij_ptr = mxGetPr(lambda_ij);

	  PairPotentials pairPot_ij = new Potentials[V[i]];
	  for (int xi=0; xi<V[i]; xi++) {
	    pairPot_ij[xi] = new Potential[V[j]];
	    for (int xj=0; xj<V[j]; xj++) {
	      pairPot_ij[xi][xj] = lambda_ij_ptr[xi + xj*V[i]];
	    }
	  }
	  mrf->assignPairPotential(i, n, pairPot_ij);
	  for (int xi=0; xi<V[i]; xi++) {
	    delete[] pairPot_ij[xi];
	  }
	  delete[] pairPot_ij;
	  pairPot_ij = 0;
	}
      }
    }
  }
}

// *****************************************************************************
// fillLambdaMat
// *****************************************************************************

void fillLambdaMat(const mxArray* mxLambda, PottsMRF* mrf) {
  
  mrf->initLambdaValues();
  
  int const& N = mrf->N;
  
  int lambda_nd = mxGetNumberOfDimensions(mxLambda);
  const int* lambda_dim = mxGetDimensions(mxLambda);
  if (lambda_nd != 2 ||
      lambda_dim[0] != 1 || lambda_dim[1] != N ||
      !mxIsCell(mxLambda)) {
    mexWarnMsgTxt("lambda must be a cell array in size: 1x(num_nodes)\nuniform potentials are taken as a default.");
    for (int i=0; i<N; i++) {
      for (int n=0; n<mrf->neighbNum(i); n++) {
	mrf->lambdaValues[i][n] = 1.0;
      }
    }
  }
  else {
    // fill lambdaMat
    for (int i=0; i<N; i++) {
      mxArray* lambda_i = mxGetCell(mxLambda,i);
      int i_nd = mxGetNumberOfDimensions(lambda_i);
      const int* i_dim = mxGetDimensions(lambda_i);
      if (i_nd != 2 || i_dim[0] != 1 || i_dim[1] != mrf->neighbNum(i)) {
	mexErrMsgTxt("each cell {i} in lambda must be a row vector: 1x(num_neighbours(i))\n");
      }

      double* lambda_i_ptr = mxGetPr(lambda_i);
      for (int n=0; n<mrf->neighbNum(i); n++) {
	mrf->lambdaValues[i][n] = lambda_i_ptr[n];
      }
    }
  }
}

// *****************************************************************************
// fillRegions
// *****************************************************************************

void fillRegions(const mxArray* mxRegions, RegionLevel& regions) {
  int reg_nd = mxGetNumberOfDimensions(mxRegions);
  const int* reg_dim = mxGetDimensions(mxRegions);
  if (reg_nd != 2 || reg_dim[0] != 1 || !mxIsCell(mxRegions)) {
    mexErrMsgTxt("regions must be a row cell array: 1x(num_regions)");
  }
  // fill regions
  int numRegs = reg_dim[1];
  for (int i=0; i<numRegs; i++) {
    mxArray* region_i = mxGetCell(mxRegions,i);
    int i_nd = mxGetNumberOfDimensions(region_i);
    const int* i_dim = mxGetDimensions(region_i);
    if (i_nd != 2 || i_dim[0] != 1) {
      mexErrMsgTxt("each cell {i} in regions must be a row vector: 1x(num_nodes_in_reg(i))\n");
    }

    int numNodes_i = i_dim[1];
    double* reg_i_nodes_ptr = mxGetPr(region_i);
    Nodes nodes;
    nodes.clear();
    for (int n=0; n<numNodes_i; n++) {
      int node = (int)(reg_i_nodes_ptr[n]) - 1;
      nodes.push_back(node);
    }
    Region reg;
    reg.assignNodes(nodes);
    if (!regions.addRegion(reg)) {
      mexErrMsgTxt("one region cannot be a subset of another\n");
    }
  }
}

// *****************************************************************************
// fillRegionLevels
// *****************************************************************************

void fillRegionLevels(const mxArray* mxAllRegs, vector<RegionLevel>& allRegions) {
  int all_regs_nd = mxGetNumberOfDimensions(mxAllRegs);
  const int* all_regs_dim = mxGetDimensions(mxAllRegs);
  if (all_regs_nd != 2 || all_regs_dim[0] != 1 || !mxIsCell(mxAllRegs)) {
    mexErrMsgTxt("given adjacency of regions, the regions themselves should be a cell array in size 1x(number of layers), where each cell is a layer in the region-graph.\n this means that each layer should be a cell array in size 1x(number of regions in layer), where each cell is a region in this layer\n");
  }
  int num_layers = all_regs_dim[1];
  for (int l=0; l<num_layers; l++) {
    mxArray* layer = mxGetCell(mxAllRegs,l);
    RegionLevel* regions = new RegionLevel();
    fillRegions(layer,*regions);
    allRegions.push_back(*regions);
    regions->clear();
    delete regions;
    regions = 0;
  }
}

// *****************************************************************************
// fillInitialAssignment
// *****************************************************************************

void fillInitialAssignment(const mxArray* mxAssign, int* startX, int num_nodes) {

  // get startX
  int startX_nd = mxGetNumberOfDimensions(mxAssign);
  const int* startX_dim = mxGetDimensions(mxAssign);
  if (startX_nd != 2 ||
      startX_dim[0] != 1 || startX_dim[1] != num_nodes) {
    mexErrMsgTxt("startX  must be a row vector: 1x(num_nodes)");
  }

  double* startXPtr = mxGetPr(mxAssign);
  // fill startX
  for (int i=0; i<num_nodes; i++) {
    startX[i] = (int)startXPtr[i] - 1;
  }
}

// *****************************************************************************
// fillRhoMat
// *****************************************************************************

void fillRhoMat(const mxArray* mxRho, const MRF* mrf, double** rho) {
  
  int const& N = mrf->N;
  
  int rho_nd = mxGetNumberOfDimensions(mxRho);
  const int* rho_dim = mxGetDimensions(mxRho);
  if (rho_nd != 2 ||
      (rho_dim[0]*rho_dim[1] != N) ||
      !mxIsCell(mxRho)) {
    mexErrMsgTxt("rho must be a cell array in size: 1x(num_nodes)\n");
  }

  // fill rhoMat
  for (int i=0; i<N; i++) {
    mxArray* rho_i = mxGetCell(mxRho, i);
    int i_nd = mxGetNumberOfDimensions(rho_i);
    const int* i_dim = mxGetDimensions(rho_i);
    if (i_nd != 2 ||
	(i_dim[0]*i_dim[1] != mrf->neighbNum(i))) {
      mexErrMsgTxt("each cell {i} in rho must be a vector in length (neighb_num(i))\n");
    }

    double* rho_i_ptr = mxGetPr(rho_i);
    for (int n=0; n<mrf->neighbNum(i); n++) {
      rho[i][n] = rho_i_ptr[n];
    }
  }
}

// *****************************************************************************
// fillInt
// *****************************************************************************

void fillInt(const mxArray* mxIntArray, int* intArray, int arrayLen) {

  int arr_nd = mxGetNumberOfDimensions(mxIntArray);
  const int* arr_dim = mxGetDimensions(mxIntArray);
  if (arr_nd != 2 ||
      arr_dim[0] != 1 || arr_dim[1] != arrayLen) {
    mexErrMsgTxt("intArray must be a row vector: 1x(arrayLen)");
  }

  double* intArrPtr = mxGetPr(mxIntArray);
  for (int i=0; i<arrayLen; i++) {
    intArray[i] = (int)intArrPtr[i];
  }
}

// *****************************************************************************
// fillDouble
// *****************************************************************************

void fillDouble(const mxArray* mxDoubleArray, double* doubleArray, int arrayLen) {

  int arr_nd = mxGetNumberOfDimensions(mxDoubleArray);
  const int* arr_dim = mxGetDimensions(mxDoubleArray);
  if (arr_nd != 2 ||
      arr_dim[0] != 1 || arr_dim[1] != arrayLen) {
    mexErrMsgTxt("doubleArray must be a row vector: 1x(arrayLen)");
  }

  double* doubleArrPtr = mxGetPr(mxDoubleArray);
  for (int i=0; i<arrayLen; i++) {
    doubleArray[i] = doubleArrPtr[i];
  }
}

