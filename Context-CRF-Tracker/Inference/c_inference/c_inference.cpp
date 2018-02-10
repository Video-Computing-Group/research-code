#include "mex.h"
#include "fillMethods.h"
#include "LogLoopy.h"
#include "LogLoopySTime.h"
#include "LogPairsGBP.h"
#include "LogGBP.h"
#include "Gibbs.h"
#include "Wolff.h"
#include "SwendsenWang.h"
#include "Metropolis.h"
#include "LogMeanField.h"
#include "GBPPreProcessor.h"
#include <iostream>

// *****************************************************************************
// enumerators
// *****************************************************************************

enum algorithmType {AT_LOOPY,AT_GBP,AT_GIBBS,AT_WOLFF,AT_SWENDSEN_WANG,
		    AT_METROPOLIS,AT_MEAN_FIELD};

// *****************************************************************************
// mexFunction
// *****************************************************************************

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  // **************************************************************************
  // variables declaration
  // **************************************************************************

  // variables for both Loopy and GBP
  double*** initMsg = 0;
  bool trw = false;
  bool full = true;
  double* countingNode = 0; // relevant for loopy, or for gbp if trw = true
  
  
  // variables for Loopy
  Strategy strategy;
  SumOrMax sumOrMax;
  double gbp_alpha;
  double** rho = 0; // relevant if trw = true
  bool saveTime = false;

  // variables for GBP
  int*** assignInd = 0;
  double* bethe = 0;
  GBPPreProcessor* processor = 0;
  MRF* reg_mrf = 0;
  RegionLevel* regions = 0;
  vector<RegionLevel>* allRegions = 0;
  Potentials* bigRegsPot = 0;
  bool allLevels = false;
  bool regBeliefs = false;
  int num_regs = 0;

  // variables for Monte-Carlo
  int burningTime, samplingInterval, num_samples;
  int* startX = 0;

  // for Loopy, Mean-Field
  bool logspace = false;
  bool logBels = false;

  // for Loopy,GBP,Mean-Field
  int maxIter;
  double threshold;

  
  // **************************************************************************
  // reading input arguments
  // **************************************************************************

  // check number of input arguments.
  // arguments should be:
  //
  // adjMat - 1xN cell array, each cell {i} is a row vector with the indices of
  //          i's neighbours
  //
  // lambda - there are 2 forms for lambda:
  //          1. in general MRF algorithms (loopy, gbp, gibbs, mean-field) :
  //             lambda should be a cell array of 1xN, each cell {i} is a cell
  //             array of 1xneighbNum(i). each cell {i}{n} is a VixVj matrix,
  //             where j is the n-th neighbour of i
  //          2. in PottsMRF alhorithms (monte-carlo algorithms which are planned
  //             for Potts model, i.e. metropolis and the cluster algorithms
  //             wolff and swendsen-wang) :
  //             here lambda should be 1xN cell array, each cell {i} is a row
  //             vector with the strength of interaction of i with each of its
  //             neighbours
  // note: Psi{i,j} = exp( [lambda(i,j), 0; 0, lambda(i,j)] )
  //
  // local - cell array of Nx1, each cell {i} is a row vector of length Vi
  //
  // algorithm - integer representing the inference algorithm to use, see the
  //             enumerator algorithmType at the top of this page
  //
  // temperature - double scalar, the temperature of the system
  //
  // model - integeger representing the model, see the enumerator in "definitions.h"
  //
  // trw - use Tree-Reweighted
  //
  // for other parameters required for each algorithm see the header of "inference.m"
  //
  //
  // note: N = number of nodes, V = number of possible values
  
  if (nrhs < 10 || nrhs > 20) {
    mexErrMsgTxt("Incorrect number of inputs.");
  }

  // get algorithm-type
  algorithmType algo_type = (algorithmType)((int)(mxGetScalar(prhs[3])));
  Model model = (Model)((int)(mxGetScalar(prhs[5])));
  bool potts_model = ((model==POTTS) ||
		      (algo_type==AT_WOLFF) ||
		      (algo_type==AT_SWENDSEN_WANG));
  bool monte_carlo = ((algo_type==AT_GIBBS) ||
		      (algo_type==AT_WOLFF) ||
		      (algo_type==AT_SWENDSEN_WANG) ||
		      (algo_type==AT_METROPOLIS));
  // check number of output arguments
  if ((nlhs > 6) || ((nlhs > 2) && (algo_type != AT_GBP) && (algo_type != AT_LOOPY))) {
    mexErrMsgTxt("Too many output arguments.");
  }

  // get number of nodes and adjMat
  vector<Nodes>* adjMat = new vector<Nodes>();
  fillAdjMat(prhs[0],*adjMat);
  int num_nodes = adjMat->size();

  // define the MRF
  MRF* mrf = 0;
  if (potts_model) {
    mrf = new PottsMRF(*adjMat);
  }
  else {
    mrf = new MRF(*adjMat);
  }
  
  // get local potentials
  fillLocalMat(prhs[2],mrf);
  
  // get pairwise potentials
  if (potts_model) {
    fillLambdaMat(prhs[1],(PottsMRF*)mrf);
  }
  else {
    fillPsiMat(prhs[1],mrf);
  }

  // For monte-carlo algorithms (gibbs, wolff, swendsen-wang), get
  // the initial state and the sampling parameters
  if (monte_carlo) {
    if (nrhs != 10) {
      mexErrMsgTxt("incorrect number of inputs");
    }

    startX = new int[num_nodes];
    fillInitialAssignment(prhs[6], startX, num_nodes);

    // get burningTime, samplingInterval, num_samples
    burningTime = (int)(mxGetScalar(prhs[7]));
    samplingInterval = (int)(mxGetScalar(prhs[8]));
    num_samples = (int)(mxGetScalar(prhs[9]));
    
  }
  else {
    // for all non-monte-carlo-algorithms (Mean-Field, BP & GBP)
    maxIter = (int)(mxGetScalar(prhs[6]));
    threshold = mxGetScalar(prhs[7]);
    // get log-space flag
    logspace = ((int)(mxGetScalar(prhs[8]))) > 0;
    mrf->logspace = logspace;
    logBels = ((int)(mxGetScalar(prhs[9]))) > 0;

    if (algo_type == AT_LOOPY) {

      // for loopy belief propagation:
      // get sum-or-max-flag and strategy

      if (nrhs != 17) {
	mexErrMsgTxt("incorrect number of inputs");
      }
	
      sumOrMax = (SumOrMax)((int)(mxGetScalar(prhs[10])));
      strategy = (Strategy)((int)(mxGetScalar(prhs[11])));
      trw = ((int)(mxGetScalar(prhs[12]))) > 0;
      if (trw) {
	rho = new double*[num_nodes];
	for (int i=0; i<num_nodes; i++) {
	  rho[i] = new double[mrf->neighbNum(i)];
	}
	fillRhoMat(prhs[13],mrf,rho);
      }
	
      // get save-time flag
      saveTime = ((int)(mxGetScalar(prhs[15]))) > 0;

      // get initial messages, if given
      int initM_nd = mxGetNumberOfDimensions(prhs[14]);
      const int* initM_dim = mxGetDimensions(prhs[14]);
      if ((initM_nd == 2) && (initM_dim[0] == 1) &&
	  (initM_dim[1] == num_nodes) && mxIsCell(prhs[14])) {
	initMsg = new double**[num_nodes];
	for (int i=0; i<num_nodes; i++) {
	  mxArray* initMsg_i = mxGetCell(prhs[14],i);
	  int Ni = mrf->neighbNum(i);
	  int len = (saveTime ? num_nodes : Ni);
	  initMsg[i] = new double*[len];
	  if (saveTime) {
	    for (int j=0; j<num_nodes; j++) {
	      initMsg[i][j] = 0;
	    }
	  }
	  for (int n=0; n<Ni; n++) {
	    int j = mrf->adjMat[i][n];
	    int nei = (saveTime ? j : n);
	    initMsg[i][nei] = new double[mrf->V[j]];
	    mxArray* initMsg_ij = mxGetCell(initMsg_i,nei);
	    fillDouble(initMsg_ij,initMsg[i][nei],mrf->V[j]);
	  }
	}
      }

      int count_nd = mxGetNumberOfDimensions(prhs[16]);
      const int* count_dim = mxGetDimensions(prhs[16]);
      if ((count_nd==2) && (count_dim[0]*count_dim[1]==num_nodes)) {
	countingNode = new double[num_nodes];
	fillDouble(prhs[16], countingNode, num_nodes);

	// incorporate local potentials into pairwise
	for (int i=0; i<num_nodes; i++) {
	  if (mrf->neighbNum(i)>0) {
	    int j = mrf->adjMat[i][0];
	    if (i<j) {
	      for (int xi=0; xi<mrf->V[i]; xi++) {
		for (int xj=0; xj<mrf->V[j]; xj++) {
		  mrf->lambdaMat[i][0][xi][xj] *= mrf->localMat[i][xi];
		}
		mrf->localMat[i][xi] = 1.0;
	      }
	    }
	    else {
	      int n = 0;
	      while (mrf->adjMat[j][n] != i) {
		n++;
	      }
	      for (int xi=0; xi<mrf->V[i]; xi++) {
		for (int xj=0; xj<mrf->V[j]; xj++) {
		  mrf->lambdaMat[j][n][xj][xi] *= mrf->localMat[i][xi];
		}
		mrf->localMat[i][xi] = 1.0;
	      }	      
	    }
	  }
	  else {
	    mexPrintf("warning: the graph is not connected\n");
	  }
	}
      }
    }

    if (algo_type == AT_GBP) {

      // for generalized belief propagation:
      // get regions, regions-adj (if given), sum-or-max-flag
      // and alpha

      if (nrhs != 20) {
	mexErrMsgTxt("incorrect number of inputs");
      }

      allLevels = (int)(mxGetScalar(prhs[11])) > 0;
      if (allLevels) {
	allRegions = new vector<RegionLevel>();
	allRegions->clear();
	fillRegionLevels(prhs[10],*allRegions);
      }
      else {
	regions = new RegionLevel();
	fillRegions(prhs[10],*regions);
      }

      sumOrMax = (SumOrMax)((int)(mxGetScalar(prhs[12])));
      gbp_alpha = mxGetScalar(prhs[13]);

      trw = ((int)(mxGetScalar(prhs[14]))) > 0;
      if (trw) {
	countingNode = new double[num_nodes];
	fillDouble(prhs[15], countingNode, num_nodes);
      }
      full = ((int)(mxGetScalar(prhs[16]))) > 0;

      // get initial messages, if given
      int initM_nd = mxGetNumberOfDimensions(prhs[17]);
      const int* initM_dim = mxGetDimensions(prhs[17]);
      if ((initM_nd == 2) && mxIsCell(prhs[17])) {
	num_regs = initM_dim[0] * initM_dim[1];
	initMsg = new double**[num_regs];
	for (int i=0; i<num_regs; i++) {

	  mxArray* initMsg_i = mxGetCell(prhs[17],i);

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

      // get potentials for the big regions, if given
      int regPot_nd = mxGetNumberOfDimensions(prhs[18]);
      const int* regPot_dim = mxGetDimensions(prhs[18]);
      if ((regPot_nd == 2) && mxIsCell(prhs[18])) {
	num_regs = regPot_dim[0] * regPot_dim[1];
	bigRegsPot = new Potentials[num_regs];
	for (int i=0; i<num_regs; i++) {

	  mxArray* regPot_i = mxGetCell(prhs[18],i);

	  int regPot_i_nd = mxGetNumberOfDimensions(regPot_i);
	  const int* regPot_i_dim = mxGetDimensions(regPot_i);
	  if (regPot_i_nd != 2) {
	    mexErrMsgTxt("each cell {i} in big-regions' potentials for GBP should be a vector in length of number of possible states for the region i\n");
	  }
	  int Vi = regPot_i_dim[0] * regPot_i_dim[1];
	  bigRegsPot[i] = new Potential[Vi];
	  fillDouble(regPot_i,bigRegsPot[i],Vi);
	}	
      }

      // if true - get the region beliefs (instead of the single beliefs)
      regBeliefs = ((int)(mxGetScalar(prhs[19]))) > 0;

    }
  }

  // get tepmerature
  double temperature = mxGetScalar(prhs[4]);
  mrf->setTemperature(temperature);
  
  // **************************************************************************
  // create the algorithm
  // **************************************************************************

  InferenceAlgorithm* algorithm = 0;
  switch (algo_type) {

    case AT_LOOPY:

      if (saveTime) {
	if (logspace) {
	  algorithm = new LogLoopySTime(mrf,sumOrMax,strategy,maxIter,rho,initMsg,logBels,threshold);
	}
	else {
	  algorithm = new LoopySTime(mrf,sumOrMax,strategy,maxIter,rho,initMsg,threshold);
	}
      }
      else {
	if (logspace) {
	  if (countingNode != 0) {
	    algorithm = new LogPairsGBP(mrf,sumOrMax,strategy,maxIter,countingNode,initMsg,logBels,threshold);
	  }
	  else {
	    algorithm = new LogLoopy(mrf,sumOrMax,strategy,maxIter,rho,initMsg,logBels,threshold);
	  }
	}
	else {
	  if (countingNode != 0) {
	    algorithm = new PairsGBP(mrf,sumOrMax,strategy,maxIter,countingNode,initMsg,threshold);
	  }
	  else {
	    algorithm = new Loopy(mrf,sumOrMax,strategy,maxIter,rho,initMsg,threshold);
	  }
	}
      }
      break;

    case AT_GBP:
      if (allLevels) {
	processor = new GBPPreProcessor(allRegions, mrf, trw, full, countingNode, bigRegsPot);
      }
      else {
	processor = new GBPPreProcessor(*regions, mrf, trw, full, countingNode, bigRegsPot);
	regions->clear();
	delete regions;
	regions = 0;
      }

      reg_mrf = processor->getRegionMRF();
      assignInd = processor->getAssignTable();
      bethe = processor->getBethe();

      if (logspace) {
	algorithm = new LogGBP(reg_mrf,assignInd,bethe,sumOrMax,gbp_alpha,maxIter,initMsg,logBels,threshold);
      }
      else {
	algorithm = new GBP(reg_mrf,assignInd,bethe,sumOrMax,gbp_alpha,maxIter,initMsg,threshold);
      }
      
      break;

    case AT_GIBBS:

      algorithm = new Gibbs(mrf,startX,burningTime,samplingInterval,num_samples);

      delete[] startX;
      startX = 0;
      
      break;

    case AT_WOLFF:

      algorithm = new Wolff((PottsMRF*)mrf,startX,burningTime,samplingInterval,num_samples);

      delete[] startX;
      startX = 0;
      
      break;

    case AT_SWENDSEN_WANG:

      algorithm = new SwendsenWang((PottsMRF*)mrf,startX,burningTime,samplingInterval,num_samples);

      delete[] startX;
      startX = 0;
      
      break;

    case AT_METROPOLIS:

      algorithm = new Metropolis(mrf,startX,burningTime,samplingInterval,num_samples);

      delete[] startX;
      startX = 0;
      
      break;

    case AT_MEAN_FIELD:

      if (logspace) {
	algorithm = new LogMeanField(mrf,maxIter,logBels,threshold);
      }
      else {
	algorithm = new MeanField(mrf,maxIter,threshold);
      }

      break;
      
    default:

      mexErrMsgTxt("invalid algorithm type. possible values are: 0-loopy, 1-gbp, 2-gibbs, 3-wolff, 4-swendswen-wang, 5-metropolis, 6-mean-field");
      break;
  }
  

  // **************************************************************************
  // make inference
  // **************************************************************************
  int converged;
  double** beliefs = algorithm->inference(&converged);

  double** singleBeliefs = 0;
  double**** pairBeliefs = 0;
  

  switch (algo_type) {
    
    case AT_LOOPY:
      
      if (nlhs > 2) {
	if (countingNode != 0) {
	  pairBeliefs = ((PairsGBP*)algorithm)->calcPairBeliefs();
	}
	else {
	  pairBeliefs = ((Loopy*)algorithm)->calcPairBeliefs();
	}
      }
      
      break;
      
    case AT_GBP:

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
      if (!regBeliefs || nlhs > 4) {
	
	singleBeliefs = new double*[num_nodes];
	for (int i=0; i<num_nodes; i++) {
	  singleBeliefs[i] = new double[mrf->V[i]];      
	}
	processor->extractSingle(beliefs,singleBeliefs,sumOrMax);
      
	if ((!regBeliefs && nlhs > 2) || (regBeliefs && nlhs > 5)) {
	  pairBeliefs = new double***[num_nodes];
	  for (int i=0; i<num_nodes; i++) {
	    pairBeliefs[i] = new double**[mrf->neighbNum(i)];
	    for (int n=0; n<mrf->neighbNum(i); n++) {
	      pairBeliefs[i][n] = 0;
	      int j = mrf->adjMat[i][n];
	      if (i<j) {
		pairBeliefs[i][n] = new double*[mrf->V[i]];
		for (int xi=0; xi<mrf->V[i]; xi++) {
		  pairBeliefs[i][n][xi] = new double[mrf->V[j]];
		}
	      }
	    }
	  }
	  processor->extractPairs(beliefs,pairBeliefs,sumOrMax);

	}
	if (!regBeliefs) {
	  beliefs = singleBeliefs;
	}
      }
      break;
      
    default:
      
      break;
  }
  
  // **************************************************************************
  // assign results to output argument (if given)
  // **************************************************************************

  if (regBeliefs) {
    int num_regs = reg_mrf->N;
    int regs_dims[2] = {1,num_regs};

    Region** all_regions = processor->getAllRegions();

    // assign: 1. regions 2. regions' adj-matrix 3. region beliefs 4. convergence flag
    plhs[0] = mxCreateCellArray(2,regs_dims);
    plhs[1] = mxCreateCellArray(2,regs_dims);
    plhs[2] = mxCreateCellArray(2,regs_dims);
    plhs[3] = mxCreateDoubleScalar(converged);

    for (int i=0; i<num_regs; i++) {

      // regions
      Region* region = all_regions[i];
      int reg_i_size = (int)(region->size());
      int reg_i_dims[2] = {1,reg_i_size};
      mxArray* reg_i = mxCreateNumericArray(2,reg_i_dims,mxDOUBLE_CLASS, mxREAL);
      double* reg_i_ptr = mxGetPr(reg_i);
      for (int n=0; n<reg_i_size; n++) {
	reg_i_ptr[n] = (double)((*region)[n] + 1);
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

      // region beliefs
      int val_dims[2] = {reg_mrf->V[i],1};
      mxArray* bel_i = mxCreateNumericArray(2,val_dims,mxDOUBLE_CLASS, mxREAL);
      double* resBelPtr = mxGetPr(bel_i);
      for (int xi=0; xi<reg_mrf->V[i]; xi++) {
	resBelPtr[xi] = beliefs[i][xi];
      }
      mxSetCell(plhs[2],i,bel_i);

    }

    if (nlhs > 4) {
      int bel_dims[2] = {1,num_nodes};

      // assign: 5. single beliefs 6. pairwise beliefs (if required)
      plhs[4] = mxCreateCellArray(2,bel_dims);
      if (nlhs > 5) {
	plhs[5] = mxCreateCellArray(2,bel_dims);
      }

      for (int i=0; i<num_nodes; i++) {

	// single beliefs
	int val_dims[2] = {mrf->V[i],1};
	mxArray* bel_i = mxCreateNumericArray(2,val_dims,mxDOUBLE_CLASS, mxREAL);
	double* resBelPtr = mxGetPr(bel_i);
	for (int xi=0; xi<mrf->V[i]; xi++) {
	  resBelPtr[xi] = singleBeliefs[i][xi];
	}
	mxSetCell(plhs[4],i,bel_i);

	// pairwise beliefs
	if (nlhs > 5) {
	  int pair_i_dims[2] = {1, mrf->neighbNum(i)};
	  mxArray* pbel_i = mxCreateCellArray(2,pair_i_dims);
	  
	  for (int n=0; n<mrf->neighbNum(i); n++) {
	    int j = mrf->adjMat[i][n];
	    if (i<j) {
	      int pval_dims[2] = {mrf->V[i], mrf->V[j]};
	      mxArray* pbel_ij = mxCreateNumericArray(2,pval_dims,mxDOUBLE_CLASS,mxREAL);
	      
	      double* resPBelPtr = mxGetPr(pbel_ij);
	      for (int xi=0; xi<mrf->V[i]; xi++) {
		for (int xj=0; xj<mrf->V[j]; xj++) {
		  resPBelPtr[xi + xj*mrf->V[i]] = pairBeliefs[i][n][xi][xj];
		}
	      }
	      mxSetCell(pbel_i, n, pbel_ij);	      
	    }
	  }
	  mxSetCell(plhs[5], i, pbel_i);
	}
      }
    }
  }
  else {
    if (nlhs > 0) {
      int bel_dims[2] = {1,num_nodes};
      plhs[0] = mxCreateCellArray(2,bel_dims);
      for (int i=0; i<num_nodes; i++) {
	int val_dims[2] = {mrf->V[i],1};
	mxArray* bel_i = mxCreateNumericArray(2,val_dims,mxDOUBLE_CLASS, mxREAL);
	double* resBelPtr = mxGetPr(bel_i);
	for (int xi=0; xi<mrf->V[i]; xi++) {
	  resBelPtr[xi] = beliefs[i][xi];
	}
	mxSetCell(plhs[0],i,bel_i);
      }
      if (nlhs > 1) {
	plhs[1] = mxCreateDoubleScalar(converged); // For matlab6.5
	//plhs[1] = mxCreateScalarDouble(converged);
	if (nlhs > 2) {
	  int pair_dims[2] = {1,num_nodes};
	  plhs[2] = mxCreateCellArray(2,pair_dims);

	  for (int i=0; i<num_nodes; i++) {
	    int pair_i_dims[2] = {1, mrf->neighbNum(i)};
	    mxArray* bel_i = mxCreateCellArray(2,pair_i_dims);
	  
	    for (int n=0; n<mrf->neighbNum(i); n++) {
	      int j = mrf->adjMat[i][n];
	      if (i<j) {
		int pval_dims[2] = {mrf->V[i], mrf->V[j]};
		mxArray* bel_ij = mxCreateNumericArray(2,pval_dims,mxDOUBLE_CLASS,mxREAL);
	      
		double* resBelPtr = mxGetPr(bel_ij);
		for (int xi=0; xi<mrf->V[i]; xi++) {
		  for (int xj=0; xj<mrf->V[j]; xj++) {
		    resBelPtr[xi + xj*mrf->V[i]] = pairBeliefs[i][n][xi][xj];
		  }
		}
		mxSetCell(bel_i, n, bel_ij);	      
	      }
	    }
	    mxSetCell(plhs[2], i, bel_i);
	  }

	  if (nlhs > 3) {

	    double*** msg = 0;
	    int msg_dims[2];
	  
	    switch (algo_type) {
    
	      case AT_LOOPY:
      
		if (countingNode != 0) {
		  msg = ((PairsGBP*)algorithm)->getMessages();
		}
		else {
		  msg = ((Loopy*)algorithm)->getMessages();
		}
	    
		msg_dims[0] = 1;
		msg_dims[1] = num_nodes;
	      
		plhs[3] = mxCreateCellArray(2,msg_dims);
		for (int i=0; i<num_nodes; i++) {
		  int Ni = mrf->neighbNum(i);
		  int msg_i_dims[2];
		  msg_i_dims[0] = 1;
		  msg_i_dims[1] = (saveTime ? num_nodes : Ni);
		  mxArray* msg_i = mxCreateCellArray(2,msg_i_dims);
		  for (int n=0; n<Ni; n++) {
		    int j = mrf->adjMat[i][n];
		    int nei = (saveTime ? j : n);
		    int msg_ij_dims[2] = {1, mrf->V[j]};
		    mxArray* msg_ij = mxCreateNumericArray(2,msg_ij_dims,mxDOUBLE_CLASS, mxREAL);
		    double* msg_ij_ptr = mxGetPr(msg_ij);
		    for (int xj=0; xj<mrf->V[j]; xj++) {
		      msg_ij_ptr[xj] = msg[i][nei][xj];
		    }
		    mxSetCell(msg_i,nei,msg_ij);
		  }
		  mxSetCell(plhs[3],i,msg_i);
		}
		

		break;
      
	      case AT_GBP:

		if (!trw) {
		  msg = ((GBP*)algorithm)->getMessages();

		  msg_dims[0] = 1;
		  msg_dims[1] = reg_mrf->N;
		  plhs[3] = mxCreateCellArray(2,msg_dims);
		  for (int i=0; i<reg_mrf->N; i++) {
		    int Ni = reg_mrf->neighbNum(i);
		    int msg_i_dims[2] = {1,Ni};
		    mxArray* msg_i = mxCreateCellArray(2,msg_i_dims);
		    for (int n=0; n<Ni; n++) {
		      int j = reg_mrf->adjMat[i][n];
		      int numStates = reg_mrf->V[max(i,j)];
		      int msg_ij_dims[2] = {1, numStates};
		      mxArray* msg_ij = mxCreateNumericArray(2,msg_ij_dims,mxDOUBLE_CLASS, mxREAL);
		      double* msg_ij_ptr = mxGetPr(msg_ij);
		      for (int xs=0; xs<numStates; xs++) {
			msg_ij_ptr[xs] = msg[i][n][xs];
		      }
		      mxSetCell(msg_i,n,msg_ij);
		    }
		    mxSetCell(plhs[3],i,msg_i);
		  }
		}
		break;

	      default:
	      
		break;
	    
	    }
	  
	  }
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
  if (pairBeliefs != 0 && algo_type == AT_GBP) {
    for (int i=0; i<num_nodes; i++) {
      for (int n=0; n<mrf->neighbNum(i); n++) {
	if (pairBeliefs[i][n] != 0) {
	  for (int xi=0; xi<mrf->V[i]; xi++) {
	    delete[] pairBeliefs[i][n][xi];
	  }
	  delete[] pairBeliefs[i][n];
	}
      }
      delete[] pairBeliefs[i];
    }
    delete[] pairBeliefs;
    pairBeliefs = 0;
  }
  if (processor != 0) {
    delete processor;
    processor = 0;
  }
  if (rho != 0) {
    for (int i=0; i<num_nodes; i++) {
      delete[] rho[i];
      rho[i] = 0;
    }
    delete[] rho;
    rho = 0;
  }
  if (countingNode != 0) {
    delete[] countingNode;
    countingNode = 0;
  }
  if (bigRegsPot != 0) {
    for (int i=0; i<num_regs; i++) {
      delete[] bigRegsPot[i];
    }
    delete[] bigRegsPot;
    bigRegsPot = 0;
  }
  
  delete mrf;
  mrf = 0;
  delete adjMat;
  adjMat = 0;
  
}


