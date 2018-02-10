#include "mex.h"
#include "fillMethods.h"
#include "GBP.h"

// *****************************************************************************
// methods declaration
// *****************************************************************************
void extractSingleBeliefs(double** regsBeliefs, double** singleBeliefs,
			  int* extractSingle, RegionLevel& regions,
			  int const* regCard, int const* nodeCard,
			  int num_nodes, SumOrMax smFlag);
void extractPairBeliefs(double** regsBeliefs, double**** pairBeliefs,
			int** extractPairs, RegionLevel& regions,
			int const* regCard, int const* nodeCard,
			int num_nodes, vector<Nodes>& adjMat,
			SumOrMax smFlag);

// *****************************************************************************
// mexFunction
// *****************************************************************************

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  // **************************************************************************
  // variables declaration
  // **************************************************************************

  MRF* reg_mrf = 0;
  RegionLevel* regions = 0;
  int*** assignInd = 0;
  double* bethe = 0;
  int* extractSingle = 0;
  int** extractPairs = 0;
  int* Vn = 0;
  SumOrMax sumOrMax;
  double alpha;
  double*** initMsg = 0;
  int maxIter;

  // **************************************************************************
  // reading input arguments
  // **************************************************************************

  // check number of input arguments.
  // arguments should be:
  // bigRegions: 1xNr cell array, each cell {i} is a row vector with the indices
  //             of nodes in the region i
  // adjMat:     the nodes' adjacencies matrix, an 1xNn cell array,
  //             each cell {i} is a row vector with the indices of i's
  //             neighbours
  // regsAdjMat: the regions' adjacencies matrix, an 1xNr cell array,
  //             each cell {i} is a row vector with the indices of i's
  //             neighbours
  // Vn:         the number of possible states for each node
  // regsLocal:  1xNr cell array, each cell {i} is a row vector of length Vr_i
  // 
  // assignInd:  1xNr cell array, each cell {i} is a cell array of 1xneighbNum(i),
  //             in which each cell {n} is:
  //             1. empty if i>j - where j is the n'th neighbour of i
  //             2. a column vector of length Vr_i with the indices of the
  //                corresponding assignment of j, if i<j
  // bethe:      a row vector of length Nr, with the data related to the double
  //             count of each region
  // extractSingle: a row vector of length Nn, where each element (i) is the
  //                index of the region from which the belief for node i should
  //                be taken
  // extractPairs:  1xNn cell array, each cell {i} is a row vector of length
  //                neighbNum(i). each element (n) in this vector is:
  //                1. -1 if i>j - where j is the n'th neighbour of i
  //                2. the index of the region from which the pairwise beliefs
  //                   of the pair <i,j> should be extracted, if i<j
  // sumOrMax:      0 - sum gbp, 1 - max gbp
  // alpha:         averaging parameter between new messages - weight alpha -
  //                and old messages - weight 1-alpha
  // initMsg:       OPTIONAL, defines the initial messages between the regions
  //                1xNr cell array, each cell {i} is a cell array of 1xneighbNum(i),
  //                in which each cell {n} is a vector of 1xnumStates(max(i,j))
  //                (which means the number of states of the son-region). this vector
  //                holds the initial messages from region i to region j.
  //
  // ====
  // note: Nr = number of regions, Nn = number of nodes,
  // ====  Vr = number of possible values for each region,
  //       Vn = number of possible values for each node
  
  if (nrhs != 13) {
    mexErrMsgTxt("Incorrect number of inputs.");
  }

  // check number of output arguments
  if (nlhs > 4) {
    mexErrMsgTxt("Too many output arguments.");
  }

  // get number of nodes and adjMat
  vector<Nodes>* adjMat = new vector<Nodes>();
  fillAdjMat(prhs[1],*adjMat);
  int num_nodes = adjMat->size();

  // get number of regions and regsAdjMat
  vector<Nodes>* regsAdjMat = new vector<Nodes>();
  fillAdjMat(prhs[2],*regsAdjMat);
  int num_regs = regsAdjMat->size();

  // get number of states for each node
  Vn = new int[num_nodes];
  fillInt(prhs[3],Vn,num_nodes);  
  
  // define the MRF
  reg_mrf = new MRF(*regsAdjMat);
  // get local potentials
  fillLocalMat(prhs[4],reg_mrf);
  
  // get regions
  regions = new RegionLevel();
  fillRegions(prhs[0],*regions);

  // get assignInd
  assignInd = new int**[num_regs];
  for (int i=0; i<num_regs; i++) {
    mxArray* ass_i = mxGetCell(prhs[5],i);
    int Ni = reg_mrf->neighbNum(i);
    assignInd[i] = new int*[Ni];
    for (int n=0; n<Ni; n++) {
      assignInd[i][n] = 0;
      int j = reg_mrf->adjMat[i][n];
      if (i<j) {
	mxArray* ass_ij = mxGetCell(ass_i,n);
	int Vi = reg_mrf->V[i];
	assignInd[i][n] = new int[Vi];
	fillInt(ass_ij,assignInd[i][n],Vi);
      }
    }
  }

  // get bethe
  bethe = new double[num_regs];
  fillDouble(prhs[6],bethe,num_regs);
  
  // get extractSingle information
  extractSingle = new int[num_nodes];
  fillInt(prhs[7],extractSingle,num_nodes);
  
  extractPairs = new int*[num_nodes];
  for (int i=0; i<num_nodes; i++) {
    mxArray* exPairs_i = mxGetCell(prhs[8],i);
    int Ni = (*adjMat)[i].size();
    extractPairs[i] = new int[Ni];
    fillInt(exPairs_i,extractPairs[i],Ni);
  }

  // get sum-or-max flag
  sumOrMax = (SumOrMax)((int)(mxGetScalar(prhs[9])));

  // get alpha
  alpha = mxGetScalar(prhs[10]);

  // get maximum number of iterations
  maxIter = (int)(mxGetScalar(prhs[11]));

  // get initial messages, if given
  int initM_nd = mxGetNumberOfDimensions(prhs[12]);
  const int* initM_dim = mxGetDimensions(prhs[12]);
  if ((initM_nd == 2) && (initM_dim[0] == 1) &&
      (initM_dim[1] == num_regs) && mxIsCell(prhs[12])) {
    initMsg = new double**[num_regs];
    for (int i=0; i<num_regs; i++) {
      mxArray* initMsg_i = mxGetCell(prhs[12],i);
      int Ni = reg_mrf->neighbNum(i);
      initMsg[i] = new double*[Ni];
      for (int n=0; n<Ni; n++) {
	int j = reg_mrf->adjMat[i][n];
	int numStates = reg_mrf->V[max(i,j)];
	initMsg[i][n] = new double[numStates];
	mxArray* initMsg_ij = mxGetCell(initMsg_i,n);
	fillDouble(initMsg_ij,initMsg[i][n],numStates);
      }
    }
  }

  
  // **************************************************************************
  // create the algorithm and make inference
  // **************************************************************************

  GBP* gbp = new GBP(reg_mrf,assignInd,bethe,sumOrMax,alpha,maxIter,initMsg);
  
  int converged;
  double** beliefs = gbp->inference(&converged);
  
  double** singleBeliefs = 0;
  double**** pairBeliefs = 0;

  if (nlhs > 0) {
    singleBeliefs = new double*[num_nodes];
    for (int i=0; i<num_nodes; i++) {
      singleBeliefs[i] = new double[Vn[i]];      
    }
    extractSingleBeliefs(beliefs,singleBeliefs,extractSingle,*regions,
			 reg_mrf->V,Vn,num_nodes,sumOrMax);

    if (nlhs > 3) {
      pairBeliefs = new double***[num_nodes];
      for (int i=0; i<num_nodes; i++) {
	int Ni = (*adjMat)[i].size();
	pairBeliefs[i] = new double**[Ni];
	for (int n=0; n<Ni; n++) {
	  pairBeliefs[i][n] = 0;
	  int j = (*adjMat)[i][n];
	  if (i<j) {
	    pairBeliefs[i][n] = new double*[Vn[i]];
	    for (int xi=0; xi<Vn[i]; xi++) {
	      pairBeliefs[i][n][xi] = new double[Vn[j]];
	    }
	  }
	}
      }
      extractPairBeliefs(beliefs,pairBeliefs,extractPairs,*regions,
			 reg_mrf->V,Vn,num_nodes,*adjMat,sumOrMax);
    }
  }
  
  // **************************************************************************
  // assign results to output argument (if given)
  // **************************************************************************
  // output arguments:
  // 1. regions' messages
  // 2. single nodes' beliefs
  // 3. cflag - number of iterations untill convergance, -1 if didn't converge
  // 4. pairwise beliefs
  //

  if (nlhs > 0) {

    // 1. regions' messages
    
    double*** msg = gbp->getMessages();
    int msg_dims[2] = {1,num_regs};
    plhs[0] = mxCreateCellArray(2,msg_dims);
    for (int i=0; i<num_regs; i++) {
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
      mxSetCell(plhs[0],i,msg_i);
    }

    if (nlhs > 1) {

      // 2. single nodes' beliefs
      
      int bel_dims[2] = {1,num_nodes};
      plhs[1] = mxCreateCellArray(2,bel_dims);
      for (int i=0; i<num_nodes; i++) {
	int val_dims[2] = {Vn[i],1};
	mxArray* bel_i = mxCreateNumericArray(2,val_dims,mxDOUBLE_CLASS, mxREAL);
	double* resBelPtr = mxGetPr(bel_i);
	for (int xi=0; xi<Vn[i]; xi++) {
	  resBelPtr[xi] = singleBeliefs[i][xi];
	}
	mxSetCell(plhs[1],i,bel_i);
      }
      
      if (nlhs > 2) {

	// 3. convergence flag
	
	plhs[2] = mxCreateDoubleScalar(converged); // For matlab6.5
	//plhs[2] = mxCreateScalarDouble(converged);
	
	if (nlhs > 3) {

	  // 4. pairwise beliefs
	  
	  int pair_dims[2] = {1,num_nodes};
	  plhs[3] = mxCreateCellArray(2,pair_dims);

	  for (int i=0; i<num_nodes; i++) {
	  
	    int Ni = (*adjMat)[i].size();
	    int pair_i_dims[2] = {1, Ni};
	    mxArray* bel_i = mxCreateCellArray(2,pair_i_dims);

	    for (int n=0; n<Ni; n++) {
	      int j = (*adjMat)[i][n];
	      if (i<j) {
		int pval_dims[2] = {Vn[i], Vn[j]};
		mxArray* bel_ij = mxCreateNumericArray(2,pval_dims,mxDOUBLE_CLASS,mxREAL);
	      
		double* resBelPtr = mxGetPr(bel_ij);
		for (int xi=0; xi<Vn[i]; xi++) {
		  for (int xj=0; xj<Vn[j]; xj++) {
		    resBelPtr[xi + xj*Vn[i]] = pairBeliefs[i][n][xi][xj];
		  }
		}
		mxSetCell(bel_i, n, bel_ij);	      
	      }
	    }
	    mxSetCell(plhs[3], i, bel_i);
	  }
	}
      }
    }
  }

  // **************************************************************************
  // free memory
  // **************************************************************************

  delete gbp;
  gbp = 0;

  if (singleBeliefs != 0) {
    for (int i=0; i<num_nodes; i++) {
      delete[] singleBeliefs[i];
    }
    delete[] singleBeliefs;    
    singleBeliefs = 0;
  }
  if (pairBeliefs != 0) {
    for (int i=0; i<num_nodes; i++) {
      int Ni = (*adjMat)[i].size();
      for (int n=0; n<Ni; n++) {
	if (pairBeliefs[i][n] != 0) {
	  for (int xi=0; xi<Vn[i]; xi++) {
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

  if (assignInd != 0) {
    for (int i=0; i<num_regs; i++) {
      int Ni = (*regsAdjMat)[i].size();
      for (int n=0; n<Ni; n++) {
	if (assignInd[i][n] != 0) {
	  delete[] assignInd[i][n];
	}
      }
      delete[] assignInd[i];
    }
    delete[] assignInd;
    assignInd = 0;
  }

  delete reg_mrf;
  reg_mrf = 0;

  delete adjMat;
  adjMat = 0;

  delete regsAdjMat;
  adjMat = 0;

  delete[] Vn;
  Vn = 0;

  delete regions;
  regions = 0;

  delete[] bethe;
  bethe = 0;

  delete[] extractSingle;
  extractSingle = 0;

  for (int i=0; i<num_nodes; i++) {
    delete[] extractPairs[i];
  }
  delete[] extractPairs;
  extractPairs = 0;
}


// *****************************************************************************
// methods implementation
// *****************************************************************************
void extractSingleBeliefs(double** regsBeliefs, double** singleBeliefs,
			  int* extractSingle, RegionLevel& regions,
			  int const* regCard, int const* nodeCard,
			  int num_nodes, SumOrMax smFlag) {
  
  Assignment assignment;
  assignment.clear();
  assignment.resize(num_nodes);


  for (int n=0; n<num_nodes; n++) {
    int reg = extractSingle[n];
    Region& region = regions[reg];
    switch (smFlag) {
      case SUM:
	for (int xn=0; xn<nodeCard[n]; xn++) {
	  singleBeliefs[n][xn] = 0.0;
	}
	break;
      case MAX:
	for (int xn=0; xn<nodeCard[n]; xn++) {
	  singleBeliefs[n][xn] = -1.0;
	}
	break;
      default:
	break;
    }

    for (int reg_assignInd=0; reg_assignInd<regCard[reg]; reg_assignInd++) {
      region.indToAssign(reg_assignInd, assignment, nodeCard);
      int xn = assignment[n];
      switch (smFlag) {
	case SUM:
	  singleBeliefs[n][xn] += regsBeliefs[reg][reg_assignInd];
	  break;
	case MAX:
	  if (regsBeliefs[reg][reg_assignInd] > singleBeliefs[n][xn]) {
	    singleBeliefs[n][xn] = regsBeliefs[reg][reg_assignInd];
	  }
	  break;
	default:
	  break;
      }
    }

    // normalize for max-gbp
    if (smFlag == MAX) {
      double sum_bel_n = 0.0;
      for (int xn=0; xn<nodeCard[n]; xn++) {
	sum_bel_n += singleBeliefs[n][xn];      
      }
      for (int xn=0; xn<nodeCard[n]; xn++) {
	if (sum_bel_n > 0.0) {
	  singleBeliefs[n][xn] /= sum_bel_n;
	}
      }    
    }
  }
  
}

void extractPairBeliefs(double** regsBeliefs, double**** pairBeliefs,
			int** extractPairs, RegionLevel& regions,
			int const* regCard, int const* nodeCard,
			int num_nodes, vector<Nodes>& adjMat,
			SumOrMax smFlag) {
  
  Assignment assignment;
  assignment.clear();
  assignment.resize(num_nodes);


  for (int i=0; i<num_nodes; i++) {
    int Ni = adjMat[i].size();
    for (int n=0; n<Ni; n++) {
      int j = adjMat[i][n];
      if (i<j) {
	int reg = extractPairs[i][n];
	Region& region = regions[reg];
	
	for (int xi=0; xi<nodeCard[i]; xi++) {
	  for (int xj=0; xj<nodeCard[j]; xj++) {
	    
	    switch (smFlag) {
	      case SUM:
		pairBeliefs[i][n][xi][xj] = 0.0;
		break;
	      case MAX:
		pairBeliefs[i][n][xi][xj] = -1.0;
		break;
	      default:
		break;
	    }
	  }
	}
	for (int reg_assignInd=0; reg_assignInd<regCard[reg]; reg_assignInd++) {
	  region.indToAssign(reg_assignInd, assignment, nodeCard);
	  int xi = assignment[i];
	  int xj = assignment[j];
	  switch (smFlag) {
	    case SUM:
	      pairBeliefs[i][n][xi][xj] += regsBeliefs[reg][reg_assignInd];
	      break;
	    case MAX:
	      if (regsBeliefs[reg][reg_assignInd] > pairBeliefs[i][n][xi][xj]) {
		pairBeliefs[i][n][xi][xj] = regsBeliefs[reg][reg_assignInd];
	      }
	      break;
	    default:
	      break;
	  }
	}
	// normalize for max-gbp
	if (smFlag == MAX) {
	  double sum_bel_ij = 0.0;
	  for (int xi=0; xi<nodeCard[i]; xi++) {
	    for (int xj=0; xj<nodeCard[j]; xj++) {
	      sum_bel_ij += pairBeliefs[i][n][xi][xj];
	    }
	  }
	  for (int xi=0; xi<nodeCard[i]; xi++) {
	    for (int xj=0; xj<nodeCard[j]; xj++) {
	      pairBeliefs[i][n][xi][xj] /= sum_bel_ij;
	    }
	  }
	}
      }
    }
  }  
}
