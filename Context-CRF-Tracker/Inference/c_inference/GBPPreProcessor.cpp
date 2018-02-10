#include "GBPPreProcessor.h"
#include "MathFunctions.h"
#include "mex.h"
#include <iostream>

using namespace std;

GBPPreProcessor::GBPPreProcessor(RegionLevel const& bigRegions, MRF const* mrf,
				 bool trw, bool full, double* countingNode, Potentials* bigRegsPot) {

  pp_mrf = mrf;
  init();

  pp_numBigRegs = bigRegions.size();
  pp_numAllRegs = bigRegions.size();
  pp_numRegsLevels = (pp_numAllRegs > 0 ? 1 : 0);  
    
  pp_allRegions = new vector<RegionLevel>();
  pp_allRegions->clear();
  pp_allRegions->push_back(bigRegions);

  pp_numRegsInLevel = new vector<int>();
  pp_numRegsInLevel->clear();
  pp_numRegsInLevel->push_back(bigRegions.size());

  pp_numRegsBeforeLevel = new vector<int>();
  pp_numRegsBeforeLevel->clear();
  pp_numRegsBeforeLevel->push_back(0);
  
  pp_regions_mrf = 0;
  pp_adjMat = new vector<Nodes>();
  pp_countingNode = 0;

  if (full) {
    buildFullRegionGraph();
    if (trw) {
      createConvexRegionsAdj();
    }
    else {
      createRegionsAdj();
    }
  }
  else {
    build2LayersRegionGraph();
    
    if (countingNode != 0) {
      pp_countingNode = new double[mrf->N];
      for (int i=0; i<mrf->N; i++) {
	pp_countingNode[i] = countingNode[i];
      }
    }
    create2LayersRegionsAdj();
  }
  calcBethe();
  calcNumOfStates();
  calcAssignTable();
  calcPotentials(bigRegsPot);
  pp_regions_mrf->logspace = pp_mrf->logspace;
  if (pp_mrf->logspace)
    pp_regions_mrf->setTemperature(pp_mrf->getTemperature());

}

GBPPreProcessor::GBPPreProcessor(vector<RegionLevel>* allRegions, MRF const* mrf,
				 bool trw, bool full, double* countingNode, Potentials* bigRegsPot) {
  pp_mrf = mrf;
  
  init();

  pp_allRegions = allRegions;
  pp_numRegsLevels = allRegions->size();
  pp_numBigRegs = (*pp_allRegions)[0].size();
  pp_numAllRegs = pp_numBigRegs;

  pp_numRegsInLevel = new vector<int>();
  pp_numRegsInLevel->clear();
  pp_numRegsInLevel->push_back(pp_numBigRegs);

  pp_numRegsBeforeLevel = new vector<int>();
  pp_numRegsBeforeLevel->clear();
  pp_numRegsBeforeLevel->push_back(0);
  
  for (int level=1; level<pp_numRegsLevels; level++) {
    pp_numRegsBeforeLevel->push_back(pp_numAllRegs);
    pp_numAllRegs += (*pp_allRegions)[level].size();
    pp_numRegsInLevel->push_back((*pp_allRegions)[level].size());
  }

  // initialize the pointers array to the regions
  pp_allRegionsPtr = new Region*[pp_numAllRegs];
  int i = 0;
  for (int level=0; level<pp_numRegsLevels; level++) {
    RegionLevel& regions = (*pp_allRegions)[level];
    for (int r=0; r<(int)(regions.size()); r++) {
      pp_allRegionsPtr[i] = &(regions[r]);
      i++;
    }
  }
  
  pp_regions_mrf = 0;
  pp_adjMat = new vector<Nodes>();
  pp_countingNode = 0;

  if (trw) {    
    if (countingNode != 0) {
      pp_countingNode = new double[mrf->N];
      for (int i=0; i<mrf->N; i++) {
	pp_countingNode[i] = countingNode[i];
      }
    }
    if (full) {
      createConvexRegionsAdj();
    }
    else {
      create2LayersRegionsAdj();
    }
  }
  else {
    createRegionsAdj();
  }
  calcBethe();
  calcNumOfStates();
  calcAssignTable();
  calcPotentials(bigRegsPot);
  pp_regions_mrf->logspace = pp_mrf->logspace;
  if (pp_mrf->logspace)
    pp_regions_mrf->setTemperature(pp_mrf->getTemperature());
}

void GBPPreProcessor::init() {

  pp_parents = 0;
  pp_doubleCount = 0;
  pp_bethe = 0;
  pp_assignTable = 0;
  pp_allRegionsPtr = 0;
  pp_extractSingle = 0;
  pp_extractPairs = 0;
  
}

GBPPreProcessor::~GBPPreProcessor() {

  pp_allRegions->clear();
  delete pp_allRegions;
  pp_allRegions = 0;

  if (pp_numRegsInLevel != 0) {
    pp_numRegsInLevel->clear();
    delete pp_numRegsInLevel;
    pp_numRegsInLevel = 0;
  }

  if (pp_numRegsBeforeLevel != 0) {
    pp_numRegsBeforeLevel->clear();
    delete pp_numRegsBeforeLevel;
    pp_numRegsBeforeLevel = 0;
  }
  
  if (pp_parents != 0) {
    delete[] pp_parents;
    pp_parents = 0;
  }
  if (pp_doubleCount != 0) {
    delete[] pp_doubleCount;
    pp_doubleCount = 0;
  }
  if (pp_bethe != 0) {
    delete[] pp_bethe;
    pp_bethe = 0;
  }

  if (pp_allRegionsPtr != 0) {
    delete[] pp_allRegionsPtr;
    pp_allRegionsPtr = 0;
  }
  
  if (pp_assignTable != 0) {
    for (int i=0; i<pp_numAllRegs; i++) {
      for (int n=0; n<pp_regions_mrf->neighbNum(i); n++) {
	if (pp_assignTable[i][n] != 0) {
	  delete[] pp_assignTable[i][n];
	}
      }
      delete[] pp_assignTable[i];
    }
    delete[] pp_assignTable;
    pp_assignTable = 0;
  }

  if (pp_regions_mrf != 0) {
    delete pp_regions_mrf;
    pp_regions_mrf = 0;
  }
  delete pp_adjMat;
  pp_adjMat = 0;
  
  if (pp_extractSingle != 0) {
    delete[] pp_extractSingle;
    pp_extractSingle = 0;
  }
  if (pp_extractPairs != 0) {
    for (int i=0; i<pp_mrf->N; i++) {
      delete[] pp_extractPairs[i];
    }
    delete[] pp_extractPairs;
    pp_extractPairs = 0;
  }

}

void GBPPreProcessor::buildFullRegionGraph() {

  // create all levels of regions
  
  RegionLevel* currLevel = &((*pp_allRegions)[0]);
  RegionLevel nextLevel;
  nextLevel.clear();
  currLevel->intersections(nextLevel);
  while (!nextLevel.empty()) {

    // add this level to the region graph
    pp_numRegsBeforeLevel->push_back(pp_numAllRegs);
    pp_numAllRegs += nextLevel.size();
    pp_numRegsLevels++;
    pp_allRegions->push_back(nextLevel);
    pp_numRegsInLevel->push_back(nextLevel.size());

    // try to create the next intersections-level
    currLevel = &((*pp_allRegions)[pp_numRegsLevels-1]);
    nextLevel.clear();
    currLevel->intersections(nextLevel);
  }

  // initialize the pointers array to the regions
  pp_allRegionsPtr = new Region*[pp_numAllRegs];
  int i=0;
  for (int level=0; level<pp_numRegsLevels; level++) {
    RegionLevel& regLevel = (*pp_allRegions)[level];
    for (int r=0; r<(*pp_numRegsInLevel)[level]; r++) {
      pp_allRegionsPtr[i] = &(regLevel[r]);
      i++;
    }
  }

}


void GBPPreProcessor::createRegionsAdj() {

  // define neighbours

  pp_adjMat->resize(pp_numAllRegs);
  
  pp_parents = new int[pp_numAllRegs]; // will be cleared after calculating behte
  for (int i=0; i<pp_numAllRegs; i++) {
    pp_parents[i] = 0;
  }

  // and calculate double counts
  
  pp_doubleCount = new double[pp_numAllRegs]; // will be cleared after calculating bethe
  for (int i=0; i<pp_numBigRegs; i++) {
    pp_doubleCount[i] = 1.0;
  }
  for (int i=pp_numBigRegs; i<pp_numAllRegs; i++) {
    pp_doubleCount[i] = 0.0;
  }

  for (int sLevel=1; sLevel<pp_numRegsLevels; sLevel++) {

    RegionLevel& sonLevel = (*pp_allRegions)[sLevel];

    for (int pLevel=0; pLevel<sLevel; pLevel++) {

      RegionLevel& parentLevel = (*pp_allRegions)[pLevel];
      bool directSon = (sLevel==(pLevel+1));

      for (int i=0; i<(*pp_numRegsInLevel)[pLevel]; i++) {

	Region& parentRegion = parentLevel[i];
	int pid = (*pp_numRegsBeforeLevel)[pLevel] + i;

	for (int j=0; j<(*pp_numRegsInLevel)[sLevel]; j++) {

	  Region& sonRegion = sonLevel[j];
	  int sid = (*pp_numRegsBeforeLevel)[sLevel] + j;
	  
	  Nodes inter_nodes;
	  parentRegion.intersection(sonRegion,inter_nodes);

	  // if the son is a subset of the parent
	  if (inter_nodes.size() == sonRegion.size()) {
	    if (directSon) { // define adjacency
	      (*pp_adjMat)[pid].push_back(sid);
	      (*pp_adjMat)[sid].push_back(pid);
	      pp_parents[sid]++;
	    }
	    // and calculate double counting
	    pp_doubleCount[sid] += pp_doubleCount[pid];
	  }
	}
      }
    }

    // fix double counting
    for (int j=0; j<(*pp_numRegsInLevel)[sLevel]; j++) {
      int sid = (*pp_numRegsBeforeLevel)[sLevel] + j;
      pp_doubleCount[sid] = 1.0 - pp_doubleCount[sid];
    }
  }

  pp_regions_mrf = new MRF(*pp_adjMat);
}

void GBPPreProcessor::createConvexRegionsAdj() {

  // define neighbours

  pp_adjMat->resize(pp_numAllRegs);
  
  pp_parents = new int[pp_numAllRegs]; // will be cleared after calculating behte
  int* sons = new int[pp_numAllRegs]; // will be cleared after calculating the counting numbers
  for (int i=0; i<pp_numAllRegs; i++) {
    pp_parents[i] = 0;
    sons[i] = 0;
  }

  // and calculate double counts
  
  pp_doubleCount = new double[pp_numAllRegs]; // will be cleared after calculating bethe
  for (int i=0; i<pp_numBigRegs; i++) {
    pp_doubleCount[i] = 1.0;
  }
  for (int i=pp_numBigRegs; i<pp_numAllRegs; i++) {
    pp_doubleCount[i] = 0.0;
  }

  for (int sLevel=1; sLevel<pp_numRegsLevels; sLevel++) {

    RegionLevel& sonLevel = (*pp_allRegions)[sLevel];

    int pLevel = sLevel - 1;
    RegionLevel& parentLevel = (*pp_allRegions)[pLevel];
    for (int i=0; i<(*pp_numRegsInLevel)[pLevel]; i++) {

      Region& parentRegion = parentLevel[i];
      int pid = (*pp_numRegsBeforeLevel)[pLevel] + i;

      for (int j=0; j<(*pp_numRegsInLevel)[sLevel]; j++) {

	Region& sonRegion = sonLevel[j];
	int sid = (*pp_numRegsBeforeLevel)[sLevel] + j;
	  
	Nodes inter_nodes;
	parentRegion.intersection(sonRegion,inter_nodes);

	// if the son is a subset of the parent
	if (inter_nodes.size() == sonRegion.size()) {
	  (*pp_adjMat)[pid].push_back(sid);
	  (*pp_adjMat)[sid].push_back(pid);
	  pp_parents[sid]++;
	  sons[pid]++;
	}
      }
    }
  }
  // calculate double counting
  for (int i=pp_numBigRegs; i<pp_numAllRegs; i++) {
    for (int n=0; n<(*pp_adjMat)[i].size(); n++) {
      int j = (*pp_adjMat)[i][n];
      if (j<i) {
	double c_ij = pp_doubleCount[j] / sons[j];
	pp_doubleCount[i] += c_ij;
      }
    }
  }
  for (int i=pp_numBigRegs; i<(*pp_numRegsBeforeLevel)[pp_numRegsLevels-1]; i++) {
    pp_doubleCount[i] = 0.0;
  }

  for (int i=(*pp_numRegsBeforeLevel)[pp_numRegsLevels-1]; i<pp_numAllRegs; i++) {
    pp_doubleCount[i] = - pp_doubleCount[i];
  }

  delete[] sons;
  sons = 0;

  pp_regions_mrf = new MRF(*pp_adjMat);
}

void GBPPreProcessor::build2LayersRegionGraph() {

  // create all levels of regions
  
  RegionLevel* currLevel = &((*pp_allRegions)[0]);
  RegionLevel nextLevel;
  nextLevel.clear();
  currLevel->intersections(nextLevel);

  // add this level to the region graph
  pp_numRegsBeforeLevel->push_back(pp_numAllRegs);
  pp_numAllRegs += nextLevel.size();
  pp_numRegsLevels++;
  pp_allRegions->push_back(nextLevel);
  pp_numRegsInLevel->push_back(nextLevel.size());

  // initialize the pointers array to the regions
  pp_allRegionsPtr = new Region*[pp_numAllRegs];
  int i=0;
  for (int level=0; level<pp_numRegsLevels; level++) {
    RegionLevel& regLevel = (*pp_allRegions)[level];
    for (int r=0; r<(*pp_numRegsInLevel)[level]; r++) {
      pp_allRegionsPtr[i] = &(regLevel[r]);
      i++;
    }
  }

}

void GBPPreProcessor::create2LayersRegionsAdj() {

  // define neighbours

  pp_adjMat->resize(pp_numAllRegs);
  
  pp_parents = new int[pp_numAllRegs]; // will be cleared after calculating bethe
  for (int i=0; i<pp_numAllRegs; i++) {
    pp_parents[i] = 0;
  }

  // and calculate double counts
  
  pp_doubleCount = new double[pp_numAllRegs]; // will be cleared after calculating behte
  for (int i=0; i<pp_numBigRegs; i++) {
    pp_doubleCount[i] = 1.0;
  }
  if (pp_countingNode != 0) {
    for (int i=pp_numBigRegs; i<pp_numAllRegs; i++) {
      Region& reg_i = *(pp_allRegionsPtr[i]);
      if (reg_i.size() > 0) {
	pp_doubleCount[i] = pp_countingNode[reg_i[0]];
      }
      else {
	mexPrintf("error: reg no. %d is empty\n",i);
	pp_doubleCount[i] = 0.0;
      }
    }
  }
  else {
    for (int i=pp_numBigRegs; i<pp_numAllRegs; i++) {
      pp_doubleCount[i] = 0.0;
    }
  }

  for (int sLevel=1; sLevel<pp_numRegsLevels; sLevel++) {

    RegionLevel& sonLevel = (*pp_allRegions)[sLevel];

    for (int pLevel=0; pLevel<sLevel; pLevel++) {

      RegionLevel& parentLevel = (*pp_allRegions)[pLevel];
      bool directSon = (sLevel==(pLevel+1));

      for (int i=0; i<(*pp_numRegsInLevel)[pLevel]; i++) {
	Region& parentRegion = parentLevel[i];
	int pid = (*pp_numRegsBeforeLevel)[pLevel] + i;

	for (int j=0; j<(*pp_numRegsInLevel)[sLevel]; j++) {

	  Region& sonRegion = sonLevel[j];
	  int sid = (*pp_numRegsBeforeLevel)[sLevel] + j;
	  
	  Nodes inter_nodes;
	  parentRegion.intersection(sonRegion,inter_nodes);

	  // if the son is a subset of the parent
	  if (inter_nodes.size() == sonRegion.size()) {
	    if (directSon) { // define adjacency
	      (*pp_adjMat)[pid].push_back(sid);
	      (*pp_adjMat)[sid].push_back(pid);
	      pp_parents[sid]++;
	    }
	  }
	}
      }
    }

  }

  pp_regions_mrf = new MRF(*pp_adjMat);
}

void GBPPreProcessor::calcParentsAndDoubleCount() {

  // calculate parents  
  pp_parents = new int[pp_numAllRegs]; // will be cleared after calculating behte
  for (int i=0; i<pp_numAllRegs; i++) {
    pp_parents[i] = 0;
    for (int n=0; n<pp_regions_mrf->neighbNum(i); n++) {
      int j = (*pp_adjMat)[i][n];
      if (j<i) {
	pp_parents[i]++;
      }
      else {
	break;
      }
    }
  }

  // and calculate double counts
  
  pp_doubleCount = new double[pp_numAllRegs]; // will be cleared after calculating behte
  for (int i=0; i<pp_numBigRegs; i++) {
    pp_doubleCount[i] = 1.0;
  }
  for (int i=pp_numBigRegs; i<pp_numAllRegs; i++) {
    pp_doubleCount[i] = 0.0;
    Region& sonRegion = *(pp_allRegionsPtr[i]);
    
    // add double counts of grand-parents
    int firstDirectParent = (*pp_adjMat)[i][0];
    for (int j=0; j<firstDirectParent; j++) {
      Region& parentRegion = *(pp_allRegionsPtr[j]);
      Nodes inter_nodes;
      parentRegion.intersection(sonRegion,inter_nodes);
      // if the son is a subset of the parent
      if (inter_nodes.size() == sonRegion.size()) {
	pp_doubleCount[i] += pp_doubleCount[j];
      }
    }

    // add double counts of direct parents
    for (int n=0; n<pp_parents[i]; n++) {
      int j = (*pp_adjMat)[i][n];
      pp_doubleCount[i] += pp_doubleCount[j];
    }

    // fix double counts
    pp_doubleCount[i] = 1.0 - pp_doubleCount[i];
  }
  
}

void GBPPreProcessor::calcBethe() {
  
  pp_bethe = new double[pp_numAllRegs];
  for (int i=0; i<pp_numBigRegs; i++) {
    pp_bethe[i] = 1.0;
  }
  for (int i=pp_numBigRegs; i<pp_numAllRegs; i++) {
    double q = (1.0 - pp_doubleCount[i]) / max(1.0, double(pp_parents[i]));
    pp_bethe[i] = 1.0 / (2.0 - q);
  }
  
  delete[] pp_doubleCount;
  pp_doubleCount = 0;
  delete[] pp_parents;
  pp_parents = 0;
}

void GBPPreProcessor::calcNumOfStates() {
  int const* V = pp_mrf->V;

  for (int i=0; i<pp_numAllRegs; i++) {
    Region* region = pp_allRegionsPtr[i];

    // calc i's number of possible states
    int i_numStates = 1;
    
    nodes_iterator reg_niter = region->begin();
    while (reg_niter != region->end()) {
      int& node = (*reg_niter);
      i_numStates *= V[node];

      reg_niter++;
    }

    pp_regions_mrf->V[i] = i_numStates;

  }
}

void GBPPreProcessor::calcAssignTable() {

  int const& N = pp_mrf->N;
  int const* V = pp_mrf->V;
  
  pp_assignTable = new int**[pp_numAllRegs];
  for (int i=0; i<pp_numAllRegs; i++) {
    pp_assignTable[i] = new int*[pp_regions_mrf->neighbNum(i)];

    Region* i_region = pp_allRegionsPtr[i];
    int i_numStates = pp_regions_mrf->V[i];
    
    Assignment assignment;
    assignment.clear();
    assignment.resize(N);

    for (int n=0; n<pp_regions_mrf->neighbNum(i); n++) {
      
      int j = (*pp_adjMat)[i][n];
      
      pp_assignTable[i][n] = 0;
      
      if (i<j) {

	Region* j_region = pp_allRegionsPtr[j];
	
	pp_assignTable[i][n] = new int[i_numStates];
	for (int i_assignInd=0; i_assignInd<i_numStates; i_assignInd++) {
	  i_region->indToAssign(i_assignInd, assignment, V);
	  int j_assignInd = j_region->assignToInd(assignment, V);
	  pp_assignTable[i][n][i_assignInd] = j_assignInd;
	}
      }
    }
    
  }
}

void GBPPreProcessor::calcPotentials(Potentials* bigRegsPot) {

  int* regCard = pp_regions_mrf->V;

  pp_regions_mrf->initLocalPotentials();
  Potentials* regPot = pp_regions_mrf->localMat;
  if (pp_mrf->logspace) {
    for (int i=0; i<pp_numAllRegs; i++) {
      for (int xi=0; xi<regCard[i]; xi++) {
	regPot[i][xi] = 0.0;
      }
    }
  }
  else {
    for (int i=0; i<pp_numAllRegs; i++) {
      for (int xi=0; xi<regCard[i]; xi++) {
	regPot[i][xi] = 1.0;
      }
    }
  }
  int const& N = pp_mrf->N;
  int const* V = pp_mrf->V;
  
  pp_extractSingle = new int[N];
  pp_extractPairs = new int*[N];
  for (int i=0; i<N; i++) {
    pp_extractSingle[i] = -1;
    pp_extractPairs[i] = new int[pp_mrf->neighbNum(i)];
    for (int n=0; n<pp_mrf->neighbNum(i); n++) {
      pp_extractPairs[i][n] = -1;
    }
  }

  int numLocalUsed = 0;
  bool* recievedPot = new bool[pp_numBigRegs];
  for (int i=0; i<pp_numBigRegs; i++) {
    recievedPot[i] = false;
  }  

  Assignment assignment;
  assignment.clear();
  assignment.resize(N);

  // deal the local potentials of the nodes between the big regions  

  for (int i=0; i<pp_numBigRegs; i++) {
    
    Region* region = pp_allRegionsPtr[i];
    const_nodes_iterator niter = region->begin();
    while (niter != region->end()) {
      int n = *niter;
      if (pp_extractSingle[n]==-1) {
	if (pp_mrf->isLocal()) {
	  for (int i_assignInd=0; i_assignInd<regCard[i]; i_assignInd++) {
	    region->indToAssign(i_assignInd, assignment, V);
	    int xn = assignment[n];
	    if (pp_mrf->logspace) {
	      // 	    cout << "i = " << i << " n = " << n << " exp(-pp_mrf->localMat[n][xn]) = " << exp(-pp_mrf->localMat[n][xn]) << endl;
	      regPot[i][i_assignInd] += pp_mrf->localMat[n][xn];
	    }
	    else {
	      // 	    cout << "i = " << i << " n = " << n << " pp_mrf->localMat[n][xn] = " << pp_mrf->localMat[n][xn] << endl;
	      regPot[i][i_assignInd] *= pp_mrf->localMat[n][xn];
	    }
	  }
	}
	numLocalUsed++;

	recievedPot[i] = true;
	pp_extractSingle[n] = i;
      }
      niter++;
    }

    if (numLocalUsed == N) {
      break;
    }
  }
  // deal the pairwise potentials of the nodes between the big regions
  for (int i=0; i<N; i++) {
    for (int n=0; n<pp_mrf->neighbNum(i); n++) {
      int j = pp_mrf->adjMat[i][n];
      if (i<j) {
	for (int reg=0; reg<pp_numBigRegs; reg++) {
	  Region* region = pp_allRegionsPtr[reg];

	  if (region->contains(i) && region->contains(j)) {
	    if (pp_mrf->isPairwise()) {
	      for (int reg_assignInd=0; reg_assignInd<regCard[reg]; reg_assignInd++) {
		region->indToAssign(reg_assignInd, assignment, V);
		int xi = assignment[i];
		int xj = assignment[j];
		if (pp_mrf->logspace) {
		  regPot[reg][reg_assignInd] += pp_mrf->pairPotential(i,n,xi,xj);
		}
		else {
		  regPot[reg][reg_assignInd] *= pp_mrf->pairPotential(i,n,xi,xj);
		}
	      }
	    }
	    
	    recievedPot[reg] = true;
	    pp_extractPairs[i][n] = reg;
	    break;
	  }
	}
      }
    }
  }
  // deal the regions potentials, if given
  if (bigRegsPot != 0) {
    for (int reg=0; reg<pp_numBigRegs; reg++) {
      for (int reg_assignInd=0; reg_assignInd<regCard[reg]; reg_assignInd++) {
	if (pp_mrf->logspace) {
	  regPot[reg][reg_assignInd] += bigRegsPot[reg][reg_assignInd];
	}
	else {
	  regPot[reg][reg_assignInd] *= bigRegsPot[reg][reg_assignInd];
	}
	recievedPot[reg] = true;
      }
    }
  }
  
  // check wether exists a big region which is redundant, and issue an
  // error message

  for (int i=0; i<pp_numBigRegs; i++) {
    if (!recievedPot[i]) {
      //      cerr << "region number " << i << " is redundant: " << (*pp_allRegionsPtr[i]);
      mexPrintf("Warning: region number %d is redundant: ",i);
      Region* region = pp_allRegionsPtr[i];
      const_nodes_iterator niter = region->begin();
      while (niter != region->end()) {
	mexPrintf(" %d",*niter);
	niter++;
      }
      mexPrintf("\n");
    }
  }
  
  delete[] recievedPot;
  recievedPot = 0;

}

void GBPPreProcessor::extractSingle(double** regsBeliefs,
				    double** singleBeliefs,
				    SumOrMax smFlag) const {

  if (pp_mrf->logspace) {
    extractSingleLog(regsBeliefs,singleBeliefs,smFlag);
  }
  else {
    int const& N = pp_mrf->N;
    int const* V = pp_mrf->V;
    int* regCard = pp_regions_mrf->V;
  
    Assignment assignment;
    assignment.clear();
    assignment.resize(N);

    for (int n=0; n<N; n++) {
      int reg = pp_extractSingle[n];
      if (reg<0) {
	for (int xn=0; xn<V[n]; xn++) {
	  singleBeliefs[n][xn] = 1.0 / V[n];
	}
      }
      else {
	Region* region = pp_allRegionsPtr[reg];
	switch (smFlag) {
	  case SUM:
	    for (int xn=0; xn<V[n]; xn++) {
	      singleBeliefs[n][xn] = 0.0;
	    }
	    break;
	  case MAX:
	    for (int xn=0; xn<V[n]; xn++) {
	      singleBeliefs[n][xn] = -1.0;
	    }
	    break;
	  default:
	    break;
	}

	for (int reg_assignInd=0; reg_assignInd<regCard[reg]; reg_assignInd++) {
	  region->indToAssign(reg_assignInd, assignment, V);
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
	  for (int xn=0; xn<V[n]; xn++) {
	    sum_bel_n += singleBeliefs[n][xn];      
	  }
	  for (int xn=0; xn<V[n]; xn++) {
	    if (sum_bel_n > 0.0) {
	      singleBeliefs[n][xn] /= sum_bel_n;
	    }
	  }    
	}
      }
    }
  }
}

void GBPPreProcessor::extractSingleLog(double** regsBeliefs,
				       double** singleBeliefs,
				       SumOrMax smFlag) const {

  int const& N = pp_mrf->N;
  int const* V = pp_mrf->V;
  int* regCard = pp_regions_mrf->V;
  
  Assignment assignment;
  assignment.clear();
  assignment.resize(N);


  for (int n=0; n<N; n++) {
    int reg = pp_extractSingle[n];
    if (reg<0) {
      for (int xn=0; xn<V[n]; xn++) {
	singleBeliefs[n][xn] = 0.0;
      }
    }
    else{
      Region* region = pp_allRegionsPtr[reg];
      switch (smFlag) {
	case SUM:
	  for (int xn=0; xn<V[n]; xn++) {
	    singleBeliefs[n][xn] = -HUGE_VAL;
	  }
	  break;
	case MAX:
	  for (int xn=0; xn<V[n]; xn++) {
	    singleBeliefs[n][xn] = HUGE_VAL;
	  }
	  break;
	default:
	  break;
      }

      for (int reg_assignInd=0; reg_assignInd<regCard[reg]; reg_assignInd++) {
	region->indToAssign(reg_assignInd, assignment, V);
	int xn = assignment[n];
	switch (smFlag) {
	  case SUM:
	    singleBeliefs[n][xn] = AddLogFactor(singleBeliefs[n][xn],regsBeliefs[reg][reg_assignInd],-pp_mrf->getTemperature());
	    break;
	  case MAX:
	    if (regsBeliefs[reg][reg_assignInd] < singleBeliefs[n][xn]) {
	      singleBeliefs[n][xn] = regsBeliefs[reg][reg_assignInd];
	    }
	    break;
	  default:
	    break;
	}
      }

      // normalize
      double norm_bel_n = HUGE_VAL;
      double sum_bel_n = 0.0;
      for (int xn=0; xn<V[n]; xn++) {
	if (singleBeliefs[n][xn] < norm_bel_n) {
	  norm_bel_n = singleBeliefs[n][xn];
	}
      }
      for (int xn=0; xn<V[n]; xn++) {
	singleBeliefs[n][xn] -= norm_bel_n;
      }
    }
  }
}

void GBPPreProcessor::extractPairs(double** regsBeliefs,
				   double**** pairBeliefs,
				   SumOrMax smFlag) const {
  if (pp_mrf->logspace) {
    extractPairsLog(regsBeliefs,pairBeliefs,smFlag);
  }
  else {
    int const& N = pp_mrf->N;
    int const* V = pp_mrf->V;
    int* regCard = pp_regions_mrf->V;
  
    Assignment assignment;
    assignment.clear();
    assignment.resize(N);


    for (int i=0; i<N; i++) {
      for (int n=0; n<pp_mrf->neighbNum(i); n++) {
	int j = pp_mrf->adjMat[i][n];
	if (i<j) {
	  int reg = pp_extractPairs[i][n];
	  if (reg<0) {
	    for (int xi=0; xi<V[i]; xi++) {
	      for (int xj=0; xj<V[j]; xj++) {
		pairBeliefs[i][n][xi][xj] = 1.0 / (V[i] * V[j]);
	      }
	    }
	  }
	  else {
	    Region* region = pp_allRegionsPtr[reg];
	
	    for (int xi=0; xi<V[i]; xi++) {
	      for (int xj=0; xj<V[j]; xj++) {
	    
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
	      region->indToAssign(reg_assignInd, assignment, V);
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
	      for (int xi=0; xi<V[i]; xi++) {
		for (int xj=0; xj<V[j]; xj++) {
		  sum_bel_ij += pairBeliefs[i][n][xi][xj];
		}
	      }
	      for (int xi=0; xi<V[i]; xi++) {
		for (int xj=0; xj<V[j]; xj++) {
		  pairBeliefs[i][n][xi][xj] /= sum_bel_ij;
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

void GBPPreProcessor::extractPairsLog(double** regsBeliefs,
				      double**** pairBeliefs,
				      SumOrMax smFlag) const {
  int const& N = pp_mrf->N;
  int const* V = pp_mrf->V;
  int* regCard = pp_regions_mrf->V;
  
  Assignment assignment;
  assignment.clear();
  assignment.resize(N);


  for (int i=0; i<N; i++) {
    for (int n=0; n<pp_mrf->neighbNum(i); n++) {
      int j = pp_mrf->adjMat[i][n];
      if (i<j) {
	int reg = pp_extractPairs[i][n];
	if (reg<0) {
	  for (int xi=0; xi<V[i]; xi++) {
	    for (int xj=0; xj<V[j]; xj++) {
	      pairBeliefs[i][n][xi][xj] = 0.0;
	    }
	  }
	}
	else {
	  Region* region = pp_allRegionsPtr[reg];
	
	  for (int xi=0; xi<V[i]; xi++) {
	    for (int xj=0; xj<V[j]; xj++) {
	    
	      switch (smFlag) {
		case SUM:
		  pairBeliefs[i][n][xi][xj] = -HUGE_VAL;
		  break;
		case MAX:
		  pairBeliefs[i][n][xi][xj] = HUGE_VAL;
		  break;
		default:
		  break;
	      }
	    }
	  }
	  for (int reg_assignInd=0; reg_assignInd<regCard[reg]; reg_assignInd++) {
	    region->indToAssign(reg_assignInd, assignment, V);
	    int xi = assignment[i];
	    int xj = assignment[j];
	    switch (smFlag) {
	      case SUM:
		pairBeliefs[i][n][xi][xj] = AddLogFactor(pairBeliefs[i][n][xi][xj],regsBeliefs[reg][reg_assignInd],
							 -pp_mrf->getTemperature());
		break;
	      case MAX:
		if (regsBeliefs[reg][reg_assignInd] < pairBeliefs[i][n][xi][xj]) {
		  pairBeliefs[i][n][xi][xj] = regsBeliefs[reg][reg_assignInd];
		}
		break;
	      default:
		break;
	    }
	  }
      
	  // normalize for max-gbp
	  double norm_bel_ij = HUGE_VAL;
	  for (int xi=0; xi<V[i]; xi++) {
	    for (int xj=0; xj<V[j]; xj++) {
	      if (pairBeliefs[i][n][xi][xj] < norm_bel_ij) {
		norm_bel_ij = pairBeliefs[i][n][xi][xj];
	      }
	    }
	  }
	  for (int xi=0; xi<V[i]; xi++) {
	    for (int xj=0; xj<V[j]; xj++) {
	      pairBeliefs[i][n][xi][xj] -= norm_bel_ij;
	    }
	  }
	}
      }
    }
  }
}

void GBPPreProcessor::printAllRegions() const {

  for (int i=0; i<pp_numAllRegs; i++) {
    cout << i << '\t' << (*pp_allRegionsPtr[i]);
  }
  cout << endl;
}

void GBPPreProcessor::printAdj() const {

  for (int i=0; i<pp_numAllRegs; i++) {
    Nodes& reg = (*pp_adjMat)[i];
    cout << i << '\t';
    for (unsigned int j=0; j<reg.size(); j++) {
      cout << reg[j] << ' ';
    }
    cout << endl;
  }
  cout << endl;
}

void GBPPreProcessor::printBethe() const {
  if (pp_bethe != 0) {
    for (int i=0; i<pp_numAllRegs; i++) {
      cout << i << " bethe: " << pp_bethe[i] << endl;
    }
    cout << endl;
  }
}
