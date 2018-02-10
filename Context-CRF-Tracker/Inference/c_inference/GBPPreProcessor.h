#include "RegionLevel.h"
#include "MRF.h"
#include <vector>

#ifndef __GBP_PRE_PROCESSOR__
#define __GBP_PRE_PROCESSOR__

class GBPPreProcessor {

  /**
     This class creates the full region-graph needed for makeing inference
     using Jonathan Yedidia's GBP algorithm
   
     Part of the c_inference package
     @version November 2004
     @author Talya Meltzer
  */

  typedef vector<RegionLevel>::iterator level_iterator;
  typedef vector<RegionLevel>::const_iterator const_level_iterator;
  
 public:

  // ctors:

  // given the nodes' MRF and the big regions,
  // builds the region-graph and calculates the MRF of the regions,
  // the assign-table, and bethe
  GBPPreProcessor(RegionLevel const& bigRegions, MRF const* mrf, bool trw = false, 
		  bool full = true, double* countingNode = 0, Potentials* bigRegsPot = 0);
  
  // given the region-graph (defined by the regions and their adjacencies)
  // and the nodes' MRF, calculates the MRF of the regions (local potentials),
  // the assign-table, and bethe
  GBPPreProcessor(vector<RegionLevel>* allRegions, MRF const* mrf, bool trw = false,
		  bool full = true, double* countingNode = 0, Potentials* bigRegsPot = 0);

  ~GBPPreProcessor(); // dtor

  // get the fields needed for GBP inference algorithm
  // (use only after "convertToRegionMRF")
  MRF* getRegionMRF() { return pp_regions_mrf; }
  double* getBethe() { return pp_bethe; }
  int*** getAssignTable() { return pp_assignTable; }
  int* getExtractSingle() { return pp_extractSingle; }
  int** getExtractPairs() { return pp_extractPairs; }
  Region** getAllRegions() { return pp_allRegionsPtr; }

  // extract the given regions' beliefs into single-node beliefs
  void extractSingle(double** regsBeliefs,
		     double** singleBeliefs,
		     SumOrMax smFlag = MAX) const;
  void extractSingleLog(double** regsBeliefs,
			double** singleBeliefs,
			SumOrMax smFlag = MAX) const;
  // extract the given regions' beliefs into pairwise beliefs
  void extractPairs(double** regsBeliefs,
		    double**** pairBeliefs,
		    SumOrMax smFlag = MAX) const;
  void extractPairsLog(double** regsBeliefs,
		       double**** pairBeliefs,
		       SumOrMax smFlag = MAX) const;
  
  // debug methods
  void printAllRegions() const;
  void printAdj() const;
  void printBethe() const;

 private:

  // data members
  
  int* pp_parents;
  double* pp_doubleCount;
  double* pp_bethe;
  int*** pp_assignTable;

  MRF const* pp_mrf;
  MRF* pp_regions_mrf;
  vector<Nodes>* pp_adjMat;

  int pp_numBigRegs;
  int pp_numAllRegs;
  int pp_numRegsLevels;
  vector<RegionLevel>* pp_allRegions;
  Region** pp_allRegionsPtr;
  vector<int>* pp_numRegsInLevel;
  vector<int>* pp_numRegsBeforeLevel;
  double* pp_countingNode;
  
  int* pp_extractSingle;
  int** pp_extractPairs;

  // private methods
  void init();
  void buildFullRegionGraph();
  void createRegionsAdj();
  void createConvexRegionsAdj();
  void build2LayersRegionGraph();
  void create2LayersRegionsAdj();
  void calcParentsAndDoubleCount();
  void calcBethe();
  void calcNumOfStates();
  void calcAssignTable();
  void calcPotentials(Potentials* bigRegsPot);
};

#endif
