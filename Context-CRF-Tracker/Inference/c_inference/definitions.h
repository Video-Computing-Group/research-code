#include <vector>

#ifndef __REGIONS_DEFINITIONS__
#define __REGIONS_DEFINITIONS__

/**
   enumerators and type-definitions used for BP algorithm (Loopy and GBP)
   
   Part of the c_inference package
   @version November 2004
   @author Talya Meltzer
*/

enum Model {GENERAL, POTTS};
enum Strategy {SEQUENTIAL, PARALLEL};
enum SumOrMax {SUM, MAX};
enum ArcNodes {FATHER, SON};

using namespace std;
typedef std::vector<int> Nodes;
typedef Nodes::iterator nodes_iterator;
typedef Nodes::const_iterator const_nodes_iterator;
typedef vector<int> Assignment;
typedef vector<int> CardVec;
typedef double Potential;
typedef Potential* Potentials;
typedef Potentials* PairPotentials;
typedef PairPotentials* PairPotentialsPtr;

class Region;
typedef vector<Region>::iterator region_iterator;
typedef vector<Region>::const_iterator const_region_iterator;

#endif
