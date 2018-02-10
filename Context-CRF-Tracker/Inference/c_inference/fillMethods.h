#include "mex.h"
#include "MRF.h"
#include "PottsMRF.h"
#include "RegionLevel.h"

void fillAdjMat(const mxArray* mxAdj, vector<Nodes>& adjMat);
void fillLocalMat(const mxArray* mxLocal, MRF* mrf);
void fillLocalPot(const mxArray* mxLocal, Potentials* local, int* V, int& N);
void fillPsiMat(const mxArray* mxLambda, MRF* mrf);
void fillLambdaMat(const mxArray* mxLambda, PottsMRF* mrf);
void fillRegions(const mxArray* mxRegions, RegionLevel& regions);
void fillRegionLevels(const mxArray* mxAllRegs, vector<RegionLevel>& allRegions);
void fillInitialAssignment(const mxArray* mxAssign, int* startX, int num_nodes);
void fillRhoMat(const mxArray* mxRho, const MRF* mrf, double** rho);
void fillInt(const mxArray* mxIntArray, int* intArray, int arrayLen);
void fillDouble(const mxArray* mxDoubleArray, double* doubleArray, int arrayLen);
