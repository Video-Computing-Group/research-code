#include "PairsGBP.h"
#include <math.h>

#ifndef __LOG_PAIRS_GBP__
#define __LOG_PAIRS_GBP__

class LogPairsGBP : public PairsGBP {

  /**
     This class makes inference using the simple form of GBP for pairs-regions
     in the log-space
   
     Part of the c_inference package
     @version August 2006
     @author Talya Meltzer
  */
  
 public:

 public:

  // ctor
  LogPairsGBP(MRF const* mrf, SumOrMax m = MAX, Strategy s = SEQUENTIAL,
	      int maxIter = 2000, double* doubleCount = 0, double*** initMsg = 0,
	      bool logBels = false, double th = log(pow(10.,-8))) :
    PairsGBP(mrf,m,s,maxIter,doubleCount,initMsg,th) {lpgbp_logBels = logBels;}

    virtual ~LogPairsGBP() {} // dtor

    virtual double** inference(int* converged);
    virtual double**** calcPairBeliefs();

 protected:
    bool lpgbp_logBels; // return the beliefs in the -log-space
};

#endif
