#include "LoopySTime.h"

#ifndef __LOG_LOOPY_SAVE_TIME__
#define __LOG_LOOPY_SAVE_TIME__

class LogLoopySTime : public LoopySTime {

  /**
     This class makes inference using Loopy Belief Propagation algorithm
     in log-space and in mode of saving time but wasting space for
     the messages
   
     Part of the c_inference package
     @version August 2005
     @author Talya Meltzer
  */
  
 public:

  // ctor
  LogLoopySTime(MRF const* mrf, SumOrMax m = MAX, Strategy s = SEQUENTIAL,
		int maxIter = 2000, double** trwRho = 0, double*** initMsg = 0,
		bool logBels = false, double th = log(pow(10.,-8))) :
    LoopySTime(mrf,m,s,maxIter,trwRho,initMsg,th) {ll_logBels = logBels;}

    virtual ~LogLoopySTime() {} // dtor

    virtual double** inference(int* converged);
    virtual double** inferenceTRBP(int* converged);
    virtual double**** calcPairBeliefs();
    virtual double**** calcPairBeliefsTRBP();

 protected:
    bool ll_logBels; // return the beliefs in the -log-space
  
};

#endif
