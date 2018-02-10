#include "GBP.h"
#include <vector>
#include <math.h>

#ifndef __LOG_GENERALIZED_BP__
#define __LOG_GENERALIZED_BP__

class LogGBP : public GBP {

  /**
     This class makes inference using Jonathan Yedidia's GBP algorithm in log-space
   
     Part of the c_inference package
     @version November 2004
     @author Talya Meltzer
  */

 public:

  // ctor
  LogGBP(MRF const* reg_mrf, int*** assignInd, double* bethe,
      SumOrMax m = MAX, double alpha = 0.5, int maxIter = 2000,
      double*** initMsg = 0, bool logBels = false, double th = pow(10.,-8)) :
    GBP(reg_mrf, assignInd, bethe, m, alpha, maxIter, initMsg, th)

    {
      lgbp_logBels = logBels;
      initMessages(initMsg);
      initBeliefs();
    }

    virtual ~LogGBP() {} // dtor

    virtual double** inference(int* converged);

 protected:

    // data members

    bool lgbp_logBels;  // return the beliefs in the -log-space
  
    // protected methods
  
    virtual void initMessages(double*** initMsg);
    virtual void initBeliefs();
    virtual void calcIncomingMessages(double* incoming_i, int i, int j);  
    void minNormalize(double* dataVec, int Vj);
    void normalizeBeliefs();
};

#endif
