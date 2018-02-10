#include "MeanField.h"

#ifndef __LOG_MEAN_FIELD__
#define __LOG_MEAN_FIELD__

class LogMeanField : public MeanField {

  /**
     This class makes inference using mean-field approximation in the log-space
   
     Part of the c_inference package
     @version March 2006
     @author Talya Meltzer
  */
  
 public:

  // ctor
  LogMeanField(MRF const* mrf, int maxIter = 2000, bool logBels = false, double th = log(pow(10.,-8))) :
    MeanField(mrf,maxIter,th) {lmf_logBels = logBels;}
  
  virtual ~LogMeanField() {} // dtor

  virtual double** inference(int* converged);

 protected:
    bool lmf_logBels; // return the beliefs in the -log-space
  
};

#endif
