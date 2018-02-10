#include "InferenceAlgorithm.h"
#include <math.h>

#ifndef __MEAN_FIELD__
#define __MEAN_FIELD__

class MeanField : public InferenceAlgorithm {

  /**
     This class makes inference using mean-field approximation
   
     Part of the c_inference package
     @version November 2004
     @author Talya Meltzer
  */
  
 public:

  // ctor
  MeanField(MRF const* mrf, int maxIter = 2000, double th = pow(10.,-8)) :
    InferenceAlgorithm(mrf), mf_maxIter(maxIter), mf_th(th) {}
  
  virtual ~MeanField() {} // dtor

  virtual double** inference(int* converged);
  
 protected:
  
  int mf_maxIter; // maximum number of iterations in inference
  double mf_th; // threshold for convergence
  
};

#endif
