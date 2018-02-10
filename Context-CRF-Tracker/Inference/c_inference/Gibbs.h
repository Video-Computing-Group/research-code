#include "MonteCarlo.h"

#ifndef __GIBBS__
#define __GIBBS__

class Gibbs : public MonteCarlo {

  /**
     This class makes inference using Gibbs sampling method,
     where in each step, only one node can change its state
   
     Part of the c_inference package
     @version November 2004
     @author Talya Meltzer
  */
  
 public:

  // ctor
  Gibbs(MRF const* mrf, int* startX, int burningTime, int samplingInterval, int S) :
    MonteCarlo(mrf,startX,burningTime,samplingInterval,S) {};
  
  virtual ~Gibbs() {}; // dtor

 protected:
  
  virtual void transition();
  
};

#endif
