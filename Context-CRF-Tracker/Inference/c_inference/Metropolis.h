#include "MonteCarlo.h"
#include "PottsMRF.h"

#ifndef __METROPOLIS__
#define __METROPOLIS__

class Metropolis : public MonteCarlo {

  /**
     This class makes inference using Metropolis sampling method,
     where in each step, only one node can change its state.
     A flip is always done if the potential of the new state is
     bigger than the potential of the current one,
     else it's done in probability new_potential / curr_potential
   
     Part of the c_inference package
     @version November 2004
     @author Talya Meltzer
  */
  
 public:

  // ctor
  Metropolis(MRF const* mrf, int* startX, int burningTime, int samplingInterval, int S) :
    MonteCarlo(mrf,startX,burningTime,samplingInterval,S) {}
  
    virtual ~Metropolis() {} // dtor

 protected:
  
    virtual void transition();

};

#endif
