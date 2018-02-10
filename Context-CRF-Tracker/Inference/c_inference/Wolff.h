#include "MonteCarlo.h"
#include "PottsMRF.h"

#ifndef __WOLFF__
#define __WOLFF__

class Wolff : public MonteCarlo {

  /**
     This class makes inference using Wolff sampling method,
     where in each steps one cluster can change its state
     This method is specific for q-state Potts model
   
     Part of the c_inference package
     @version November 2004
     @author Talya Meltzer
  */
  
 public:

  // ctor
  Wolff(PottsMRF const* mrf, int* startX, int burningTime, int samplingInterval, int S) :
    MonteCarlo(mrf,startX,burningTime,samplingInterval,S) {}
  
  virtual ~Wolff() {} // dtor

 protected:

  virtual void transition();
  int cluster(int i, int i_nei, int new_xi);
  bool accept(int i, int i_nei);
  
};

#endif
