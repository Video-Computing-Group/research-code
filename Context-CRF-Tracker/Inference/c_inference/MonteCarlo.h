#include "InferenceAlgorithm.h"

#ifndef __MONTE_CARLO__
#define __MONTE_CARLO__

class MonteCarlo : public InferenceAlgorithm {

  /**
     This class defines the interface for making inference using
     Monte-Carlo sampling method
   
     Part of the c_inference package
     @version November 2004
     @author Talya Meltzer
  */
  
 public:

  // ctor
  MonteCarlo(MRF const* mrf, int* startX, int burningTime, int samplingInterval, int S) :
    InferenceAlgorithm(mrf),
    mc_burning(burningTime), mc_interval(samplingInterval),
    mc_S(S)
    { mc_currX = 0; init(startX); }
  
  virtual ~MonteCarlo(); // dtor
  
  virtual void reachEquilibrium();
  virtual double** inference(int* converged);

  void getCurrentState(int* currX) const;
  
  void init(int* startX);
  void initState(int* startX);
  void setBurningTime(int burningTime) { mc_burning = burningTime; }
  void setSamplingInterval(int samplingInterval) { mc_interval = samplingInterval; }
  void setSamplesNum(int S) { mc_S = S; }

 protected:
  
  int mc_burning; // burning time to wait before sampling
  int mc_interval; // time to wait between samples
  int mc_S; // number of samples to take

  int* mc_currX; // the current state

  double* mc_Pi;
  double** mc_Pxi;

  virtual void transition() = 0;
  
  void freeState();
  int chooseInteger(double* q, int n);
};

#endif
