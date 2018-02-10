#include "MonteCarlo.h"
#include "PottsMRF.h"

#ifndef __SWENDSEN_WANG__
#define __SWENDSEN_WANG__

class SwendsenWang : public MonteCarlo {

  /**
     This class makes inference using Swendsen-Wang sampling method,
     where in each steps number of clusters can change their states
     This method is specific for q-state Potts model
   
     Part of the c_inference package
     @version November 2004
     @author Talya Meltzer
  */
  
 public:

  // ctor
  SwendsenWang(PottsMRF const* mrf, int* startX, int burningTime, int samplingInterval, int S) :
    MonteCarlo(mrf,startX,burningTime,samplingInterval,S)
    { init(); }
  
  virtual ~SwendsenWang(); // dtor

 protected:

  bool** sw_bondFrozen;
  int* sw_cluster;
  int* sw_labelLabel;
  int* sw_new;
  bool* sw_newChosen;

  virtual void transition();
  void init();
  void initializeClusterVariables();
  void freezeOrMeltBonds();
  int properLabel(int label);
  void labelClusters();
  void flipClusterNodes();

};

#endif
