#include "definitions.h"
#include "MRF.h"

#ifndef __INFERENCE_ALGORITHM__
#define __INFERENCE_ALGORITHM__

class InferenceAlgorithm {
  
  /**
     This class defines the interface for making inference
   
     Part of the c_inference package
     @version November 2004
     @author Talya Meltzer
  */
  
 public:
  // ctor
  InferenceAlgorithm(MRF const* mrf);
  virtual ~InferenceAlgorithm(); // dtor

  virtual double** inference(int* converged) = 0;
  void initBeliefs();
  
  void getMaximalStates(vector<int>& max_states, int i, double epsilon) const;
  
 protected:
  
  MRF const* ia_mrf; // holds all data about the variables: cardinality, neighbours, potentials
  
  double** ia_beliefs; // will contain the probabilities calculated for
  //                       xi=0 and xi=1 (for i=0,..,N-1)

  void freeBeliefs();
};

#endif
