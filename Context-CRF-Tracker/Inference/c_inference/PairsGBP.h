#include "InferenceAlgorithm.h"
#include <math.h>

#ifndef __PAIRS_GBP__
#define __PAIRS_GBP__

class PairsGBP : public InferenceAlgorithm {

  /**
     This class makes inference using the simple form of GBP for pairs-regions
   
     Part of the c_inference package
     @version August 2006
     @author Talya Meltzer
  */
  
 public:

  // ctor
  PairsGBP(MRF const* mrf, SumOrMax m = MAX, Strategy s = SEQUENTIAL,
	   int maxIter = 2000, double* doubleCount = 0,
	   double*** initMsg = 0, double th = pow(10.,-8)) :
    InferenceAlgorithm(mrf),
    pgbp_strategy(s), pgbp_sumOrMax(m), pgbp_maxIter(maxIter), pgbp_th(th)
    { pgbp_messages = 0; pgbp_pairBeliefs = 0; pgbp_beta = 0;
      initMessages(initMsg); initPairBeliefs(); initBeta(doubleCount); }

  virtual ~PairsGBP(); // dtor

  virtual double** inference(int* converged);
  void initMessages(double*** initMsg);
  void initPairBeliefs();
  void initBeta(double* doubleCount);
  virtual double**** calcPairBeliefs();
  
  double*** getMessages() { return pgbp_messages; }
  
 protected:
  
  Strategy pgbp_strategy; // strategy of updating
  double*** pgbp_messages; // the messages from node to node
  double**** pgbp_pairBeliefs; // the pairwise beliefs
  SumOrMax pgbp_sumOrMax; // use sum or max 
  int pgbp_maxIter; // maximum number of iterations in inference
  double pgbp_th; // threshold for convergence
  double* pgbp_beta; // beta_i = 1 - d_i - c_i

  void freeMessages();
  void freePairBeliefs();
  void freeBeta();
};

#endif
