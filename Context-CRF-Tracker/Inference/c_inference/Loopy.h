#include "InferenceAlgorithm.h"
#include <math.h>

#ifndef __LOOPY__
#define __LOOPY__

class Loopy : public InferenceAlgorithm {

  /**
     This class makes inference using Loopy Belief Propagation algorithm
   
     Part of the c_inference package
     @version November 2004
     @author Talya Meltzer
  */
  
 public:

  // ctor
  Loopy(MRF const* mrf, SumOrMax m = MAX, Strategy s = SEQUENTIAL,
	int maxIter = 2000, double** trwRho = 0, double*** initMsg = 0,
	double th = pow(10.,-8)) :
    InferenceAlgorithm(mrf),
    l_strategy(s), l_sumOrMax(m), l_maxIter(maxIter), l_th(th)
    { l_messages = 0; l_pairBeliefs = 0; l_trwRho = 0;
      initMessages(initMsg); initPairBeliefs(); initTRWRho(trwRho); }

  virtual ~Loopy(); // dtor

  virtual double** inference(int* converged);
  virtual double** inferenceTRBP(int* converged);
  void initMessages(double*** initMsg);
  void initPairBeliefs();
  void initTRWRho(double** trwRho);
  virtual double**** calcPairBeliefs();
  virtual double**** calcPairBeliefsTRBP();
  
  double*** getMessages() { return l_messages; }
  
 protected:
  
  Strategy l_strategy; // strategy of updating
  double*** l_messages; // the messages from node to node
  double**** l_pairBeliefs; // the pairwise beliefs
  SumOrMax l_sumOrMax; // use sum or max 
  int l_maxIter; // maximum number of iterations in inference
  double l_th; // threshold for convergence
  double** l_trwRho; // rho values for each edge in trw mode

  void freeMessages();
  void freePairBeliefs();
  void freeTRWRho();
};

#endif
