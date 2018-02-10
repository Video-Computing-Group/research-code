#include "InferenceAlgorithm.h"
#include <vector>
#include <math.h>

#ifndef __GENERALIZED_BP__
#define __GENERALIZED_BP__

class GBP : public InferenceAlgorithm {

  /**
     This class makes inference using Jonathan Yedidia's GBP algorithm
   
     Part of the c_inference package
     @version November 2004
     @author Talya Meltzer
  */

 public:

  // ctor
  GBP(MRF const* reg_mrf, int*** assignInd, double* bethe,
      SumOrMax m = MAX, double alpha = 0.5, int maxIter = 2000,
      double*** initMsg = 0, double th = pow(10.,-8)) :

    InferenceAlgorithm(reg_mrf), gbp_sumOrMax(m), gbp_alpha(alpha), gbp_maxIter(maxIter), gbp_th(th)

    {
      gbp_bethe = bethe;
      gbp_assignInd = assignInd;
      gbp_messages = 0;
      initMessages(initMsg);
      initBeliefs();
    }

  virtual ~GBP(); // dtor

  virtual double** inference(int* converged);

  double*** getMessages() { return gbp_messages; }
  bool isMaxMarg(double epsilon) const;
  bool isSumMarg(double epsilon) const;
  
 protected:

  // data members
  
  SumOrMax gbp_sumOrMax; // use sum or max 
  double gbp_alpha; // averging parameter: alpha*newMessage + (1-alpha)*lastMessage
  int gbp_maxIter; // maximum number of iterations in inference
  double gbp_th; // threshold for convergence

  double* gbp_bethe; // relates to the double-counting of each region
  int*** gbp_assignInd; // convertion table between fathers' and sons' assignment indices
  double*** gbp_messages; // the messages from node to node
  vector<int>* gbp_arcsOrder; // the updating order while passing messages
  int gbp_numArcs;
  
  // protected methods
  
  virtual void initMessages(double*** initMsg);
  virtual void initBeliefs();
  void defineUpdatingOrder();
  virtual void calcIncomingMessages(double* incoming_i, int i, int j);  
  void normalize(double* dataVec, int Vj, double epsilon);
  void freeMessages();
};

#endif
