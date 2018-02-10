#include "MRF.h"

#ifndef __POTTS_MODEL_MARKOV_RANDOM_FIELD__
#define __POTTS_MODEL_MARKOV_RANDOM_FIELD__

class PottsMRF : public MRF {

  /**
     This class holds all the data defining the Markov Random Field
     for Potts Model (holding only lambda_ij, where
     pairPotential(i,j) = exp(-lambda_ij) for xi=xj,
     and exp(+lambda_ij) for xi!=xj
   
     Part of the c_inference package
     @version November 2004
     @author Talya Meltzer
  */
  
 public:
  
  PottsMRF(vector<Nodes>& adj) : MRF(adj) { lambdaValues = 0; }
  
  virtual ~PottsMRF();

  virtual double pairPotential (int i, int n, int xi, int xj) const;
  virtual double pairEnergy (int i, int n, int xi, int xj) const;
  void initLambdaValues();
  virtual void setTemperature(double sT);
  virtual double getEnergy(int const* assignment) const;
  virtual bool isPairwise() const { return lambdaValues != 0; }
  
  double** lambdaValues; // defines the potentials between
  //                        2 neighbouring variables
  
}; //PottsMRF;

#endif
