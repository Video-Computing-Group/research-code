#include "definitions.h"

#ifndef __MARKOV_RANDOM_FIELD__
#define __MARKOV_RANDOM_FIELD__

class MRF {

  /**
     This class holds all the data defining the Markov Random Field
   
     Part of the c_inference package
     @version November 2004
     @author Talya Meltzer
  */
  
 public:
  
  MRF(vector<Nodes>& adj) :
    adjMat(adj) {
    
    N = adjMat.size();
    V = new int[N];
    for (int i=0; i<N; i++) {
      V[i] = 0;
    }
    lambdaMat = 0;
    localMat = 0;

    logspace = false;
    T = 1.0;
  }
  
  virtual ~MRF();

  inline int neighbNum(int i) const { return adjMat[i].size(); }
  virtual double pairPotential (int i, int n, int xi, int xj) const;
  virtual double pairEnergy (int i, int n, int xi, int xj) const;
  void assignPairPotential(int i, int nj, PairPotentials& pairPot);
  void initPairPotentials();
  void initLocalPotentials();
  virtual void setTemperature(double sT);
  double getTemperature() const { return T;}
  virtual double getEnergy(int const* assignment) const;
  bool isLocal() const { return localMat != 0; }
  virtual bool isPairwise() const { return lambdaMat != 0; }
  
  int N; // no. of variables
  int* V; // cardinality - no. of possible values for each variable

  vector<Nodes>& adjMat; // neighbouring matrix
  PairPotentials** lambdaMat; // defines the potentials between
  //                             2 neighbouring variables
  Potentials* localMat; // defines the local potentials of the variables
  bool logspace;

 protected:
  double T; // temperature
  
}; //MRF;

#endif
