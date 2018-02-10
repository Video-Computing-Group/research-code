#include "definitions.h"
#include <vector>
#include <ostream>

#ifndef __REGION__
#define __REGION__

class Region {

  /**
     This class defines a region in the graph as a collection of nodes
   
     Part of the c_inference package
     @version November 2004
     @author Talya Meltzer
  */
  
 public:

  // ctor
  Region();

  void assignNodes(Nodes const& nodes);
  
  unsigned int size() const { return r_size; }
  void intersection(Region const& otherReg, Nodes& inter_nodes) const;
  bool contains(int node) const;

  Nodes& nodes() { return r_nodes; }

  nodes_iterator begin() { return r_nodes.begin(); }
  nodes_iterator end() { return r_nodes.end(); }
  const_nodes_iterator begin() const { return r_nodes.begin(); }
  const_nodes_iterator end() const { return r_nodes.end(); }

  // indexing the possible assignments for the region according their cardinality:
  
  // convert index into assignment
  void indToAssign(int const& index, Assignment& assign, int const* card) const;
  // convert assignment into index
  int assignToInd(Assignment const& assign, int const* card) const;

  // operators
  int operator[](int i) const { return r_nodes[i]; }
  int& operator[](int i) { return r_nodes[i]; }
  bool operator==(Region const& otherReg) const;
  bool operator==(Nodes const& nodes) const;
  friend ostream& operator<<(ostream& os, Region const& reg);

 private:
  
  Nodes r_nodes;
  unsigned int r_size;
};

#endif
