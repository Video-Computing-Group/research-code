#include "Region.h"

Region::Region() {
  r_nodes.clear();
  r_size = 0;
}

void Region::assignNodes(Nodes const& nodes) {
  r_nodes.clear();
  const_nodes_iterator niter = nodes.begin();
  while (niter != nodes.end()) {
    r_nodes.push_back(*niter);
    niter++;
  }
  r_size = r_nodes.size();
}

void Region::intersection(Region const& otherReg, Nodes& inter_nodes) const {
  // find the intersection - assuming the nodes in each region are sorted
  inter_nodes.clear();
  
  const_nodes_iterator this_node = begin();
  const_nodes_iterator other_node = otherReg.begin();
  
  while (1) {
    while ((this_node != end()) &&
	   ((*this_node) < (*other_node))) {
      this_node++;
    }
    if (this_node == end()) {
      break;
    }
    while ((other_node != otherReg.end()) &&
	   ((*other_node) < (*this_node))) {
      other_node++;
    }
    if (other_node == otherReg.end()) {
      break;
    }
    if ((*this_node) == (*other_node)) {
      inter_nodes.push_back(*this_node);
      this_node++;
      other_node++;
    }	
  }
      
}

bool Region::contains(int node) const {

  const_nodes_iterator iter = begin();
  while (iter != end()) {
    if ((*iter)==node) {
      return true;
    }
    iter++;
  }
  return false;
}

void Region::indToAssign(int const& index, Assignment& assign, int const* card) const {
  // assign should be initialized outside to the size of number of nodes
  int curr_ind = index;
  for (int i=r_size-1; i>=0; i--) {
    int factor = card[r_nodes[i]];
    int i_assign = curr_ind % factor;
    curr_ind = curr_ind / factor;
    assign[r_nodes[i]] = i_assign;
  }
}

int Region::assignToInd(Assignment const& assign, int const* card) const {
  int index = 0 ;
  int factor = 1;
  for (int i=r_size-1 ; i>=0 ; i--) {
    index += (factor*assign[r_nodes[i]]);
    factor*= card[r_nodes[i]];
  }
  return index;  
}

bool Region::operator==(Region const& otherReg) const {

  // assuming both regions are sorted
  
  if (r_size != otherReg.size()) {
    return false;
  }
  
  for (unsigned int i=0; i<r_size; i++) {
    if (r_nodes[i] != otherReg[i]) {
      return false;
    }
  }

  return true;
}

bool Region::operator==(Nodes const& nodes) const {

  // assuming both regions are sorted
  
  if (r_size != nodes.size()) {
    return false;
  }
  
  for (unsigned int i=0; i<r_size; i++) {
    if (r_nodes[i] != nodes[i]) {
      return false;
    }
  }

  return true;
}

ostream& operator<<(ostream& os, Region const& reg) {
  for (unsigned int i=0; i<reg.size(); i++) {
    os << reg.r_nodes[i] << ' ';
  }
  os << endl;
  return os;
}
