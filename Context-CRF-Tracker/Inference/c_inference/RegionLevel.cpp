#include "RegionLevel.h"

bool RegionLevel::addRegion(Region& region) {

  region_iterator riter = begin();
  while (riter != end()) {
    Nodes inter_nodes;
    riter->intersection(region,inter_nodes);
    // check if the given region is a subset of an existing one
    if (inter_nodes.size() == region.size()) {
      return false;
    }
    // check if the existing region is a subset of the given one
    if (inter_nodes.size() == riter->size()) {
      // if so - remove this region
      riter = removeRegion(riter);
    }
    else {
      riter++;
    }
  }

  rl_regions.push_back(region);

  return true;
}

region_iterator RegionLevel::removeRegion(region_iterator pos) {
  return rl_regions.erase(pos);
}

void RegionLevel::intersections(RegionLevel& nextLevel) const {

  if (empty()) {
    return;
  }
    
  const_region_iterator riter = begin();
  while (riter != end()) {
    const_region_iterator piter = riter;
    piter++;
    while (piter != end()) {
      Nodes inter_nodes;
      riter->intersection((*piter),inter_nodes);
      if (!inter_nodes.empty()) {
	Region region;
	region.assignNodes(inter_nodes);
	nextLevel.addRegion(region);
      }

      piter++;
    }

    riter++;
  }
}

void RegionLevel::operator= (RegionLevel const& otherLevel) {
  clear();
  const_region_iterator riter = otherLevel.begin();
  while (riter != otherLevel.end()) {
    rl_regions.push_back(*riter);
    riter++;
  }
}

ostream& operator<<(ostream& os, RegionLevel const& regLevel) {
  const_region_iterator riter = regLevel.begin();
  while (riter != regLevel.end()) {
    os << (*riter);
    riter++;
  }
  return os;
}
