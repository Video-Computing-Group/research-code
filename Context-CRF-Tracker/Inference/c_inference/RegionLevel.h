#include "Region.h"
#include <vector>

#ifndef __REGION_LEVEL__
#define __REGION_LEVEL__

class RegionLevel {

  /**
     This class defines a level in the region-graph (used for GBP), as a layer
     of regions, where none of them is a subset of another
   
     Part of the c_inference package
     @version November 2004
     @author Talya Meltzer
  */
  
 public:

  // ctor
  RegionLevel() { rl_regions.clear(); }

  // accessors

  unsigned int size() const { return rl_regions.size(); }
  bool empty() const { return rl_regions.empty(); }
  
  region_iterator begin() { return rl_regions.begin(); }
  region_iterator end() { return rl_regions.end(); }
  const_region_iterator begin() const { return rl_regions.begin(); }
  const_region_iterator end() const { return rl_regions.end(); }

  Region operator[](int i) const { return rl_regions[i]; }
  Region& operator[](int i) { return rl_regions[i]; }

  vector<Region>& regions() { return rl_regions; }

  // create a new RegionLevel by intersections of regions in this RegionLevel
  void intersections(RegionLevel& nextLevel) const;

  // mutators

  // region will not be added if it's a subset of an existing one
  bool addRegion(Region& region);

  // removes the element at position pos, returns iterator to the next element
  region_iterator removeRegion(region_iterator pos);

  // clear all regions
  void clear() { rl_regions.clear(); }

  void operator= (RegionLevel const& otherLevel);
  friend ostream& operator<<(ostream& os, RegionLevel const& regLevel);

 private:
  vector<Region> rl_regions;
};

#endif
