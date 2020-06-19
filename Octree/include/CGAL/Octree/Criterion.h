
#ifndef OCTREE_CRITERION_H
#define OCTREE_CRITERION_H

#include <CGAL/Octree/Octree_node.h>

/*

typedef Octree_node<Kernel, PointRange> Node;

// Possible criterions
struct Stop_at_max_depth {
  std::size_t max_depth;

  Stop_at_max_depth(const std::size_t &max_depth) : max_depth(max_depth) {}

  bool operator()(const Node &n) const {
    return n.depth() == max_depth; // not sure you can know that from node only,
    // otherwise your criterion could also take a
    // reference to the full octree as parameter
  }
};

struct Stop_at_max_number_of_points {
  std::size_t max_nb_points;

  Stop_at_max_number_of_points(const std::size_t &max_nb_points)
          : max_nb_points(max_nb_points) {}

  bool operator()(const Node &n) const {
    return n.number_of_points() // internally, the node can use std::distance(begin, end)
           < max_nb_points;
  }
};

// Just for an example using outside info (like a normal map)
// The normals remain unknown to the octree but are used for construction
struct Stop_at_normal_deviation {
  Normal_map normal_map;
  FT max_dev;

  Stop_at_normal_deviation(Normal_map normal_map,
                           FT max_dev)
          : normal_map(normal_map), max_dev(max_dev) {}

  bool operator()(const Node &n) const {
    FT dev = 0;

    for (Iterator it = n.begin(); it != n.end(); ++it)
      dev += compute_deviation(get(normal_map, *it)); // whatever compute_deviation is :)

    // if your node defines begin() and end(), you can also use a C++11 loop:
    // for (const Range_type& r : n)
    //   dev += compute_deviation(get (normal_map, r));

    return dev < max_dev;
  }
};

*/

#endif //OCTREE_CRITERION_H
