
#ifndef OCTREE_SPLIT_CRITERION_H
#define OCTREE_SPLIT_CRITERION_H

#include <iostream>

namespace CGAL {

  struct Split_to_max_depth_or_bucket_size {

    std::size_t m_max_depth, m_bucket_size;

    Split_to_max_depth_or_bucket_size(std::size_t max_depth, std::size_t bucket_size) :
            m_max_depth(max_depth), m_bucket_size(bucket_size) {}

    template<class Node>
    bool operator()(const Node &n) const {
      return (n.num_points() > m_bucket_size && n.depth() < m_max_depth);
    }
  };

  struct Split_to_max_depth {

    std::size_t m_max_depth;

    Split_to_max_depth(std::size_t max_depth) : m_max_depth(max_depth) {}

    template<class Node>
    bool operator()(const Node &n) const {
      return n.depth() < m_max_depth;
    }
  };
}

/*

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

#endif //OCTREE_SPLIT_CRITERION_H
