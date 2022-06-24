
#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Manhattan_distance_iso_box_point.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Incremental_neighbor_search.h>


using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_2;
using TreeTraits = CGAL::Search_traits_2<K>;
using Splitter = CGAL::Sliding_midpoint<TreeTraits>;
using Distance = CGAL::Manhattan_distance_iso_box_point<TreeTraits>;

using IncrNN = CGAL::Incremental_neighbor_search<TreeTraits, Distance>;
using Tree = IncrNN::Tree;

int main() {

  Tree tree;
  tree.insert({1,1});
  tree.insert({2,2});

  Point pQuery(0, 0);

  IncrNN nn(tree, {pQuery, pQuery});

  for (IncrNN::iterator it = nn.begin(); it != nn.end(); ++it) {
    std::cout << it->first << " dist: " << it->second   << std::endl;;
  }

  return 0;

}
