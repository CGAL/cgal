#ifndef EDGE_LENGTH_H
#define EDGE_LENGTH_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/property_map.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Number_type = Kernel::FT;

template <typename Arrangement> struct Edge_length {
  // Boost property-type definitions.
  using category = boost::readable_property_map_tag;
  using value_type = Number_type;
  using reference = value_type;
  using key_type = typename Arrangement::Halfedge_handle;

  value_type operator()(typename Arrangement::Halfedge_handle e) const {
    const auto diff_x = e->target()->point().x() - e->source()->point().x();
    const auto diff_y = e->target()->point().y() - e->source()->point().y();
    return diff_x * diff_x + diff_y * diff_y;
  }

  friend value_type get(const Edge_length& edge_length, key_type key)
  { return edge_length(key); }
};

#endif
