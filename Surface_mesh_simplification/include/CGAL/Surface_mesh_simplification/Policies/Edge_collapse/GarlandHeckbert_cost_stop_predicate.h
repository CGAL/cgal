#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_COST_STOP_PREDICATE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_COST_STOP_PREDICATE_H

#include <CGAL/license/Surface_mesh_simplification.h>
#include <CGAL/squared_distance_3.h>
#include <iostream>

namespace CGAL {
namespace Surface_mesh_simplification {

template <class FT>
class GarlandHeckbert_cost_stop_predicate
{
  FT m_gh_cost_threshold;

public:
  GarlandHeckbert_cost_stop_predicate(FT gh_cost_threshold)
    : m_gh_cost_threshold(gh_cost_threshold)
  {}

  template <typename F, typename Profile>
  bool operator()(const F& aCurrentCost,
                  const Profile&,
                  std::size_t,
                  std::size_t) const
  {
    return aCurrentCost >= m_gh_cost_threshold;
  }
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_COST_STOP_PREDICATE_H
