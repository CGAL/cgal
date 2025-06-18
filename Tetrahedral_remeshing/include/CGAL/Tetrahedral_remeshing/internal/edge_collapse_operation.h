#ifndef CGAL_TETRAHEDRAL_REMESHING_EDGE_COLLAPSE_OPERATION_H
#define CGAL_TETRAHEDRAL_REMESHING_EDGE_COLLAPSE_OPERATION_H

#include <CGAL/license/Tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/internal/elementary_operations.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Iterator_range.h>

#include <boost/container/flat_set.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/bimap.hpp>

#include <utility>
#include <optional>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <iterator>

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

template<typename C3t3, typename SizingFunction, typename CellSelector>
class EdgeCollapseOperation 
  : public ElementaryOperation<C3t3, 
                          typename C3t3::Triangulation::Finite_edges_iterator,
                          typename C3t3::Triangulation::Cell_handle>
{
public:
  using Base = ElementaryOperation<C3t3, 
                              typename C3t3::Triangulation::Finite_edges_iterator,
                              typename C3t3::Triangulation::Cell_handle>;
  using Triangulation = typename C3t3::Triangulation;
  using Edge = typename Triangulation::Edge;
  using ElementIteratorType= typename Base::ElementIteratorType;
  using Cell_handle = typename Triangulation::Cell_handle;
  using Vertex_handle = typename Triangulation::Vertex_handle;
  using Point = typename Triangulation::Point;
  using Lock_zone = typename Base::Lock_zone;
  using Facet = typename Triangulation::Facet;
  using Subdomain_index = typename C3t3::Subdomain_index;
  using Surface_patch_index = typename C3t3::Surface_patch_index;
  using Curve_index = typename C3t3::Curve_index;
  using EdgeIteratorType= typename C3t3::Triangulation::Finite_edges_iterator;

  EdgeCollapseOperation(C3t3& c3t3,
                       const SizingFunction& sizing,
                       const CellSelector& cell_selector,
                       const bool protect_boundaries)
    : m_c3t3(c3t3)
    , m_sizing(sizing)
    , m_cell_selector(cell_selector)
    , m_protect_boundaries(protect_boundaries)
  {}

  bool should_process_element(const EdgeIteratorType& e, const C3t3& c3t3) const override {
    std::cout<<"EdgeCollapseOperation::should_process_element"<<std::endl;
    return true;
  }

  typename CGAL::Iterator_range<EdgeIteratorType> get_element_iterators(const C3t3& c3t3) const override {
    return c3t3.triangulation().finite_edges();
  }

  bool can_apply_operation(const EdgeIteratorType& e, const C3t3& c3t3) const override {
    std::cout<<"EdgeCollapseOperation::can_apply_operation"<<std::endl;
    return true;
  }

  Lock_zone get_lock_zone(const EdgeIteratorType& e, const C3t3& c3t3) const override {
    Lock_zone zone;
    std::cout<<"EdgeCollapseOperation::get_lock_zone"<<std::endl;
    return zone;
  }

  bool execute_pre_operation(const EdgeIteratorType& e, C3t3& c3t3) override {
    std::cout<<"EdgeCollapseOperation::execute_pre_operation"<<std::endl;
    return true;
  }

  bool execute_operation(const EdgeIteratorType& e, C3t3& c3t3) override {
    std::cout<<"EdgeCollapseOperation::execute_operation"<<std::endl;
    return true;
  }

  bool execute_post_operation(const EdgeIteratorType& e, C3t3& c3t3) override {
    std::cout<<"EdgeCollapseOperation::execute_post_operation"<<std::endl;
    return true;
  }

  std::string operation_name() const override {
    return "Edge Collapse";
  }

private:
  C3t3& m_c3t3;
  const SizingFunction& m_sizing;
  const CellSelector& m_cell_selector;
  bool m_protect_boundaries;
};

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_EDGE_COLLAPSE_OPERATION_H 