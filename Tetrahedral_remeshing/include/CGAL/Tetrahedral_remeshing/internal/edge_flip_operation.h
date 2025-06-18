#ifndef CGAL_TETRAHEDRAL_REMESHING_EDGE_FLIP_OPERATION_H
#define CGAL_TETRAHEDRAL_REMESHING_EDGE_FLIP_OPERATION_H

#include <CGAL/license/Tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/internal/elementary_operations.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Triangulation_3.h>

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

template<typename C3t3, typename CellSelector>
class EdgeFlipOperation 
  : public ElementaryOperation<C3t3, 
                          typename C3t3::Triangulation::Finite_edges_iterator,
                          typename C3t3::Triangulation::Cell_handle> {
public:
  using Base = ElementaryOperation<C3t3, 
                              typename C3t3::Triangulation::Finite_edges_iterator,
                              typename C3t3::Triangulation::Cell_handle>;
  using Complex = C3t3;
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
  using EdgeIteratorType= typename C3t3::Triangulation::Finite_edges_iterator;

  EdgeFlipOperation(C3t3& c3t3,
                   const CellSelector& cell_selector,
                   const bool protect_boundaries)
    : m_c3t3(c3t3)
    , m_cell_selector(cell_selector)
    , m_protect_boundaries(protect_boundaries)
  {}

  bool should_process_element(const Edge& e, const C3t3& c3t3) const override {
    std::cout<<"EdgeFlipOperation::should_process_element"<<std::endl;
    return true;
  }

  typename CGAL::Iterator_range<EdgeIteratorType> get_element_iterators(const C3t3& c3t3) const override {
    return c3t3.triangulation().finite_edges();
  }

  bool can_apply_operation(const EdgeIteratorType& e, const C3t3& c3t3) const override {
    std::cout<<"EdgeFlipOperation::can_apply_operation"<<std::endl;
    return true;
  }

  Lock_zone get_lock_zone(const EdgeIteratorType& e, const C3t3& c3t3) const override {
    Lock_zone zone;
    std::cout<<"EdgeFlipOperation::get_lock_zone"<<std::endl;
    return zone;
  }

  bool execute_pre_operation(const EdgeIteratorType& e, C3t3& c3t3) override {
    std::cout<<"EdgeFlipOperation::execute_pre_operation"<<std::endl;
    return true;
  } 

  bool execute_operation(const EdgeIteratorType& e, C3t3& c3t3) override {
    std::cout<<"EdgeFlipOperation::execute_operation"<<std::endl;
    return true;
  }

  bool execute_post_operation(const EdgeIteratorType& e, C3t3& c3t3) override {
    std::cout<<"EdgeFlipOperation::execute_post_operation"<<std::endl;
    return true;
  }

  std::string operation_name() const override {
    return "Edge Flip";
  }

private:
  C3t3& m_c3t3;
  const CellSelector& m_cell_selector;
  bool m_protect_boundaries;
};

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_EDGE_FLIP_OPERATION_H 