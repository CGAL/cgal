#ifndef CGAL_TETRAHEDRAL_REMESHING_VERTEX_SMOOTH_OPERATION_H
#define CGAL_TETRAHEDRAL_REMESHING_VERTEX_SMOOTH_OPERATION_H

#include <CGAL/license/Tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/internal/elementary_operations.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Tetrahedral_remeshing/internal/FMLS.h>

#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive_3.h>
#include <CGAL/AABB_segment_primitive_3.h>
#include <CGAL/use.h>

#include <vector>
#include <string>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

enum class SmoothingDomain {
  INTERNAL_VERTICES,
  COMPLEX_EDGES,
  SURFACE_VERTICES
};

template<typename C3t3, typename SizingFunction, typename CellSelector, SmoothingDomain Domain>
class VertexSmoothOperation 
  : public ElementaryOperation<C3t3, 
                          typename C3t3::Triangulation::Vertex_handle,
                          typename C3t3::Triangulation::Cell_handle>
{
public:
  using Base = ElementaryOperation<C3t3, 
                              typename C3t3::Triangulation::Vertex_handle,
                              typename C3t3::Triangulation::Cell_handle>;
  using ElementType = typename Base::ElementType;
  using Lock_zone = typename Base::Lock_zone;

  typedef typename C3t3::Triangulation       Tr;
  typedef typename C3t3::Surface_patch_index Surface_patch_index;
  typedef typename Tr::Cell_handle           Cell_handle;
  typedef typename Tr::Vertex_handle         Vertex_handle;
  typedef typename Tr::Edge                  Edge;
  typedef typename Tr::Facet                 Facet;

  typedef typename Tr::Geom_traits           Gt;
  typedef typename Gt::Vector_3              Vector_3;
  typedef typename Gt::Point_3               Point_3;
  typedef typename Gt::FT                    FT;

  using Triangle_vec = std::vector<typename Tr::Triangle>;
  using Triangle_iter = typename Triangle_vec::iterator;
  using Triangle_primitive = CGAL::AABB_triangle_primitive_3<Gt, Triangle_iter>;
  using AABB_triangle_traits = CGAL::AABB_traits_3<Gt, Triangle_primitive>;
  using AABB_triangle_tree = CGAL::AABB_tree<AABB_triangle_traits>;

  using Segment_vec = std::vector<typename Gt::Segment_3>;
  using Segment_iter = typename Segment_vec::iterator;
  using Segment_primitive = CGAL::AABB_segment_primitive_3<Gt, Segment_iter>;
  using AABB_segment_traits = CGAL::AABB_traits_3<Gt, Segment_primitive>;
  using AABB_segment_tree = CGAL::AABB_tree<AABB_segment_traits>;

  VertexSmoothOperation(C3t3& c3t3,
                       const SizingFunction& sizing,
                       const CellSelector& cell_selector,
                       const bool protect_boundaries,
                       const bool smooth_constrained_edges)
    : m_c3t3(c3t3)
    , m_sizing(sizing)
    , m_cell_selector(cell_selector)
    , m_protect_boundaries(protect_boundaries)
    , m_smooth_constrained_edges(smooth_constrained_edges)
  {}

  virtual bool should_process_element(const ElementType& v, const C3t3& c3t3) const override {
    if constexpr (Domain == SmoothingDomain::SURFACE_VERTICES) {
      return false;
    } else if constexpr (Domain == SmoothingDomain::COMPLEX_EDGES) {
      return false;
    } else if constexpr (Domain == SmoothingDomain::INTERNAL_VERTICES) {
      return false;
    } else {
      return false;
    }
  }

  std::vector<ElementType> get_element_source(const C3t3& c3t3) const override {
    std::vector<ElementType> vertices;
    
    if constexpr (Domain == SmoothingDomain::SURFACE_VERTICES) {
    } else if constexpr (Domain == SmoothingDomain::COMPLEX_EDGES) {
    } else if constexpr (Domain == SmoothingDomain::INTERNAL_VERTICES) {
    } else {
    }
    std::cerr<<"VertexSmoothOperation::get_element_source not implemented"<<std::endl;
    return c3t3.triangulation().finite_vertices();
  }

    return vertices;
  }

  bool can_apply_operation(const ElementType& v, const C3t3& c3t3) const override {
    if constexpr (Domain == SmoothingDomain::SURFACE_VERTICES) {
      return true;
    } else if constexpr (Domain == SmoothingDomain::COMPLEX_EDGES) {
      return true;
    } else if constexpr (Domain == SmoothingDomain::INTERNAL_VERTICES) {
      return true;
    } else {
      return true;
    }
  }

  Lock_zone get_lock_zone(const ElementType& v, const C3t3& c3t3) const override {
    Lock_zone zone;
      return Lock_zone();
    } else if constexpr (Domain == SmoothingDomain::COMPLEX_EDGES) {
      return Lock_zone();
    } else if constexpr (Domain == SmoothingDomain::INTERNAL_VERTICES) {
      return Lock_zone();
    } else {
      return Lock_zone();
    }
    
    return zone;
  }

  bool execute_pre_operation(const ElementType& v, C3t3& c3t3) override {
      std::cout<<"VertexSmoothOperation::execute_pre_operation"<<std::endl;
    if constexpr (Domain == SmoothingDomain::SURFACE_VERTICES) {
      return true;
    } else if constexpr (Domain == SmoothingDomain::COMPLEX_EDGES) {
      return true;
    } else if constexpr (Domain == SmoothingDomain::INTERNAL_VERTICES) {
      return true;
    } else {
      return true;
    }
  }

  bool execute_operation(const ElementType& v, C3t3& c3t3) override {
      std::cout<<"VertexSmoothOperation::execute_operation"<<std::endl;
    if constexpr (Domain == SmoothingDomain::SURFACE_VERTICES) {
      return true;
    } else if constexpr (Domain == SmoothingDomain::COMPLEX_EDGES) {
      return true;
    } else if constexpr (Domain == SmoothingDomain::INTERNAL_VERTICES) {
      return true;
    } else {
      return true;
    }
  }

  bool execute_post_operation(const ElementType& v, C3t3& c3t3) override {
      std::cout<<"VertexSmoothOperation::execute_post_operation"<<std::endl;
    if constexpr (Domain == SmoothingDomain::SURFACE_VERTICES) {
      return true;
    } else if constexpr (Domain == SmoothingDomain::COMPLEX_EDGES) {
      return true;
    } else if constexpr (Domain == SmoothingDomain::INTERNAL_VERTICES) {
      return true;
    } else {
      return true;
    }
  }

  std::string operation_name() const override {
    if constexpr (Domain == SmoothingDomain::SURFACE_VERTICES) {
      return "Vertex Smooth (Surface Vertices)";
    } else if constexpr (Domain == SmoothingDomain::COMPLEX_EDGES) {
      return "Vertex Smooth (Complex Edges) - Not implemented";
    } else if constexpr (Domain == SmoothingDomain::INTERNAL_VERTICES) {
      return "Vertex Smooth (Internal Vertices) - Not implemented";
    } else {
      return "Vertex Smooth (Unknown Domain)";
    }
  }

private:

  void collect_vertices_surface_indices(const C3t3& c3t3) {
    // ... implementation ...
  }

  void compute_vertices_normals(const C3t3& c3t3) {
    // ... implementation ...
  }

  void build_aabb_trees(const C3t3& c3t3) {
    // ... implementation ...
  }

  void collect_surface_vertex_moves(const C3t3& c3t3) {
    // ... implementation ...
  }

  std::vector<Surface_patch_index> get_surface_indices(const Vertex_handle& v, const C3t3& c3t3) const {
    // ... implementation ...
  }

  std::size_t vertex_id(const Vertex_handle& v) const {
    // ... implementation ...
  }

  Point_3 project_on_tangent_plane(const Point_3& p, const Point_3& origin, const Vector_3& normal) const {
    // ... implementation ...
  }

  C3t3& m_c3t3;
  const SizingFunction& m_sizing;
  const CellSelector& m_cell_selector;
  bool m_protect_boundaries;
  bool m_smooth_constrained_edges;
};

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_VERTEX_SMOOTH_OPERATION_H