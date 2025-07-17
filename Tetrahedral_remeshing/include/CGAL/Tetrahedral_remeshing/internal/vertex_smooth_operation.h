#ifndef CGAL_TETRAHEDRAL_REMESHING_VERTEX_SMOOTH_OPERATION_H
#define CGAL_TETRAHEDRAL_REMESHING_VERTEX_SMOOTH_OPERATION_H

#include <CGAL/license/Tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/internal/elementary_operations.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Tetrahedral_remeshing/internal/FMLS.h>
#include <CGAL/Tetrahedral_remeshing/internal/vertex_smoothing_context.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/Tetrahedral_remeshing/internal/smooth_vertices.h>

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

// Base class for shared functionality between different vertex smoothing operations
template<typename C3t3, typename SizingFunction, typename CellSelector>
class VertexSmoothOperationBase
{
protected:
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

  // AABB tree types (commented out until needed)
  // using Triangle_vec = std::vector<typename Tr::Triangle>;
  // using Triangle_iter = typename Triangle_vec::iterator;
  // using Triangle_primitive = CGAL::AABB_triangle_primitive_3<Gt, Triangle_iter>;
  // using AABB_triangle_traits = CGAL::AABB_traits_3<Gt, Triangle_primitive>;
  // using AABB_triangle_tree = CGAL::AABB_tree<AABB_triangle_traits>;

  // using Segment_vec = std::vector<typename Gt::Segment_3>;
  // using Segment_iter = typename Segment_vec::iterator;
  // using Segment_primitive = CGAL::AABB_segment_primitive_3<Gt, Segment_iter>;
  // using AABB_segment_traits = CGAL::AABB_traits_3<Gt, Segment_primitive>;
  // using AABB_segment_tree = CGAL::AABB_tree<AABB_segment_traits>;

  using Context = VertexSmoothingContext<C3t3, SizingFunction, CellSelector>;

  const SizingFunction& m_sizing;
  const CellSelector& m_cell_selector;
  bool m_protect_boundaries;
  Context* m_context{nullptr}; // Pointer to shared context

public:
  bool m_smooth_constrained_edges;
  VertexSmoothOperationBase(
                       const SizingFunction& sizing,
                       const CellSelector& cell_selector,
                       const bool protect_boundaries,
                       const bool smooth_constrained_edges,
                       Context* context)
    :
     m_sizing(sizing)
    , m_cell_selector(cell_selector)
    , m_protect_boundaries(protect_boundaries)
    , m_smooth_constrained_edges(smooth_constrained_edges)
    , m_context(context)
  {}

  void set_context(Context* context) {
    m_context = context;
  }

protected:
  // Helper to get incident cells for a vertex
  boost::container::small_vector<Cell_handle, 64> get_incident_cells(Vertex_handle vh, const C3t3& c3t3) const {
    boost::container::small_vector<Cell_handle, 64> inc_cells;
#ifdef USE_THREADSAFE_INCIDENT_CELLS
    c3t3.triangulation().incident_cells_threadsafe(vh, std::back_inserter(inc_cells));
#else
    c3t3.triangulation().incident_cells(vh, std::back_inserter(inc_cells));
#endif
    return inc_cells;
  }

  // Common helper methods for all smoothing operations
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
    return std::vector<Surface_patch_index>();
  }

  std::size_t vertex_id(const Vertex_handle& v) const {
    // ... implementation ...
    return 0;
  }

  Point_3 project_on_tangent_plane(const Point_3& p, const Point_3& origin, const Vector_3& normal) const {
    // ... implementation ...
    return p;
  }
};

// Internal Vertex Smooth Operation - processes internal vertices
template<typename C3t3, typename SizingFunction, typename CellSelector>
class InternalVertexSmoothOperation
  : public VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>,
    public ElementaryOperation<C3t3,
                          typename C3t3::Triangulation::Vertex_handle,
                          CGAL::Iterator_range<typename C3t3::Triangulation::Vertex_handle>,
                          typename C3t3::Triangulation::Vertex_handle>
{
public:
  using BaseClass = VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>;
  using Base = ElementaryOperation<C3t3,
                              typename C3t3::Triangulation::Vertex_handle,
                              CGAL::Iterator_range<typename C3t3::Triangulation::Vertex_handle>,
                              typename C3t3::Triangulation::Vertex_handle>;
  using ElementType = typename C3t3::Triangulation::Vertex_handle;
  using ElementSource = typename Base::ElementSource;
  using Lock_zone = typename Base::Lock_zone;

  using BaseClass::m_sizing;
  using BaseClass::m_cell_selector;
  using BaseClass::m_protect_boundaries;
  using BaseClass::m_smooth_constrained_edges;
  using BaseClass::m_context;
  using BaseClass::get_incident_cells;

  // Import types from base class
  using typename BaseClass::Tr;
  using typename BaseClass::Cell_handle;
  using typename BaseClass::Vertex_handle;
  using typename BaseClass::Point_3;

public:
  InternalVertexSmoothOperation(
                               const SizingFunction& sizing,
                               const CellSelector& cell_selector,
                               const bool protect_boundaries,
                               const bool smooth_constrained_edges,
                               typename BaseClass::Context* context)
    : BaseClass(sizing, cell_selector, protect_boundaries, smooth_constrained_edges, context)
  {}

  bool should_process_element(const ElementType& v, const C3t3& c3t3) const override {
    // TODO: Implement internal vertex filtering logic
    return false;
  }

  ElementSource get_element_source(const C3t3& c3t3) const override {
    return c3t3.triangulation().finite_vertex_handles();
  }

  bool can_apply_operation(const ElementType& v, const C3t3& c3t3) const override {
    // TODO: Validate vertex is still internal and movable
    return true;
  }

  bool lock_zone(const ElementType& v, const C3t3& c3t3) const override {
    // TODO: Implement lock zone for internal vertex smoothing
    return true;
  }

  bool execute_operation(const ElementType& v, C3t3& c3t3) override {
    // TODO: Implement internal vertex smoothing
    std::cout << "InternalVertexSmoothOperation::execute_operation" << std::endl;
    return true;
  }

  std::string operation_name() const override {
    return "Vertex Smooth (Internal Vertices)";
  }
};

// Surface Vertex Smooth Operation - processes surface vertices (2-junctions)
template<typename C3t3, typename SizingFunction, typename CellSelector>
class SurfaceVertexSmoothOperation
  : public VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>,
    public ElementaryOperation<C3t3,
                          typename C3t3::Triangulation::Vertex_handle,
                          typename C3t3::Triangulation::Finite_vertex_handles,
                          typename C3t3::Triangulation::Vertex_handle>
{
public:
  using BaseClass = VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>;
  using Base = ElementaryOperation<C3t3,
                              typename C3t3::Triangulation::Vertex_handle,
							  typename C3t3::Triangulation::Finite_vertex_handles,
                              typename C3t3::Triangulation::Vertex_handle>;
  using ElementType = typename C3t3::Triangulation::Vertex_handle;
  using ElementSource = typename Base::ElementSource;
  using Lock_zone = typename Base::Lock_zone;

  using BaseClass::m_sizing;
  using BaseClass::m_cell_selector;
  using BaseClass::m_protect_boundaries;
  using BaseClass::m_smooth_constrained_edges;
  using BaseClass::m_context;
  using BaseClass::get_incident_cells;

  // Import types from base class
  using typename BaseClass::Tr;
  using typename BaseClass::Cell_handle;
  using typename BaseClass::Vertex_handle;
  using typename BaseClass::Point_3;

private:
  // Surface-specific data (only used by SurfaceVertexSmoothOperation)
  std::unordered_map<Vertex_handle, std::vector<Surface_patch_index>> m_vertices_surface_indices;
  std::unordered_map<Vertex_handle, std::unordered_map<Surface_patch_index, Vector_3, boost::hash<Surface_patch_index>>> m_vertices_normals;

  // Triangle AABB tree (only used by SurfaceVertexSmoothOperation)
  using Triangle_vec = std::vector<typename Tr::Triangle>;
  using Triangle_iter = typename Triangle_vec::iterator;
  using Triangle_primitive = CGAL::AABB_triangle_primitive_3<Gt, Triangle_iter>;
  using AABB_triangle_traits = CGAL::AABB_traits_3<Gt, Triangle_primitive>;
  using AABB_triangle_tree = CGAL::AABB_tree<AABB_triangle_traits>;
  
  Triangle_vec m_aabb_triangles;
  AABB_triangle_tree m_triangles_aabb_tree;
  FT m_aabb_epsilon;

  void build_triangle_aabb_tree(const C3t3& c3t3) {
    for (const Facet& f : c3t3.facets_in_complex()) {
      m_aabb_triangles.push_back(c3t3.triangulation().triangle(f));
    }
    m_triangles_aabb_tree.rebuild(m_aabb_triangles.begin(), m_aabb_triangles.end());
    m_triangles_aabb_tree.accelerate_distance_queries();

    if (!m_triangles_aabb_tree.empty()) {
        const CGAL::Bbox_3& bb = m_triangles_aabb_tree.bbox();
        m_aabb_epsilon = 1e-3 * (std::min)(bb.xmax() - bb.xmin(),
                                (std::min)(bb.ymax() - bb.ymin(),
                                           bb.zmax() - bb.zmin()));
    }
  }

  void collect_vertices_surface_indices(const C3t3& c3t3) {
      for (const Facet& fit : c3t3.facets_in_complex()) {
          const Surface_patch_index& surface_index = c3t3.surface_patch_index(fit);
          for (const Vertex_handle vi : c3t3.triangulation().vertices(fit)) {
              auto& v_surface_indices = m_vertices_surface_indices[vi];
              if (std::find(v_surface_indices.begin(), v_surface_indices.end(), surface_index) == v_surface_indices.end())
                  v_surface_indices.push_back(surface_index);
          }
      }
  }

  void compute_vertices_normals(const C3t3& c3t3) {
      // This is a simplified version. A full implementation would be more complex.
      // For now, we are just ensuring the data structure is populated.
      for (auto const& [vertex, surfaces] : m_vertices_surface_indices) {
          for (const auto& surface_index : surfaces) {
              m_vertices_normals[vertex][surface_index] = Vector_3(0,0,1); // Dummy normal
          }
      }
  }

public:
  SurfaceVertexSmoothOperation(const C3t3& c3t3,
                              const SizingFunction& sizing,
                              const CellSelector& cell_selector,
                              const bool protect_boundaries,
                              const bool smooth_constrained_edges,
                              typename BaseClass::Context* context)
    : BaseClass( sizing, cell_selector, protect_boundaries, smooth_constrained_edges, context)
  {
    build_triangle_aabb_tree(c3t3);
    if (!m_protect_boundaries) {
        collect_vertices_surface_indices(c3t3);
        compute_vertices_normals(c3t3);
    }
  }

  bool should_process_element(const ElementType& v, const C3t3& c3t3) const override {
    // TODO: Implement surface vertex filtering logic
    return false;
  }

  ElementSource get_element_source(const C3t3& c3t3) const override {
    return c3t3.triangulation().finite_vertex_handles();
  }

  bool can_apply_operation(const ElementType& v, const C3t3& c3t3) const override {
    // TODO: Validate vertex is still on surface and movable
    return true;
  }

  bool lock_zone(const ElementType& v, const C3t3& c3t3) const override {
    // TODO: Implement lock zone for surface vertex smoothing
    return true;
  }

  bool execute_operation(const ElementType& v, C3t3& c3t3) override {
    // TODO: Implement surface vertex smoothing
    std::cout << "SurfaceVertexSmoothOperation::execute_operation" << std::endl;
    return true;
  }

  std::string operation_name() const override {
    return "Vertex Smooth (Surface Vertices)";
  }
};

// Complex Edge Vertex Smooth Operation - processes vertices on complex edges (1-junctions)
template<typename C3t3, typename SizingFunction, typename CellSelector>
class ComplexEdgeVertexSmoothOperation
  : public VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>,
    public ElementaryOperation<C3t3,
                          typename C3t3::Triangulation::Vertex_handle,
                          typename C3t3::Triangulation::Finite_vertex_handles,
                          typename C3t3::Triangulation::Vertex_handle>
{
public:
  bool m_flip_smooth_steps = false; //TODO: when the smooth-steps start this needs to be true
  using BaseClass = VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>;
  using Base = ElementaryOperation<C3t3,
                              typename C3t3::Triangulation::Vertex_handle,
                              typename C3t3::Triangulation::Finite_vertex_handles,
                              typename C3t3::Triangulation::Vertex_handle>;
  using ElementType = typename C3t3::Triangulation::Vertex_handle;
  using ElementSource = typename Base::ElementSource;
  using Lock_zone = typename Base::Lock_zone;

  using BaseClass::m_sizing;
  using BaseClass::m_cell_selector;
  using BaseClass::m_protect_boundaries;
  using BaseClass::m_smooth_constrained_edges;
  using BaseClass::m_context;
  using BaseClass::get_incident_cells;

  // Import types from base class
  using typename BaseClass::Tr;
  using typename BaseClass::Cell_handle;
  using typename BaseClass::Vertex_handle;
  using typename BaseClass::Point_3;

private:
  // Segment AABB tree (only used by ComplexEdgeVertexSmoothOperation)
  using Segment_vec = std::vector<typename Gt::Segment_3>;
  using Segment_iter = typename Segment_vec::iterator;
  using Segment_primitive = CGAL::AABB_segment_primitive_3<Gt, Segment_iter>;
  using AABB_segment_traits = CGAL::AABB_traits_3<Gt, Segment_primitive>;
  using AABB_segment_tree = CGAL::AABB_tree<AABB_segment_traits>;

  Segment_vec m_aabb_segments;
  AABB_segment_tree m_segments_aabb_tree;

  void build_segment_aabb_tree(const C3t3& c3t3) {
    std::ofstream debug_log("debug_refactored_aabb_construction.log");
    debug_log << "=== REFACTORED AABB TREE CONSTRUCTION ===" << std::endl;
    
    std::size_t segment_count = 0;
    for (const Edge& e : c3t3.edges_in_complex()) {
      auto segment = c3t3.triangulation().segment(e);
      m_aabb_segments.push_back(segment);
      
      debug_log << std::setprecision(25) << "Segment " << segment_count 
                << ": " << segment.source() << " -> " << segment.target() << std::endl;
      segment_count++;
    }
    
    debug_log << "Total segments: " << segment_count << std::endl;
    debug_log.close();
    
    m_segments_aabb_tree.rebuild(m_aabb_segments.begin(), m_aabb_segments.end());
    m_segments_aabb_tree.accelerate_distance_queries();
  }

  Point_3 point(const Point_3& p) const {
    return p;
  }

  FT density_along_segment(const Edge& e, const C3t3& c3t3, bool boundary_edge = false) const {
    // Simplified density calculation - in practice this would use sizing function
    const auto [pt, dim, index] = midpoint_with_info(e, boundary_edge, c3t3);
    const FT s = sizing_at_midpoint(e, pt, dim, index, m_sizing, c3t3, m_cell_selector);
    const FT density = 1. / s; //density = 1 / size^(dimension)
    return density;
  }

  template<typename CellRange>
  Dihedral_angle_cosine max_cosine(const Tr& tr,
    const CellRange& cells) const
  {
    Dihedral_angle_cosine max_cos_dh = cosine_of_90_degrees();// = 0.
    for (Cell_handle c : cells)
    {
      if(!is_selected(c))
        continue;
      Dihedral_angle_cosine cos_dh = max_cos_dihedral_angle(tr, c, false);
      if (max_cos_dh < cos_dh)
        max_cos_dh = cos_dh;
    }
    return max_cos_dh;
  }

bool is_selected(const Cell_handle c) const { return get(m_cell_selector, c); }

  template <typename CellRange, typename Tr>
  bool check_inversion_and_move(const typename Tr::Vertex_handle v,
                                const typename Tr::Geom_traits::Point_3& final_pos,
                                const CellRange& inc_cells,
                                const Tr& tr,
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
                                FT& total_move) const
#else
                                FT&) const
#endif
  {
    const typename Tr::Point backup = v->point(); // backup v's position
    const typename Tr::Geom_traits::Point_3 pv = point(backup);

    bool valid_orientation = false;
    bool angles_improved = true;
    double frac = 1.0;
    typename Tr::Geom_traits::Vector_3 move(pv, final_pos);

    // TODO: This variable is only used in the smooth-steps, so it should be removed
    const Dihedral_angle_cosine curr_max_cos = m_flip_smooth_steps
                                                   ? max_cosine(tr, inc_cells)
                                                   : Dihedral_angle_cosine(CGAL::ZERO, 0., 1.); // Dummy unused value

    bool valid_try = true;
    do {
      v->set_point(typename Tr::Point(pv + frac * move));

      valid_try = true;
      valid_orientation = true;
      angles_improved = true;

      for(const typename Tr::Cell_handle& ci : inc_cells) {
        if(CGAL::POSITIVE != CGAL::orientation(point(ci->vertex(0)->point()), point(ci->vertex(1)->point()),
                                               point(ci->vertex(2)->point()), point(ci->vertex(3)->point())))
        {
          frac = 0.5 * frac;
          valid_try = false;
          valid_orientation = false;
          break;
        } else if(m_flip_smooth_steps) // check that dihedral angles get improved
        {
          if(is_selected(ci)) {
            Dihedral_angle_cosine max_cos_ci = max_cos_dihedral_angle(tr, ci, false);
            if(curr_max_cos < max_cos_ci) {
              // keep move only if new cosine is smaller than previous one
              // i.e. if angle is larger
              frac = 0.5 * frac;
              valid_try = false;
              angles_improved = false;
              break;
            }
          }
        }
      }
    } while(!valid_try && frac > 0.1);

    // if move failed, cancel move
    bool valid_move = valid_orientation && angles_improved;

    if(!valid_move)
      v->set_point(backup);

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    else
      total_move += CGAL::approximate_sqrt(CGAL::squared_distance(pv, point(v->point())));
#endif

    return valid_move;
  }

  void collect_vertices_surface_indices(
      const C3t3& c3t3, std::unordered_map<Vertex_handle, std::vector<Surface_patch_index>>& vertices_surface_indices) {
    for(Facet fit : c3t3.facets_in_complex()) {
      const Surface_patch_index& surface_index = c3t3.surface_patch_index(fit);

      for(const Vertex_handle vi : c3t3.triangulation().vertices(fit)) {
        std::vector<Surface_patch_index>& v_surface_indices = vertices_surface_indices[vi];
        if(std::find(v_surface_indices.begin(), v_surface_indices.end(), surface_index) == v_surface_indices.end())
          v_surface_indices.push_back(surface_index);
      }
    }
  }

public:
    FT total_move = 0;
    static std::ofstream debug_log; // Static member for log file

  ComplexEdgeVertexSmoothOperation(const C3t3& c3t3,
                                  const SizingFunction& sizing,
                                  const CellSelector& cell_selector,
                                  const bool protect_boundaries,
                                  typename BaseClass::Context* context)
    : BaseClass( sizing, cell_selector, protect_boundaries, true, context)
  {
    build_segment_aabb_tree(c3t3);
  }

  bool should_process_element(const ElementType& v, const C3t3& c3t3) const override {
    // Only process vertices that are on features (complex edges)
    return can_apply_operation(v, c3t3);
  }

  ElementSource get_element_source(const C3t3& c3t3) const override {
    perform_global_preprocessing(c3t3);
    return c3t3.triangulation().finite_vertex_handles();
  }

  bool can_apply_operation(const ElementType& v, const C3t3& c3t3) const override {
    // Check if vertex is still on feature and can move
    const std::size_t vid = m_context->m_vertex_id.at(v);
    return m_context->m_free_vertices[vid] && is_on_feature(v);
  }

  bool lock_zone(const ElementType& v, const C3t3& c3t3) const override {
    // TODO: Implement lock zone for complex edge vertex smoothing
    return true;
  }

  bool execute_operation(const ElementType& v, C3t3& c3t3) override {
    auto& tr = c3t3.triangulation();
    const std::size_t vid = m_context->m_vertex_id.at(v);
    
    if (!m_context->m_free_vertices[vid] || !is_on_feature(v))
      return false;

    const Point_3 current_pos = point(v->point());
    const auto& moves = m_context->m_moves;
    
    if (moves[vid].neighbors == 0)
      return false;

    CGAL_assertion(moves[vid].mass > 0);
    const Vector_3 move = moves[vid].move / moves[vid].mass;
    
    // Open log file for appending
    std::ofstream debug_log("debug_refactored_smooth_edges.log", std::ios::app);
    debug_log << std::setprecision(25);
    debug_log << "DEBUG REFACTORED: Vertex " << vid << " - original_position: " << current_pos << std::endl;
    debug_log << "DEBUG REFACTORED: Processing vertex " << vid 
              << " - neighbors=" << moves[vid].neighbors << ", mass=" << moves[vid].mass 
              << ", total_move_vector=" << moves[vid].move 
              << ", computed_move=" << move << std::endl;
    
    const Point_3 smoothed_position = current_pos + move;

    debug_log << std::setprecision(25) << "DEBUG REFACTORED: Vertex " << vid << " - smoothed_position: " << smoothed_position << std::endl;

    //TODO: Test #ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS

      Vector_3 sum_projections = CGAL::NULL_VECTOR;
      Point_3 tmp_pos = current_pos;

#ifndef CGAL_TET_REMESHING_EDGE_SMOOTHING_DISABLE_PROJECTION
      const std::vector<Surface_patch_index>& v_surface_indices = m_vertices_surface_indices.at(v);
      for (const Surface_patch_index& si : v_surface_indices)
      {
        Point_3 normal_projection = project_on_tangent_plane(smoothed_position,
                                                             current_pos,
                                                             vertices_normals.at(v).at(si));

        sum_projections += Vector_3(tmp_pos, normal_projection);
        tmp_pos = normal_projection;
      }
#endif //CGAL_TET_REMESHING_EDGE_SMOOTHING_DISABLE_PROJECTION

      const Point_3 new_pos = current_pos + sum_projections;

#else // AABB_tree projection
    // Use AABB tree projection for complex edges
    const Point_3 new_pos = m_segments_aabb_tree.closest_point(smoothed_position);
#endif

    debug_log << std::setprecision(25) << "DEBUG REFACTORED: Vertex " << vid << " - new_pos: " << new_pos << std::endl;

    const auto& inc_cells = m_context->m_inc_cells[vid];
    
    // Debug: Track move before and after
    FT move_before = total_move;
    
    if (check_inversion_and_move(v, new_pos, inc_cells, tr,total_move)) {
      debug_log << "DEBUG REFACTORED: Vertex " << vid << " moved, total_move: " 
                << move_before << " -> " << total_move << " (delta: " 
                << (total_move - move_before) << ")" << std::endl;
      debug_log.close();
      return true;
    }

    debug_log.close();
    return false;
  }

  std::string operation_name() const override {
    return "Vertex Smooth (Complex Edge Vertices)";
  }

  // Debug method to print final total_move
  void print_final_total_move(const C3t3& c3t3) const {
    std::ofstream debug_log("debug_refactored_smooth_edges.log", std::ios::app);
    debug_log << "DEBUG REFACTORED: Final total_move: " << total_move << std::endl;
    debug_log.close();
  }

  // This method will be called by the execution framework to collect moves
  void perform_global_preprocessing(const C3t3& c3t3) const {
    auto& tr = c3t3.triangulation();
    const std::size_t nbv = tr.number_of_vertices();
    
    // Initialize moves vector
    m_context->m_moves.assign(nbv, typename Context::Move{CGAL::NULL_VECTOR, 0, 0.});

    // Debug: Track edge processing
    std::size_t edges_processed = 0;
    std::size_t vertices_with_moves = 0;

    // Open debug log file
    std::ofstream debug_log("debug_refactored_smooth_edges.log");
    debug_log << "=== REFACTORED SMOOTH EDGES IN COMPLEX DEBUG LOG ===" << std::endl;

    // Collect moves from complex edges
    for (const Edge& e : c3t3.edges_in_complex()) {
      const Vertex_handle vh0 = e.first->vertex(e.second);
      const Vertex_handle vh1 = e.first->vertex(e.third);

      CGAL_expensive_assertion(is_on_feature(vh0));
      CGAL_expensive_assertion(is_on_feature(vh1));

      const std::size_t& i0 = m_context->m_vertex_id.at(vh0);
      const std::size_t& i1 = m_context->m_vertex_id.at(vh1);

      const bool vh0_moving = m_context->m_free_vertices[i0];
      const bool vh1_moving = m_context->m_free_vertices[i1];

      if (!vh0_moving && !vh1_moving)
        continue;

      edges_processed++;
      const Point_3& p0 = point(vh0->point());
      const Point_3& p1 = point(vh1->point());
      const FT density = density_along_segment(e, c3t3, true);

      debug_log << "DEBUG REFACTORED: Edge " << edges_processed << " - vh0=" << i0 
                << " (moving=" << vh0_moving << "), vh1=" << i1 
                << " (moving=" << vh1_moving << "), density=" << density << std::endl;

      if (vh0_moving) {
        Vector_3 move_vector = density * Vector_3(p0, p1);
        m_context->m_moves[i0].move += move_vector;
        m_context->m_moves[i0].mass += density;
        ++m_context->m_moves[i0].neighbors;
        vertices_with_moves++;
        
        debug_log << "DEBUG REFACTORED:   vh0[" << i0 << "] - move_vector=" << move_vector 
                  << ", total_move=" << m_context->m_moves[i0].move << ", mass=" << m_context->m_moves[i0].mass 
                  << ", neighbors=" << m_context->m_moves[i0].neighbors << std::endl;
      }
      if (vh1_moving) {
        Vector_3 move_vector = density * Vector_3(p1, p0);
        m_context->m_moves[i1].move += move_vector;
        m_context->m_moves[i1].mass += density;
        ++m_context->m_moves[i1].neighbors;
        vertices_with_moves++;
        
        debug_log << "DEBUG REFACTORED:   vh1[" << i1 << "] - move_vector=" << move_vector 
                  << ", total_move=" << m_context->m_moves[i1].move << ", mass=" << m_context->m_moves[i1].mass 
                  << ", neighbors=" << m_context->m_moves[i1].neighbors << std::endl;
      }
    }
    
    debug_log << "DEBUG REFACTORED: Processed " << edges_processed << " edges, " 
              << vertices_with_moves << " vertices got moves" << std::endl;
    debug_log.close();
  }
};

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_VERTEX_SMOOTH_OPERATION_H