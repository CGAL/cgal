// Copyright (c) 2013  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Xiang Gao <gaox@ethz.ch>
//

#ifndef CGAL_MEAN_CURVATURE_SKELETON_H
#define CGAL_MEAN_CURVATURE_SKELETON_H

#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <CGAL/Default.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>

#include <boost/iterator/transform_iterator.hpp>

#include <CGAL/boost/graph/iterator.h>

// Compute cotangent Laplacian
#include <CGAL/internal/Mean_curvature_skeleton/Weights.h>

// Compute the vertex normal
#include <CGAL/internal/Mean_curvature_skeleton/get_normal.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/internal/Mean_curvature_skeleton/Edge_minimum_length_stop_predicate.h>

// Non-default cost and placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>

// Skip the geometric test
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Geometric_test_skipper.h>

// Curve skeleton data structure
#include <CGAL/internal/Mean_curvature_skeleton/Curve_skeleton.h>

// For Voronoi diagram
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

// For debugging macro
#include <CGAL/internal/Mean_curvature_skeleton/Debug.h>

// Some helper functions
#include <CGAL/internal/Mean_curvature_skeleton/Utility.h>

// For correspondence tracking
#include <CGAL/internal/Mean_curvature_skeleton/Track_correspondence_visitor.h>

// For Fixed_edge_map
#include <CGAL/internal/Mean_curvature_skeleton/Fixed_edge_map.h>

// For is_collapse_ok
#include <CGAL/internal/Mean_curvature_skeleton/Collapse.h>

// For detect_degenarcy
#include <CGAL/internal/Mean_curvature_skeleton/Detect_degeneracy.h>

// Inside mesh test
#include <CGAL/Point_inside_polyhedron_3.h>

// Compute bounding box
#include <CGAL/Bbox_3.h>

#include <queue>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

// for default parameters
#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>  // for sparse linear system solver
#endif

namespace SMS = CGAL::Surface_mesh_simplification;

namespace CGAL {

/// @cond CGAL_DOCUMENT_INTERNAL

/// \ingroup PkgMeanCurvatureSkeleton3
///@brief Edge collapse algorithm tag.
enum Collapse_algorithm_tag
{
  SIMPLIFICATION, /**< algorithm from simplification package */
  LINEAR          /**< iterative linear search */
};

/// \ingroup PkgMeanCurvatureSkeleton3
///@brief Degeneracy detection algorithm tag.
enum Degeneracy_algorithm_tag
{
  HEURISTIC, /**< a simple heuristic */
  EULER      /**< counting the euler characteristic */
};

/// @endcond

/// \ingroup PkgMeanCurvatureSkeleton3
/// @brief Class providing the functionalities for extracting
///        the skeleton of a triangulated surface mesh.
///
/// \todo index pmap should also be writable :/
///
/// @tparam Traits
///         a model of `SurfaceMeshSkeletonizationTraits`
///
/// @tparam HalfedgeGraph
///         a model of `HalfedgeGraph`
/// @tparam VertexIndexMap
///         a model of `ReadablePropertyMap`
///         with `boost::graph_traits<HalfedgeGraph>::vertex_descriptor` as key and
///         `unsigned int` as value type.<br>
///         <b>%Default:</b> `boost::property_map<HalfedgeGraph, boost::vertex_index_t>::%type`.
/// @tparam HalfedgeIndexMap
///         a model of `ReadablePropertyMap`</a>
///         with `boost::graph_traits<HalfedgeGraph>::halfedge_descriptor` as key and
///         `unsigned int` as value type.<br>
///         <b>%Default:</b> `boost::property_map<HalfedgeGraph, boost::halfedge_index_t>::%type`.
/// @tparam VertexPointMap
///         a model of `ReadWritePropertyMap`
///         with `boost::graph_traits<HalfedgeGraph>::vertex_descriptor` as key and
///         `Traits::Point_3` as value type.<br>
///         <b>%Default:</b> `boost::property_map<HalfedgeGraph, boost::vertex_point_t>::%type`.
/// @tparam SparseLinearAlgebraTraits_d
///         a model of `SparseLinearAlgebraTraitsWithFactor_d`.
///         If \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available
///         and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits` is provided as default parameter.\n
/// \code
///     CGAL::Eigen_solver_traits<
///         Eigen::SparseLU<
///            CGAL::Eigen_sparse_matrix<double>::EigenType,
///            Eigen::COLAMDOrdering<int> >  >
/// \endcode
///
/// @cond CGAL_DOCUMENT_INTERNAL
/// @tparam Collapse_algorithm_tag
///         tag for selecting the edge collapse algorithm
/// @tparam Degeneracy_algorithm_tag
///         tag for selecting the degeneracy detection algorithm
/// @endcond
#ifdef DOXYGEN_RUNNING
template <class Traits,
          class HalfedgeGraph,
          class VertexIndexMap = Default,
          class HalfedgeIndexMap = Default,
          class VertexPointMap = Default,
          class SparseLinearAlgebraTraits_d = Default>
#else
template <class Traits,
          class HalfedgeGraph,
          class VertexIndexMap_ = Default,
          class HalfedgeIndexMap_  = Default,
          class VertexPointMap_ = Default,
          class SparseLinearAlgebraTraits_d_ = Default,
          Collapse_algorithm_tag Collapse_tag = LINEAR,
          Degeneracy_algorithm_tag Degeneracy_tag = EULER>
#endif
class Mean_curvature_flow_skeletonization
{
// Public types
public:
/// \name Types
/// @{
  // Template parameters
  #ifndef DOXYGEN_RUNNING
    typedef typename Default::Get<
      VertexIndexMap_,
      typename boost::property_map<HalfedgeGraph, boost::vertex_index_t>::type
    >::type VertexIndexMap;
    typedef typename Default::Get<
      HalfedgeIndexMap_,
      typename boost::property_map<HalfedgeGraph, boost::halfedge_index_t>::type
    >::type HalfedgeIndexMap;
  typedef typename Default::Get<
    VertexPointMap_,
    typename boost::property_map<HalfedgeGraph, boost::vertex_point_t>::type
  >::type VertexPointMap;
  #endif

  #ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    SparseLinearAlgebraTraits_d_,
  #if defined(CGAL_EIGEN3_ENABLED)
      CGAL::Eigen_solver_traits<
          Eigen::SparseLU<
            CGAL::Eigen_sparse_matrix<double>::EigenType,
            Eigen::COLAMDOrdering<int> >  >
  #else
    SparseLinearAlgebraTraits_d_ // no parameter provided, and Eigen is not enabled: so don't compile!
  #endif
  >::type SparseLinearAlgebraTraits_d;
  #else
  ///
  typedef SparseLinearAlgebraTraits_d_ SparseLinearAlgebraTraits_d;
  #endif
 
/// @}

  typedef typename Traits::Vector_3                                             Vector;


  // Repeat HalfedgeGraph types
  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_iterator         vertex_iterator;
  typedef typename boost::graph_traits<HalfedgeGraph>::halfedge_iterator       halfedge_iterator;
  typedef typename boost::graph_traits<HalfedgeGraph>::edge_descriptor         edge_descriptor;
  typedef typename boost::graph_traits<HalfedgeGraph>::edge_iterator           edge_iterator;
  typedef typename boost::graph_traits<HalfedgeGraph>::in_edge_iterator        in_edge_iterator;
  typedef typename boost::graph_traits<HalfedgeGraph>::out_edge_iterator       out_edge_iterator;

  typedef typename boost::graph_traits<HalfedgeGraph>::face_iterator           Facet_iterator;
  typedef Halfedge_around_face_circulator<HalfedgeGraph>            Halfedge_facet_circulator;

  // Cotangent weight calculator
  typedef typename internal::Cotangent_weight<HalfedgeGraph,
  internal::Cotangent_value_minimum_zero<HalfedgeGraph,
  internal::Cotangent_value_Meyer_secure<HalfedgeGraph> > >                    Weight_calculator;

  typedef internal::Curve_skeleton<HalfedgeGraph,
                                   VertexIndexMap,
                                   HalfedgeIndexMap,
                                   VertexPointMap>                     Skeleton;

  // Mesh simplification types
  typedef SMS::Edge_profile<HalfedgeGraph>                                     Profile;

  // Repeat Triangulation types
  typedef CGAL::Exact_predicates_exact_constructions_kernel                    K;
  typedef K::Vector_3                                                          Exact_vector;
  typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>             Vb;
  typedef CGAL::Triangulation_data_structure_3<Vb>                             Tds;
  typedef CGAL::Delaunay_triangulation_3<K, Tds>                               Delaunay;
  typedef Delaunay::Point                                                      Exact_point;
  typedef Delaunay::Cell_handle                                                Cell_handle;
  typedef Delaunay::Vertex_handle                                              TriVertex_handle;
  typedef Delaunay::Locate_type                                                Locate_type;
  typedef Delaunay::Finite_vertices_iterator                                   Finite_vertices_iterator;
  typedef Delaunay::Finite_edges_iterator                                      Finite_edges_iterator;
  typedef Delaunay::Finite_facets_iterator                                     Finite_facets_iterator;
  typedef Delaunay::Finite_cells_iterator                                      Finite_cells_iterator;

// Data members
private:

  /** a reference to the input surface mesh */
  HalfedgeGraph& m_pmesh;

  /** Storing indices of all vertices. */
  VertexIndexMap m_vertex_id_pmap;
  /** Storing indices of all edges. */
  HalfedgeIndexMap m_hedge_id_pmap;
  /** Storing the point for HalfedgeGraph vertex_descriptor. */
  VertexPointMap m_pmesh_point_pmap;

  /** Controling the velocity of movement and approximation quality. */
  double m_omega_H;
  /** Controling the smoothness of the medial approximation. */
  double m_omega_P;
  /** Edges with length less than `min_edge_length` will be collapsed. */
  double m_min_edge_length;
  /** Triangles with angle greater than `alpha_TH` will be split. */
  double m_alpha_TH;
  /** Value very close to zero. */
  double m_zero_TH;
  /** `contract_until_convergence` will stop if the change of area in one iteration
   *  is less than `delta_area`. */
  double m_delta_area;
  /** Surface area of original surface mesh. */
  double m_original_area;
  /** Maximum number of iterations. */
  int m_max_iterations;
  /** Should the skeleton be medially centered? */
  bool m_is_medially_centered;
  /** Are poles computed? */
  bool m_are_poles_computed;

  /** Cotangent weight calculator. */
  Weight_calculator m_weight_calculator;
  /** Storing the weights for edges. */
  std::vector<double> m_edge_weight;
  /** The sparse solver. */
  SparseLinearAlgebraTraits_d m_solver;

  /** Assign a unique id to a new vertex. */
  int m_vertex_id_count;
  /** The maximum id for original surface. vertices with ids
   *  greater than `m_max_id` are created during split,
   *  thus will not be considered in correspondence tracking. */
  int m_max_id;
  /** Used when assembling the matrix. */
  std::map<int, int> m_new_id;

  /** Store the id of fixed vertices. */
  std::map<size_t, bool> m_is_vertex_fixed_map;

  /** The incident angle for a halfedge. */
  std::vector<double> m_halfedge_angle;

  /** Record the correspondence between final surface
   *  and original surface points. */
  std::map<int, std::vector<int> > m_correspondence;

  /** Record the corresponding pole of a point. */
  std::map<int, int> m_poles;
  /** The normal of surface points. */
  std::vector<Vector> m_normals;
  /** The dual of a cell in Triangulation(a Voronoi point). */
  std::vector<Point> m_cell_dual;

// Private functions and classes

struct Vertex_to_point{
  VertexPointMap ppmap;
  Vertex_to_point(VertexPointMap ppmap): ppmap(ppmap){}
  typedef typename boost::property_traits<VertexPointMap>::reference result_type;
  result_type
  operator()(vertex_descriptor vd) const{
    return get(ppmap, vd);
  }
};

double diagonal_length(const Bbox_3& bbox)
{
  double dx = bbox.xmax() - bbox.xmin();
  double dy = bbox.ymax() - bbox.ymin();
  double dz = bbox.zmax() - bbox.zmin();

  double diag = dx * dx + dy * dy + dz * dz;
  return sqrt(diag);
}


double init_min_edge_length()
{
  vertex_iterator vb, ve;
  boost::tie(vb, ve) = vertices(m_pmesh);
  Vertex_to_point v_to_p(m_pmesh_point_pmap);
  Bbox_3 bbox = CGAL::bbox_3(boost::make_transform_iterator(vb, v_to_p),
                             boost::make_transform_iterator(ve, v_to_p));
  return 0.002 * diagonal_length(bbox);
}

// Public methods
public:

  /// \name Constructor
  ///@{
  #ifdef DOXYGEN_RUNNING
  /**
   * The constructor of a skeletonization object.
   *
   * The algorithm parameters are initialized such that:
   * - `omega_H()` returns 0.1
   * - `omega_P()` returns 0.2
   * - `delta_area()` returns 0.0001
   * - `max_iterations()` returns 500
   * - `is_medially_centered()` returns `true`
   * - `min_edge_length()` returns  0.002 * the length of the diagonal of the bounding box of `pmesh`
   *
   * @pre `pmesh` is a triangulated surface mesh without borders and has exactly one connected component.
   * @param pmesh 
   *        input surface mesh.
   * @param vertex_index_map 
   *        property map which associates an id to each vertex, from `0` to `num_vertices(pmesh)-1`.
   * @param halfedge_index_map
   *        property map which associates an id to each halfedge, from `0` to `num_halfedges(pmesh)-1`.
   * @param vertex_point_map 
   *        property map which associates a point to each vertex of the graph.
   * @attention  `pmesh` will be modified through edge collapse and triangle split operations while creating the skeleton.
   */
  Mean_curvature_flow_skeletonization(HalfedgeGraph& pmesh,
                                      VertexIndexMap vertex_index_map =get(boost::vertex_index, pmesh),
                                      HalfedgeIndexMap halfedge_index_map = get(boost::halfedge_index, pmesh),
                                      VertexPointMap vertex_point_map = get(boost::vertex_point, pmesh) );
  #else
  Mean_curvature_flow_skeletonization(HalfedgeGraph& pmesh,
                                      VertexIndexMap vertex_index_map,
                                      HalfedgeIndexMap halfedge_index_map,
                                      VertexPointMap vertex_point_map)
    : m_pmesh(pmesh)
    , m_vertex_id_pmap(vertex_index_map)
    , m_hedge_id_pmap(halfedge_index_map)
    , m_pmesh_point_pmap(vertex_point_map)
  {
    init_args();
    init();
  }

  Mean_curvature_flow_skeletonization(HalfedgeGraph& pmesh,
                                      VertexIndexMap vertex_index_map,
                                      HalfedgeIndexMap halfedge_index_map)
    : m_pmesh(pmesh)
    , m_vertex_id_pmap(vertex_index_map)
    , m_hedge_id_pmap(halfedge_index_map)
    , m_pmesh_point_pmap(get(boost::vertex_point, m_pmesh))
  {
    init_args();
    init();
  }

  Mean_curvature_flow_skeletonization(HalfedgeGraph& pmesh, VertexIndexMap vertex_index_map)
    : m_pmesh(pmesh)
    , m_vertex_id_pmap(vertex_index_map)
    , m_hedge_id_pmap(get(boost::halfedge_index, m_pmesh))
    , m_pmesh_point_pmap(get(boost::vertex_point, m_pmesh))
  {
    init_args();
    init();
  }

  Mean_curvature_flow_skeletonization(HalfedgeGraph& pmesh)
    : m_pmesh(pmesh)
    , m_vertex_id_pmap(get(boost::vertex_index, m_pmesh))
    , m_hedge_id_pmap(get(boost::halfedge_index, m_pmesh))
    , m_pmesh_point_pmap(get(boost::vertex_point, m_pmesh))
  {
    init_args();
    init();
  }
  #endif
  /// @} Constructor

  /// \name Algorithm Parameters
  /// @{

  /// Controls the velocity of movement and approximation quality:
  /// increasing `omega_H` makes the mean curvature flow based contraction converge
  /// faster, but results in a skeleton of lower quality.
  double omega_H()
  {
    return m_omega_H;
  }

  void set_omega_H(double value)
  {
    m_omega_H = value;
  }

  /// Controls the smoothness of the medial approximation:
  /// increasing `omega_P` results in a skeleton closer
  /// to the medial axis, but slows down the speed of contraction.
  double omega_P()
  {
    return m_omega_P;
  }

  void set_omega_P(double value) 
  {
    m_omega_P = value;
  }
  
  /// Threshold used to select edges collapsed during contraction.
  double min_edge_length()
  {
    return m_min_edge_length;
  }

  void set_min_edge_length(double value)
  {
    m_min_edge_length = value;
  }

  /// Stop criterium used by `contract_until_convergence()` to detect convergence.
  double delta_area()
  {
    return m_delta_area;
  }

  void set_delta_area(double value)
  {
    m_delta_area = value;
  }

  /// If `true`, the result skeleton is medially centered.
  bool is_medially_centered()
  {
    return m_is_medially_centered;
  }

  void set_is_medially_centered(bool value)
  {
    m_is_medially_centered = value;
  }
  
  /// Maximum number of iterations performed by `contract_until_convergence()`.
  int max_iterations()
  {
    return m_max_iterations;
  }

  void set_max_iterations(int value)
  {
    m_max_iterations = value;
  }
  
  /// Reference to the input surface mesh.
  HalfedgeGraph& halfedge_graph()
  {
    return m_pmesh;
  }

  /// \cond SKIP_FROM_MANUAL

  void set_alpha_TH(double value)
  {
    m_alpha_TH = value;
  }

  double alpha_TH()
  {
    return m_alpha_TH;
  }

  void set_zero_TH(double value)
  {
    m_zero_TH = value;
  }

  double zero_TH()
  {
    return m_zero_TH;
  }

  /// \endcond

  /// @cond CGAL_DOCUMENT_INTERNAL

  /**
   * Get the positions of fixed(degenerate) points.
   *
   * @param fixed_points
   *        return the positions of fixed points
   */
  void fixed_points(std::vector<Point>& fixed_points)
  {
    fixed_points.clear();
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = vertices(m_pmesh); vb != ve; ++vb)
    {
      int id = get(m_vertex_id_pmap, *vb);
      if (m_is_vertex_fixed_map.find(id) != m_is_vertex_fixed_map.end())
      {
        vertex_descriptor vd = *vb;
        fixed_points.push_back(get(m_pmesh_point_pmap, vd));
      }
    }
  }

  /**
   * Get the positions of non-fixed(non-degenerate) points.
   *
   * @param non_fixed_points
   *        return the positions of non-fixed points
   */
  void non_fixed_points(std::vector<Point>& non_fixed_points)
  {
    non_fixed_points.clear();
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = vertices(m_pmesh); vb != ve; ++vb)
    {
      int id = get(m_vertex_id_pmap, *vb);
      if (m_is_vertex_fixed_map.find(id) == m_is_vertex_fixed_map.end())
      {
          vertex_descriptor vd = *vb;
          non_fixed_points.push_back(get(m_pmesh_point_pmap, vd));
      }
    }
  }

  /**
   * Get the Voronoi pole for the polygonal mesh.
   *
   * @param max_poles
   *        for each surface mesh vertex, record its correspondent Voronoi pole position
   */
  void poles(std::vector<Point>& max_poles)
  {
    max_poles.resize(num_vertices(m_pmesh));
    vertex_iterator vb, ve;
    int cnt = 0;
    for (boost::tie(vb, ve) = vertices(m_pmesh); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int vid = get(m_vertex_id_pmap, v);
      max_poles[cnt++] = m_cell_dual[m_poles[vid]];
    }
  }

  /// @endcond

  /// @} Setter and Getter


  /// \name High Level Function
  /// @{


  /**
   * Creates the curve skeleton: the input surface mesh is iteratively
   * contracted until convergence, and then turned into a curve skeleton.
   *
   * This is equivalent to calling `contract_until_convergence()` and `convert_to_skeleton()`.
   *
   * @param skeleton
   *        graph that will contain the skeleton of the input surface mesh
   * @param skeleton_points
   *        property map containing the location of the vertices of the graph `skeleton`
   * @param skeleton_to_pmesh_vertices property map associating a vertex `v` of the graph `skeleton`
   *         to the indices of the vertices of the input surface mesh corresponding to `v`.
   *         Indices are retrieved using the vertex index property map given in the constructor.
   * @tparam Graph
   *         an instantiation of <A href="http://www.boost.org/libs/graph/doc/adjacency_list.html">`boost::adjacency_list`</a>
   *         as data structure for the skeleton curve.
   * @tparam GraphVertexPointMap
   *         a model of `ReadWritePropertyMap`
   *         with `boost::graph_traits<Graph>::vertex_descriptor` as key type and
   *         `Traits::Point_3` as value type  
   * @tparam GraphVertexIndicesMap
   *         a model of `ReadWritePropertyMap`
   *         with `boost::graph_traits<Graph>::vertex_descriptor` as key type and
   *         `std::vector<size_type>` as value type,
   *         where `size_type` is the value type of `VertexIndexMap`
   */
  template <class Graph, class GraphVertexPointMap, class GraphVertexIndicesMap>
  void operator()(Graph& skeleton,
                  GraphVertexPointMap& skeleton_points,
                  GraphVertexIndicesMap& skeleton_to_pmesh_vertices)
  {
    contract_until_convergence();
    convert_to_skeleton(skeleton, skeleton_points, skeleton_to_pmesh_vertices);
  }
  /// @}
  
  /// \name Low Level Functions
  /// The following functions enable the user to run the mean curvature flow skeleton algorithm step by step.
  /// @{

  /**
   * Contract the surface mesh by mean curvature flow.
   */
  void contract_geometry()
  {
    MCFSKEL_DEBUG(std::cerr << "before contract geometry";)

    update_vertex_id();

    compute_edge_weight();

    int nver = num_vertices(m_pmesh);
    int nrows;
    if (m_is_medially_centered)
    {
      nrows = nver * 3;
      if (!m_are_poles_computed)
      {
        compute_voronoi_pole();
      }
    }
    else
    {
      nrows = nver * 2;
    }
    // Assemble linear system At * A * X = At * B
    typename SparseLinearAlgebraTraits_d::Matrix A(nrows, nver);
    assemble_LHS(A);

    typename SparseLinearAlgebraTraits_d::Vector X(nver), Bx(nrows);
    typename SparseLinearAlgebraTraits_d::Vector Y(nver), By(nrows);
    typename SparseLinearAlgebraTraits_d::Vector Z(nver), Bz(nrows);
    assemble_RHS(Bx, By, Bz);

    MCFSKEL_DEBUG(std::cerr << "before solve\n";)

    // solve "At * A * X = At * B".
    double D;
    m_solver.pre_factor_non_symmetric(A, D);
    m_solver.linear_solver_non_symmetric(A, Bx, X);
    m_solver.linear_solver_non_symmetric(A, By, Y);
    m_solver.linear_solver_non_symmetric(A, Bz, Z);

    MCFSKEL_DEBUG(std::cerr << "after solve\n";)

    // copy to surface mesh
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = vertices(m_pmesh); vb != ve; ++vb)
    {
      vertex_descriptor vi = *vb;
      int id = get(m_vertex_id_pmap, vi);
      int i = m_new_id[id];
      Point p(X[i], Y[i], Z[i]);
      put(m_pmesh_point_pmap, vi, p);
    }

    MCFSKEL_DEBUG(std::cerr << "leave contract geometry\n";)
  }

  /**
   * Collapses edges with length less than `min_edge_length()` and returns the number of edges collapsed.
   */
  int collapse_edges()
  {
    internal::Fixed_edge_map<HalfedgeGraph> fixed_edge_map(m_pmesh);
    init_fixed_edge_map(fixed_edge_map);

    int num_collapses = 0;
    while (true)
    {
      int cnt;
      if (Collapse_tag == SIMPLIFICATION)
      {
        cnt = collapse_edges_simplification();
      }
      else if (Collapse_tag == LINEAR)
      {
        cnt = collapse_edges_linear(fixed_edge_map);
      }

      if (cnt == 0)
      {
        break;
      }
      else
      {
        num_collapses += cnt;
      }
    }
    return num_collapses;
  }

  /**
   * Splits faces having one angle greater than `alpha_TH()` and returns the number of faces split.
   * \todo define alpha_TH.
   */
  int split_faces()
  {
    MCFSKEL_DEBUG(std::cerr << "before split\n";)

    int num_splits = 0;
    while (true)
    {
      if (num_vertices(m_pmesh) <= 3)
      {
        break;
      }
      int cnt = split_flat_triangle();
      if (cnt == 0)
      {
        break;
      }
      else
      {
        num_splits += cnt;
      }
    }

    MCFSKEL_DEBUG(std::cerr << "after split\n";)

    return num_splits;
  }

#ifndef DOXYGEN_RUNNING
  /**
   * Sequentially calls `collapse_edges()` and `split_faces()` and returns the number of edges collapsed and faces split.
   */
  int remesh()
  {
    MCFSKEL_DEBUG(std::cerr << "before collapse edges\n";)

    int num_collapses = collapse_edges();
    MCFSKEL_INFO(std::cerr << "collapse " << num_collapses << " edges.\n";)

    int num_splits = split_faces();
    MCFSKEL_INFO(std::cerr << "split " << num_splits << " edges.\n";)

    return num_collapses + num_splits;
  }
#endif

  /**
   * Fixes the position of degenerate vertices and returns the number of newly fixed vertices.
   * \todo int -> size_t
   */
  int detect_degeneracies()
  {
    if (Degeneracy_tag == HEURISTIC)
    {
      return detect_degeneracies_heuristic();
    }
    else if (Degeneracy_tag == EULER)
    {
      return detect_degeneracies_in_disk();
    }
  }

#ifndef DOXYGEN_RUNNING
  void contract()
  {
    contract_geometry();
    remesh();
    detect_degeneracies();

    MCFSKEL_DEBUG(print_edges();)

    MCFSKEL_INFO(double area = internal::get_surface_area(m_pmesh, m_pmesh_point_pmap);)
    MCFSKEL_INFO(std::cout << "area " << area << "\n";)
  }
#endif

  /**
   * Iteratively calls the sequence `contract_geometry()`,  `collapse_edges()`, `split_faces()`, and `detect_degeneracies()`
   * until the change of surface area during one
   * iteration is less than `delta_area()` * original surface area.
   */
  void contract_until_convergence()
  {
    double last_area = 0;
    int num_iteration = 0;
    while (true)
    {
      MCFSKEL_INFO(std::cout << "iteration " << num_iteration + 1 << "\n";)

      contract_geometry();
      remesh();
      detect_degeneracies();

      double area = internal::get_surface_area(m_pmesh, m_pmesh_point_pmap);
      double area_ratio = fabs(last_area - area) / m_original_area;

      MCFSKEL_INFO(std::cout << "area " << area << "\n";)
      MCFSKEL_INFO(std::cout << "|area - last_area| / original_area "
                             << area_ratio << "\n";)

      if (area_ratio < m_delta_area)
      {
        break;
      }
      last_area = area;

      num_iteration++;
      if (num_iteration >= m_max_iterations)
      {
        break;
      }
    }
  }

  /**
   * Converts the contracted surface mesh to a skeleton curve.
   * @param skeleton
   *       graph that will contain the skeleton of `pmesh`
   * @param skeleton_points
   *        property map containing the location of the vertices of the graph `skeleton`
   *@param skeleton_to_pmesh_vertices property map associating a vertex `v` of the graph `skeleton`
   *       to the set of vertices of `pmesh` corresponding to `v`.
   * @tparam Graph
   *         an instantiation of <A href="http://www.boost.org/libs/graph/doc/adjacency_list.html">`boost::adjacency_list`</a>
   *         as a data structure for the skeleton curve.
   * @tparam GraphVertexIndicesMap
   *         a model of `ReadWritePropertyMap`
   *         with `boost::graph_traits<Graph>::vertex_descriptor as` key type and
   *         `std::vector<size_type>` as value type,
   *         where `size_type` is the value type of `VertexIndexMap`
   * @tparam GraphPointMap
   *         a model of `ReadWritePropertyMap`
   *         with `boost::graph_traits<Graph>::vertex_descriptor` as key type and
   *         `Traits::Point_3` as value type  
   */
  template <class Graph, class GraphVertexPointMap, class GraphVertexIndicesMap>
  void convert_to_skeleton(Graph& skeleton, GraphVertexPointMap& skeleton_points, GraphVertexIndicesMap& skeleton_to_pmesh_vertices)
  {
    Skeleton skeletonization(m_pmesh, m_vertex_id_pmap, m_hedge_id_pmap, m_pmesh_point_pmap);
    std::map<typename boost::graph_traits<Graph>::vertex_descriptor, std::vector<int> > skeleton_to_surface_map;
    
    skeletonization.extract_skeleton(skeleton, skeleton_points, skeleton_to_surface_map);
    correspondent_vertices(skeleton_to_surface_map, skeleton_to_pmesh_vertices);
  }

  /// @} Public Algorithm API


private:

  // --------------------------------------------------------------------------
  // Initialization
  // --------------------------------------------------------------------------

  /// Initialize the parameters for Mean_curvature_flow_skeletonization
  void init_args()
  {
    m_omega_H = 0.1;
    m_omega_P = 0.2;
    m_delta_area = 0.0001;
    m_max_iterations = 500;
    m_is_medially_centered = true;
    m_min_edge_length =  init_min_edge_length();
    m_alpha_TH = 110;
    m_zero_TH = 1e-7;
  }

  /// Initialize some global data structures such as vertex id.
  void init()
  {
    m_are_poles_computed = false;

    m_alpha_TH *= (CGAL_PI / 180.0);
    m_original_area = internal::get_surface_area(m_pmesh, m_pmesh_point_pmap);

    m_vertex_id_count = num_vertices(m_pmesh);
    m_max_id = m_vertex_id_count;

    m_is_vertex_fixed_map.clear();
    m_correspondence.clear();

    if (m_is_medially_centered)
      compute_voronoi_pole();
  }

  // --------------------------------------------------------------------------
  // Contraction
  // --------------------------------------------------------------------------

  /// Compute cotangent weights of all edges.
  void compute_edge_weight()
  {
    m_edge_weight.clear();
    m_edge_weight.reserve(2 * num_edges(m_pmesh));
    halfedge_iterator eb, ee;
    for(boost::tie(eb, ee) = halfedges(m_pmesh); eb != ee; ++eb)
    {
      m_edge_weight.push_back(this->m_weight_calculator(*eb, m_pmesh));
    }
  }

  /// Assemble the left hand side.
  void assemble_LHS(typename SparseLinearAlgebraTraits_d::Matrix& A)
  {
    MCFSKEL_DEBUG(std::cerr << "start LHS\n";)

    int nver = num_vertices(m_pmesh);

    Point_inside_polyhedron_3<HalfedgeGraph, Traits> test_inside(m_pmesh);

    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = vertices(m_pmesh); vb != ve; ++vb)
    {
      int id = get(m_vertex_id_pmap, *vb);

      int i = m_new_id[id];
      // if the vertex is fixed
      if (m_is_vertex_fixed_map.find(id) != m_is_vertex_fixed_map.end())
      {
        A.set_coef(i + nver, i, 1.0 / m_zero_TH, true);
      }
      else
      {
        A.set_coef(i + nver, i, m_omega_H, true);
        if (m_is_medially_centered)
        {
          if (id < m_max_id)
          {
            if (test_inside(m_cell_dual[m_poles[id]]) == CGAL::ON_BOUNDED_SIDE)
            {
              A.set_coef(i + nver * 2, i, m_omega_P, true);
            }
          }
        }
      }
    }

    for (boost::tie(vb, ve) = vertices(m_pmesh); vb != ve; ++vb)
    {
      int id = get(m_vertex_id_pmap, *vb);
      int i = m_new_id[id];
      double L = 1.0;
      // if the vertex is fixed
      if (m_is_vertex_fixed_map.find(id) != m_is_vertex_fixed_map.end())
      {
        L = 0;
      }
      double diagonal = 0;
      in_edge_iterator e, e_end;
      for (boost::tie(e, e_end) = in_edges(*vb, m_pmesh); e != e_end; ++e)
      {
        vertex_descriptor vj = source(*e, m_pmesh);
        double wij = m_edge_weight[get(m_hedge_id_pmap, halfedge(*e, m_pmesh))] * 2.0;
        int jd = get(m_vertex_id_pmap, vj);
        int j = m_new_id[jd];
        A.set_coef(i, j, wij * L, true);
        diagonal += -wij;
      }
      A.set_coef(i, i, diagonal, true);
    }

    MCFSKEL_DEBUG(std::cerr << "end LHS\n";)
  }

  /// Assemble the right hand side.
  void assemble_RHS(typename SparseLinearAlgebraTraits_d::Vector& Bx,
                    typename SparseLinearAlgebraTraits_d::Vector& By,
                    typename SparseLinearAlgebraTraits_d::Vector& Bz)
  {
    MCFSKEL_DEBUG(std::cerr << "start RHS\n";)

    Point_inside_polyhedron_3<HalfedgeGraph, Traits> test_inside(m_pmesh);

    // assemble right columns of linear system
    int nver = num_vertices(m_pmesh);
    vertex_iterator vb, ve;
    for (int i = 0; i < nver; ++i)
    {
      Bx[i] = 0;
      By[i] = 0;
      Bz[i] = 0;
    }

    for (boost::tie(vb, ve) = vertices(m_pmesh); vb != ve; ++vb)
    {
      vertex_descriptor vi = *vb;
      int id = get(m_vertex_id_pmap, vi);
      int i = m_new_id[id];

      double oh, op = 0.0;
      if (m_is_vertex_fixed_map.find(id) != m_is_vertex_fixed_map.end())
      {
        oh = 1.0 / m_zero_TH;
      }
      else
      {
        oh = m_omega_H;
        if (m_is_medially_centered)
        {
          if (id < m_max_id)
          {
            if (test_inside(m_cell_dual[m_poles[id]]) == CGAL::ON_BOUNDED_SIDE)
            {
              op = m_omega_P;
            }
          }
        }
      }
      Bx[i + nver] = get(m_pmesh_point_pmap, vi).x() * oh;
      By[i + nver] = get(m_pmesh_point_pmap, vi).y() * oh;
      Bz[i + nver] = get(m_pmesh_point_pmap, vi).z() * oh;
      if (m_is_medially_centered)
      {
        double x = to_double(m_cell_dual[m_poles[id]].x());
        double y = to_double(m_cell_dual[m_poles[id]].y());
        double z = to_double(m_cell_dual[m_poles[id]].z());
        Bx[i + nver * 2] = x * op;
        By[i + nver * 2] = y * op;
        Bz[i + nver * 2] = z * op;
      }
    }

    MCFSKEL_DEBUG(std::cerr << "end RHS\n";)
  }

  /// The order of vertex id is the same as the traverse order.
  void update_vertex_id()
  {
    m_new_id.clear();
    int cnt = 0;
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = vertices(m_pmesh); vb != ve; ++vb)
    {
      int id = get(m_vertex_id_pmap, *vb);
      m_new_id[id] = cnt++;
    }
  }

  // --------------------------------------------------------------------------
  // Edge collapse
  // --------------------------------------------------------------------------

  /// Collapse short edges using simplification package.
  int collapse_edges_simplification()
  {
    internal::Fixed_edge_map<HalfedgeGraph> fixed_edge_map(m_pmesh);

    init_fixed_edge_map(fixed_edge_map);


    int edge_id = -1;
    halfedge_iterator hb, he;
    for (boost::tie(hb, he) = halfedges(m_pmesh); hb != he; ++hb)
    {
      put(m_hedge_id_pmap, *hb, ++edge_id);
    }

    // This is a stop predicate (defines when the algorithm terminates).
    // The simplification stops when the length of all edges is greater
    // than the minimum threshold.
    CGAL::internal::Minimum_length_predicate<HalfedgeGraph> stop(m_min_edge_length);

    // midpoint placement without geometric test
    SMS::Geometric_test_skipper< SMS::Midpoint_placement<HalfedgeGraph> > placement;

    internal::Track_correspondence_visitor<HalfedgeGraph, VertexPointMap> vis;
    if (m_is_medially_centered)
    {
      vis = internal::Track_correspondence_visitor<HalfedgeGraph, VertexPointMap>
            (&m_pmesh_point_pmap, &m_correspondence, &m_poles, &m_cell_dual, m_max_id);
    }
    else
    {
      vis = internal::Track_correspondence_visitor<HalfedgeGraph, VertexPointMap>
            (&m_pmesh_point_pmap, &m_correspondence, m_max_id);
    }

    int r = SMS::edge_collapse
                (m_pmesh
                ,stop
                ,CGAL::get_cost(SMS::Edge_length_cost<HalfedgeGraph>())
                      .get_placement(placement)
                      .visitor(vis)
                      .edge_is_constrained_map(fixed_edge_map)
                );

    return r;
  }

  /// Track correspondent original surface points during collapse.
  void track_correspondence(vertex_descriptor v0, vertex_descriptor v1,
                            vertex_descriptor v)
  {
    int id0 = get(m_vertex_id_pmap, v0);
    int id1 = get(m_vertex_id_pmap, v1);
    int vid = get(m_vertex_id_pmap, v);
    int from, to;
    if (id0 == vid)
    {
      from = id1;
      to = id0;
    }
    else if (id1 == vid)
    {
      from = id0;
      to = id1;
    }

    if (m_correspondence.find(to) == m_correspondence.end())
    {
      m_correspondence[to] = std::vector<int>();
    }
    // only track vertex in original surface mesh
    if (from < m_max_id)
    {
      m_correspondence[to].push_back(from);
    }
    std::map<int, std::vector<int> >::iterator iter = m_correspondence.find(from);
    if (iter != m_correspondence.end())
    {
      for (size_t i = 0; i < (iter->second).size(); ++i)
      {
        m_correspondence[to].push_back((iter->second)[i]);
      }
      (iter->second).clear();
      m_correspondence.erase(iter);
    }

    if (m_is_medially_centered)
    {
      Point pole0 = Point(to_double(m_cell_dual[m_poles[id0]].x()),
                          to_double(m_cell_dual[m_poles[id0]].y()),
                          to_double(m_cell_dual[m_poles[id0]].z()));
      Point pole1 = Point(to_double(m_cell_dual[m_poles[id1]].x()),
                          to_double(m_cell_dual[m_poles[id1]].y()),
                          to_double(m_cell_dual[m_poles[id1]].z()));
      Point p1 = get(m_pmesh_point_pmap, v1);
      double dis_to_pole0 = sqrt(squared_distance(pole0, p1));
      double dis_to_pole1 = sqrt(squared_distance(pole1, p1));
      if (dis_to_pole0 < dis_to_pole1)
      {
        m_poles[id1] = m_poles[id0];
      }
      std::map<int, int>::iterator pole_iter = m_poles.find(id0);
      m_poles.erase(pole_iter);
    }
  }

  /// Collapse short edges by iteratively linear search.
  int collapse_edges_linear(internal::Fixed_edge_map<HalfedgeGraph>& fixed_edge_map)
  {
    std::vector<edge_descriptor> all_edges;
    all_edges.reserve(num_edges(m_pmesh));
    edge_iterator eb, ee;

    boost::tie(eb, ee) = edges(m_pmesh);
    std::copy(eb, ee, std::back_inserter(all_edges));

    int cnt = 0;
    for (size_t i = 0; i < all_edges.size(); ++i)
    {
      halfedge_descriptor h = halfedge(all_edges[i], m_pmesh);
      if (fixed_edge_map.is_fixed(h))
      {
        continue;
      }

      vertex_descriptor vi = source(h, m_pmesh);
      vertex_descriptor vj = target(h, m_pmesh);
      double edge_length = sqrt(squared_distance(get(m_pmesh_point_pmap, vi),
                                                 get(m_pmesh_point_pmap, vj)));
      if (internal::is_collapse_ok(m_pmesh, h) && edge_length < m_min_edge_length)
      {
        Point p = midpoint(
          get(vertex_point, m_pmesh, source(h, m_pmesh)),
          get(vertex_point, m_pmesh, target(h, m_pmesh)));

        // invalidate the edges that will be collapsed
        // since the surface mesh is closed, 6 halfedges will be collapsed
        // (opposite is automatically added)
        fixed_edge_map.set_is_fixed(h, true);
        fixed_edge_map.set_is_fixed(prev(h, m_pmesh), true);
        fixed_edge_map.set_is_fixed(prev(opposite(h, m_pmesh), m_pmesh), true);

        vertex_descriptor v = Euler::collapse_edge(edge(h, m_pmesh), m_pmesh);
        put(m_pmesh_point_pmap, v, p);

        track_correspondence(vi, vj, v);

        cnt++;
      }
    }

    return cnt;
  }

  /// Fix an edge if both incident vertices are degenerate.
  void init_fixed_edge_map(internal::Fixed_edge_map<HalfedgeGraph>& fixed_edge_map)
  {
    edge_iterator eb, ee;
    for (boost::tie(eb, ee) = edges(m_pmesh); eb != ee; ++eb)
    {
      halfedge_descriptor h = halfedge(*eb, m_pmesh);
      vertex_descriptor vi = source(h, m_pmesh);
      vertex_descriptor vj = target(h, m_pmesh);
      size_t vi_idx = get(m_vertex_id_pmap, vi);
      size_t vj_idx = get(m_vertex_id_pmap, vj);

      if (m_is_vertex_fixed_map.find(vi_idx) != m_is_vertex_fixed_map.end()
       && m_is_vertex_fixed_map.find(vj_idx) != m_is_vertex_fixed_map.end())
      {
        fixed_edge_map.set_is_fixed(h, true); // opposite is automatically added
      }
    }
  }

  // --------------------------------------------------------------------------
  // Triangle split
  // --------------------------------------------------------------------------

  /// Compute the incident angles for all the halfedges.
  void compute_incident_angle()
  {
    m_halfedge_angle.clear();
    int ne = 2 * num_edges(m_pmesh);
    m_halfedge_angle.resize(ne, 0);

    halfedge_iterator eb, ee;
    int idx = 0;
    for (boost::tie(eb, ee) = halfedges(m_pmesh); eb != ee; ++eb)
    {
      put(m_hedge_id_pmap, *eb, idx++);
    }

    for (boost::tie(eb, ee) = halfedges(m_pmesh); eb != ee; ++eb)
    {
      int e_id = get(m_hedge_id_pmap, *eb);
      halfedge_descriptor ed = *eb;

      if (is_border(ed, m_pmesh))
      {
        m_halfedge_angle[e_id] = -1;
      }
      else
      {
        vertex_descriptor vi = source(ed, m_pmesh);
        vertex_descriptor vj = target(ed, m_pmesh);
        halfedge_descriptor ed_next = next(ed, m_pmesh);
        vertex_descriptor vk = target(ed_next, m_pmesh);
        Point pi = get(m_pmesh_point_pmap, vi);
        Point pj = get(m_pmesh_point_pmap, vj);
        Point pk = get(m_pmesh_point_pmap, vk);

        double dis2_ij = squared_distance(pi, pj);
        double dis2_ik = squared_distance(pi, pk);
        double dis2_jk = squared_distance(pj, pk);
        double dis_ij = sqrt(dis2_ij);
        double dis_ik = sqrt(dis2_ik);
        double dis_jk = sqrt(dis2_jk);

        // A degenerate triangle will never undergo a split (but rather a collapse...)
        if (dis_ij < m_zero_TH || dis_ik < m_zero_TH || dis_jk < m_zero_TH)
        {
          m_halfedge_angle[e_id] = -1;
        }
        else
        {
          m_halfedge_angle[e_id] =
              acos((dis2_ik + dis2_jk - dis2_ij) / (2.0 * dis_ik * dis_jk));
        }
      }
    }
  }

  /// Project the vertex `vk` to the line of `vs` and `vt`.
  Point project_vertex(const vertex_descriptor vs,
                       const vertex_descriptor vt,
                       const vertex_descriptor vk)
  {
    Point ps = get(m_pmesh_point_pmap, vs);
    Point pt = get(m_pmesh_point_pmap, vt);
    Point pk = get(m_pmesh_point_pmap, vk);
    CGAL::internal::Vector vec_st = CGAL::internal::Vector(ps, pt);
    CGAL::internal::Vector vec_sk = CGAL::internal::Vector(ps, pk);

    vec_st.normalize();
    double t = vec_st.dot(vec_sk);
    Point st = Point(vec_st[0] * t, vec_st[1] * t, vec_st[2] * t);
    Point pn = Point(ps[0] + st[0], ps[1] + st[1], ps[2] + st[2]);

    // project the pole
    if (m_is_medially_centered)
    {
      int sid = get(m_vertex_id_pmap, vs);
      int tid = get(m_vertex_id_pmap, vt);
      Point pole_s = m_cell_dual[m_poles[sid]];
      Point pole_t = m_cell_dual[m_poles[tid]];
      Vector pole_st = pole_t - pole_s;
      Vector p_projector = pole_st / sqrt(pole_st.squared_length());
      Point pole_n = pole_s + p_projector * t;
      m_poles[m_vertex_id_count] = m_cell_dual.size();
      m_cell_dual.push_back(pole_n);
    }
    return pn;
  }

  /// Split triangles with an angle greater than `alpha_TH`.
  int split_flat_triangle()
  {
    int ne = 2 * num_edges(m_pmesh);
    compute_incident_angle();

    int cnt = 0;
    halfedge_iterator eb, ee;
    /// \todo this is unsafe, we loop over a sequence that we modify!!!
    for (boost::tie(eb, ee) = halfedges(m_pmesh); eb != ee; ++eb)
    {
      halfedge_descriptor ei = *eb;
      halfedge_descriptor ej = opposite(ei, m_pmesh);
      int ei_id = get(m_hedge_id_pmap, ei);
      int ej_id = get(m_hedge_id_pmap, ej);
      if (ei_id < 0 || ei_id >= ne
       || ej_id < 0 || ej_id >= ne)
      {
        continue;
      }

      vertex_descriptor vs = source(ei, m_pmesh);
      vertex_descriptor vt = target(ei, m_pmesh);

      double angle_i = m_halfedge_angle[ei_id];
      double angle_j = m_halfedge_angle[ej_id];
      if (angle_i < m_alpha_TH || angle_j < m_alpha_TH)
      {
        continue;
      }

      halfedge_descriptor ek;
      if (angle_i > angle_j)
      {
        ek = next(ei, m_pmesh);
      }
      else
      {
        ek = next(ej, m_pmesh);
      }
      vertex_descriptor vk = target(ek, m_pmesh);
      Point pn = project_vertex(vs, vt, vk);
      halfedge_descriptor en = internal::mesh_split(m_pmesh, m_pmesh_point_pmap, ei, pn);
      // set id for new vertex
      put(m_vertex_id_pmap, target(en,m_pmesh), m_vertex_id_count++);
      cnt++;
    }
    return cnt;
  }

  // --------------------------------------------------------------------------
  // Degeneracy detection
  // --------------------------------------------------------------------------

  /// Test degeneracy of a vertex by counting the euler characteristic of
  /// its local neighborhood disk.
  int detect_degeneracies_in_disk()
  {
    int num_fixed = 0;
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = vertices(m_pmesh); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int idx = get(m_vertex_id_pmap, v);

      if (m_is_vertex_fixed_map.find(idx) == m_is_vertex_fixed_map.end())
      {
        bool willbefixed = internal::is_vertex_degenerate(m_pmesh, m_pmesh_point_pmap,
                                                          v, m_min_edge_length);
        if (willbefixed)
        {
          m_is_vertex_fixed_map[idx] = willbefixed;
          num_fixed++;
        }
      }
    }

    MCFSKEL_INFO(std::cerr << "fixed " << num_fixed << " vertices.\n";)

    return num_fixed;
  }

  /// Test degeneracy of a vertex by a simple heuristic looking for a
  /// triangular cross section.
  int detect_degeneracies_heuristic()
  {
    int num_fixed = 0;
    double elength_fixed = m_min_edge_length;
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = vertices(m_pmesh); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int idx = boost::get(m_vertex_id_pmap, v);
      if (m_is_vertex_fixed_map.find(idx) == m_is_vertex_fixed_map.end())
      {
        bool willbefixed = false;
        int bad_counter = 0;

        in_edge_iterator eb, ee;
        for (boost::tie(eb, ee) = in_edges(v, m_pmesh); eb != ee; ++eb)
        {
          halfedge_descriptor edge = halfedge(*eb, m_pmesh);
          vertex_descriptor v0 = source(edge, m_pmesh);
          vertex_descriptor v1 = target(edge, m_pmesh);
          double length = sqrt(squared_distance(get(m_pmesh_point_pmap, v0),
                                                get(m_pmesh_point_pmap, v1)));
          if (length < elength_fixed)
          {
            if (!internal::is_collapse_ok(m_pmesh, edge))
            {
              bad_counter++;
            }
          }
        }
        willbefixed = (bad_counter >= 2);
        if (willbefixed)
        {
          m_is_vertex_fixed_map[idx] = willbefixed;
          num_fixed++;
        }
      }
    }

    MCFSKEL_INFO(std::cerr << "fixed " << num_fixed << " vertices.\n";)

    return num_fixed;
  }

  // --------------------------------------------------------------------------
  // Voronoi pole
  // --------------------------------------------------------------------------

  /// Compute the Voronoi pole for surface vertices. The pole is the furthest
  /// vertex in the Voronoi cell containing the given vertex.
  void compute_voronoi_pole()
  {
    MCFSKEL_DEBUG(std::cout << "start compute_voronoi_pole\n";)
    compute_vertex_normal();

    std::vector<std::pair<Exact_point, unsigned> > points;
    std::vector<std::vector<int> > point_to_pole;

    points.clear();
    m_cell_dual.clear();
    point_to_pole.clear();
    point_to_pole.resize(num_vertices(m_pmesh));

    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = vertices(m_pmesh); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int vid = get(m_vertex_id_pmap, v);
      Exact_point tp((get(m_pmesh_point_pmap, v)).x(),
                     (get(m_pmesh_point_pmap, v)).y(),
                     (get(m_pmesh_point_pmap, v)).z());
      points.push_back(std::make_pair(tp, vid));
    }

    Delaunay T(points.begin(), points.end());

    Finite_cells_iterator cit;
    int cell_id = 0;
    for (cit = T.finite_cells_begin(); cit != T.finite_cells_end(); ++cit)
    {
      Cell_handle cell = cit;
      Exact_point point = T.dual(cell);
      Point pt(to_double(point.x()), to_double(point.y()), to_double(point.z()));
      m_cell_dual.push_back(pt);
      // each cell has 4 incident vertices
      for (int i = 0; i < 4; ++i)
      {
        TriVertex_handle vt = cell->vertex(i);
        int id = vt->info();
        point_to_pole[id].push_back(cell_id);
      }
      cell_id++;
    }

    m_poles.clear();
    for (size_t i = 0; i < point_to_pole.size(); ++i)
    {
      Point surface_point = Point(to_double(points[i].first.x()),
                                  to_double(points[i].first.y()),
                                  to_double(points[i].first.z()));

      double max_neg_t = 1;
      int max_neg_i = -1;

      for (size_t j = 0; j < point_to_pole[i].size(); ++j)
      {
        int pole_id = point_to_pole[i][j];
        Point cell_point = m_cell_dual[pole_id];
        Vector vt = cell_point - surface_point;
        Vector n = m_normals[i];

        double t = vt * n;

        // choose the one with maximum distance along the normal
        if (t < 0 && t < max_neg_t)
        {
          max_neg_i = pole_id;
          max_neg_t = t;
        }
      }

      m_poles[i] = max_neg_i;
    }
    m_are_poles_computed = true;
  }

  /// Compute an approximate vertex normal for all vertices.
  void compute_vertex_normal()
  {
    m_normals.resize(num_vertices(m_pmesh));

    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = vertices(m_pmesh); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int vid = get(m_vertex_id_pmap, v);
      m_normals[vid] = internal::get_vertex_normal<typename HalfedgeGraph::Vertex,Traits>(*v);
    }
  }

  // --------------------------------------------------------------------------
  // Vertex info association
  // --------------------------------------------------------------------------

  template <class Skeleton_vertex_descriptor, class GraphVerticesMap>
  void correspondent_vertices(std::map<Skeleton_vertex_descriptor, std::vector<int> >& skeleton_to_surface_map,
                              GraphVerticesMap& skeleton_to_surface)
  {
    typename std::map<Skeleton_vertex_descriptor, std::vector<int> >::iterator iter;
    for (iter = skeleton_to_surface_map.begin();
         iter != skeleton_to_surface_map.end(); ++iter)
    {
      Skeleton_vertex_descriptor i = iter->first;

      skeleton_to_surface[i] = std::vector<int>();
      for (size_t j = 0; j < skeleton_to_surface_map[i].size(); ++j)
      {
        int id = skeleton_to_surface_map[i][j];
        if (m_correspondence.find(id) != m_correspondence.end())
        {
          skeleton_to_surface[i].insert(skeleton_to_surface[i].end(),
                                        m_correspondence[id].begin(),
                                        m_correspondence[id].end());
        }

        if (id < m_max_id)
        {
          skeleton_to_surface[i].push_back(id);
        }
      }
    }
  }


  // --------------------------------------------------------------------------
  // Debug
  // --------------------------------------------------------------------------

  void print_edges()
  {
    halfedge_iterator eb, ee;

    std::map<halfedge_descriptor, bool> visited;

    for (boost::tie(eb, ee) = halfedges(m_pmesh); eb != ee; ++eb)
    {
      if (!visited[*eb])
      {
        vertex_descriptor vi = source(*eb, m_pmesh);
        vertex_descriptor vj = target(*eb, m_pmesh);
        size_t vi_idx = get(m_vertex_id_pmap, vi);
        size_t vj_idx = get(m_vertex_id_pmap, vj);
        std::cout << vi_idx << " " << vj_idx << "\n";

        visited[*eb] = true;
        visited[opposite(*eb,m_pmesh)] = true;
      }
    }
  }
};

} //namespace CGAL

#endif // CGAL_MEAN_CURVATURE_SKELETON_H
