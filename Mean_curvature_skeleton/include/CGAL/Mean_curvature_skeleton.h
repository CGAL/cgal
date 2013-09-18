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

/**
 * @file Mean_curvature_skeleton.h
 * @brief The class `MCF_Skeleton` containing the API to extract
 * curve skeleton for a closed triangular mesh.
 */

#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <CGAL/Default.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>

#include <boost/iterator/transform_iterator.hpp>

// Compute cotangent Laplacian
#include <CGAL/internal/Mean_curvature_skeleton/Weights.h>

// Compute the vertex normal
#include <CGAL/internal/Mean_curvature_skeleton/get_normal.h>

// Low level collapse function
#include <CGAL/Surface_mesh_simplification/halfedge_collapse_Polyhedron_3.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/internal/Mean_curvature_skeleton/Edge_minimum_length_stop_predicate.h>

// Non-default cost and placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>

// Skip the geometric test
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Geometric_test_skipper.h>

// Adaptor for Polyhedron_3
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>

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
///@brief Define the default sparse linear systems solver type.
template <class FT>
struct MCF_default_solver
{
  /// Default solver type.
  typedef CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT<typename CGAL::Eigen_sparse_matrix<FT>::EigenType> > type;
};

/// \ingroup PkgMeanCurvatureSkeleton3
///@brief Define the default HalfedgeGraphPointPMap type.
template <class Polyhedron>
struct MCF_default_halfedge_graph_pmap
{
  /// Default HalfedgeGraphPointPMap type.
  typedef typename boost::property_map<Polyhedron, CGAL::vertex_point_t>::type type;
};

/// \ingroup PkgMeanCurvatureSkeleton3
/// @brief Class providing necessary parameters for `MCF_Skeleton`.
///
/// @tparam HalfedgeGraph
///         a model of `HalfedgeGraph`
template<class HalfedgeGraph>
class MCF_skel_args
{
public:
  typedef CGAL::Bbox_3                                                   Bbox;
  typedef typename HalfedgeGraph::Traits                                 Kernel;
  typedef typename Kernel::Point_3                                       Point;
  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_iterator   vertex_iterator;

private:
  /// Functor used to transform the vertex iterator to point iterator.
  static Point v_to_p(vertex_descriptor v)
  {
    return v->point();
  }

  /// Compute the diagonal length for a bounding box.
  double diagonal(Bbox& bbox)
  {
    double dx = bbox.xmax() - bbox.xmin();
    double dy = bbox.ymax() - bbox.ymin();
    double dz = bbox.zmax() - bbox.zmin();

    double diag = dx * dx + dy * dy + dz * dz;
    return sqrt(diag);
  }

public:

  /**
   * The constructor of a MCF_skel_args object. The constructor
   * will set the `edgelength_TH` to 0.002 * diagonal of the bounding box
   * of the given mesh. The other parameters are set to their default values:
   *
   * omega_H = 0.1
   *
   * omega_P = 0.2
   *
   * delta_area = 0.0001
   *
   * max_iterations = 500
   *
   * is_medially_centered = true
   *
   * @param P
   *        Triangulated surface mesh used to extract skeleton.
   */
  MCF_skel_args(HalfedgeGraph& P) : omega_H(0.1), omega_P(0.2),
    delta_area(0.0001), max_iterations(500), is_medially_centered(true)
  {
    vertex_iterator vb, ve;
    boost::tie(vb, ve) = boost::vertices(P);
    Bbox bbox = CGAL::bbox_3(boost::make_transform_iterator(vb, v_to_p),
                             boost::make_transform_iterator(ve, v_to_p));
    double diag = diagonal(bbox);
    edgelength_TH = 0.002 * diag;
  }

  /** Control the velocity of movement and approximation quality.
   *  Increasing `omega_H` will make the MCF based contraction converges
   *  faster, but result in a skeleton of worse quality. */
  double omega_H;
  /** Control the smoothness of the medial approximation.
   *  Increasing `omega_P` will result a skeleton more closely attached
   *  to the medial axis, but slow down the speed of contraction. */
  double omega_P;
  /** Edges with length less than `edgelength_TH` will be collapsed */
  double edgelength_TH;
  /** `run_to_converge` function of `MCF_Skeleton` class will stop if
   *  the change of area in one iteration is less than `delta_area` */
  double delta_area;
  /** Maximum number of iterations. Used to prevent the algorithm running
      too many iterations but still does not satisfy the stopping criteria
      defined by `delta_area` */
  int max_iterations;
  /** If set to true, the result skeleton is medially centered.
      Otherwise set to false. */
  bool is_medially_centered;
};

/// \ingroup PkgMeanCurvatureSkeleton3
/// @brief Class providing the functionalities for extracting
///        the skeleton of a triangulated surface mesh.
///
/// @tparam HalfedgeGraph
///         a model of `HalfedgeGraph`
/// @tparam Graph
///         a model of boost::adjacency_list
///         It is a data structure for the skeleton curve.
/// @tparam VertexIndexMap
///         a model of `ReadWritePropertyMap`</a>
///         with MCF_Skeleton::vertex_descriptor as key and
///         `unsigned int` as value type
/// @tparam EdgeIndexMap
///         a model of `ReadWritePropertyMap`</a>
///         with MCF_Skeleton::edge_descriptor as key and
///         `unsigned int` as value type
/// @tparam GraphCorrelationPMap
///         a model of `ReadWritePropertyMap`/a>
///         with Graph::vertex_descriptor as key and
///         `std::vector<int>` as value type
/// @tparam GraphPointPMap
///         a model of `ReadWritePropertyMap`</a>
///         with Graph::vertex_descriptor as key and
///         MCF_Skeleton::Point as value type
/// @tparam HalfedgeGraphPointPMap
///         a model of `ReadWritePropertyMap`</a>
///         with HalfedgeGraph::vertex_descriptor as key and
///         MCF_Skeleton::Point as value type
/// @tparam SparseLinearAlgebraTraits_d
///         a model of `SparseLinearAlgebraTraitsWithPreFactor_d`
/// @cond CGAL_DOCUMENT_INTERNAL
/// @tparam Collapse_algorithm_tag
///         tag for selecting the edge collapse algorithm
/// @tparam Degeneracy_algorithm_tag
///         tag for selecting the degeneracy detection algorithm
/// @endcond
#ifdef DOXYGEN_RUNNING
template <class HalfedgeGraph,
          class Graph,
          class VertexIndexMap,
          class EdgeIndexMap,
          class GraphCorrelationPMap,
          class HalfedgeGraphPointPMap,
          class GraphPointPMap,
          class SparseLinearAlgebraTraits_d>
#else
template <class HalfedgeGraph,
          class Graph,
          class VertexIndexMap,
          class EdgeIndexMap,
          class GraphCorrelationPMap,
          class GraphPointPMap,
          class HalfedgeGraphPointPMap,
          class SparseLinearAlgebraTraits_d,
          Collapse_algorithm_tag Collapse_tag = LINEAR,
          Degeneracy_algorithm_tag Degeneracy_tag = EULER>
#endif
class MCF_Skeleton
{
// Public types
public:

  // Geometric types
  typedef typename HalfedgeGraph::Traits         Kernel;
  typedef typename Kernel::Vector_3              Vector;
  typedef typename Kernel::Point_3               Point;

  // Repeat HalfedgeGraph types
  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_iterator         vertex_iterator;
  typedef typename HalfedgeGraph::Vertex_handle                                Vertex_handle;

  typedef typename boost::graph_traits<HalfedgeGraph>::edge_descriptor         edge_descriptor;
  typedef typename boost::graph_traits<HalfedgeGraph>::edge_iterator           edge_iterator;
  typedef typename boost::graph_traits<HalfedgeGraph>::in_edge_iterator        in_edge_iterator;
  typedef typename boost::graph_traits<HalfedgeGraph>::out_edge_iterator	     out_edge_iterator;

  typedef typename HalfedgeGraph::Face_handle                                  Face_handle;
  typedef typename HalfedgeGraph::Facet_iterator                               Facet_iterator;
  typedef typename HalfedgeGraph::Halfedge_around_facet_circulator             Halfedge_facet_circulator;

  // Cotangent weight calculator
  typedef typename internal::Cotangent_weight<HalfedgeGraph,
  internal::Cotangent_value_minimum_zero<HalfedgeGraph,
  internal::Cotangent_value_Meyer_secure<HalfedgeGraph> > >                    Weight_calculator;

  typedef internal::Curve_skeleton<HalfedgeGraph, Graph,
  VertexIndexMap, EdgeIndexMap,
  HalfedgeGraphPointPMap, GraphPointPMap>                                      Skeleton;

  // Repeat Graph types
  typedef typename boost::graph_traits<Graph>::vertex_descriptor               vertex_desc;
  typedef typename boost::graph_traits<Graph>::vertex_iterator                 vertex_iter;
  typedef typename boost::graph_traits<Graph>::in_edge_iterator                in_edge_iter;
  typedef typename boost::graph_traits<Graph>::out_edge_iterator               out_edge_iter;
  typedef typename boost::graph_traits<Graph>::edge_iterator                   edge_iter;
  typedef typename boost::graph_traits<Graph>::edge_descriptor                 edge_desc;

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

  /** Source triangulated surface mesh for skeletonization */
  HalfedgeGraph* mesh;
  /** A copy of source mesh.
      All the modifications are operated on it. */
  HalfedgeGraph* polyhedron;

  /** If the object owns the copy of the polyhedron.
   *  If true, the polyhedron will be deleted when the destructor is called. */
  bool owns_polyhedron;

  /** Storing indices of all vertices. */
  VertexIndexMap vertex_id_pmap;
  /** Storing indices of all edges. */
  EdgeIndexMap edge_id_pmap;
  /** Storing the point for HalfedgeGraph vertex_descriptor. */
  HalfedgeGraphPointPMap hg_point_pmap;

  /** Controling the velocity of movement and approximation quality. */
  double omega_H;
  /** Controling the smoothness of the medial approximation. */
  double omega_P;
  /** Edges with length less than `edgelength_TH` will be collapsed. */
  double edgelength_TH;
  /** Triangles with angle greater than `alpha_TH` will be split. */
  double alpha_TH;
  /** Value very close to zero. */
  double zero_TH;
  /** `run_to_converge` will stop if the change of area in one iteration
   *  is less than `delta_area`. */
  double delta_area;
  /** Surface area of original mesh. */
  double original_area;
  /** Maximum number of iterations. */
  int max_iterations;
  /** Should the skeleton be medially centered? */
  bool is_medially_centered;
  /** Are poles computed? */
  bool are_poles_computed;

  /** Cotangent weight calculator. */
  Weight_calculator weight_calculator;
  /** Storing the weights for edges. */
  std::vector<double> edge_weight;
  /** The sparse solver. */
  SparseLinearAlgebraTraits_d m_solver;

  /** Assign a unique id to a new vertex. */
  int vertex_id_count;
  /** The maximum id for original surface. vertices with ids
   *  greater than `max_id` are created during split,
   *  thus will not be considered in correspondence tracking. */
  int max_id;
  /** Used when assembling the matrix. */
  std::map<int, int> new_id;

  /** Store the id of fixed vertices. */
  std::map<size_t, bool> is_vertex_fixed_map;

  /** The incident angle for a halfedge. */
  std::vector<double> halfedge_angle;

  /** Record the correspondence between final surface
   *  and original surface points. */
  std::map<int, std::vector<int> > correspondence;
  /** Record the correspondence between skeletal points
   *  and original surface points. */
  std::map<vertex_desc, std::vector<int> > skeleton_to_surface_map;

  /** Record the corresponding pole of a point. */
  std::map<int, int> poles;
  /** The normal of surface points. */
  std::vector<Vector> normals;
  /** The dual of a cell in Triangulation(a Voronoi point). */
  std::vector<Point> cell_dual;

// Public methods
public:

  /// \name Constructor and Destructor
  /// @{

  /**
   * The constructor of a MCF_Skeleton object.
   *
   * @pre the polyhedron is a watertight triangular mesh
   * @param P
   *        triangulated surface mesh used to extract skeleton
   * \note  The algorithm will make a copy of the source mesh and
   *        keep the source mesh intact.
   * @param Vertex_index_map
   *        property map for associating an id to each vertex
   * @param Edge_index_map
   *        property map for associating an id to each edge
   * @param Skeleton_args
   *        parameters for MCF_Skeleton algorithm
   */
  MCF_Skeleton(HalfedgeGraph& P,
              VertexIndexMap Vertex_index_map,
              EdgeIndexMap Edge_index_map,
              MCF_skel_args<HalfedgeGraph> Skeleton_args
              )
    :mesh(&P), polyhedron(new HalfedgeGraph(P)),
     vertex_id_pmap(Vertex_index_map),
     edge_id_pmap(Edge_index_map),
     hg_point_pmap(boost::get(vertex_point, P))
  {
    owns_polyhedron = true;
    init_args(Skeleton_args);
    init();
  }

  /**
   * The constructor of a MCF_Skeleton object.
   *
   * @pre the polyhedron is a watertight triangular mesh
   * @param P
   *        triangulated surface mesh used to extract skeleton
   * \note  The source mesh will be modified by the algorithm.
   * @param Vertex_index_map
   *        property map for associating an id to each vertex
   * @param Edge_index_map
   *        property map for associating an id to each edge
   * @param Skeleton_args
   *        parameters for MCF_Skeleton algorithm
   */
  MCF_Skeleton(HalfedgeGraph* P,
              VertexIndexMap Vertex_index_map,
              EdgeIndexMap Edge_index_map,
              MCF_skel_args<HalfedgeGraph> Skeleton_args
              )
    :mesh(NULL), polyhedron(P),
     vertex_id_pmap(Vertex_index_map),
     edge_id_pmap(Edge_index_map),
     hg_point_pmap(boost::get(vertex_point, *P))
  {
    owns_polyhedron = false;
    init_args(Skeleton_args);
    init();
  }

  ~MCF_Skeleton()
  {
    if (owns_polyhedron)
    {
      delete polyhedron;
    }
  }

  /// @} Constructor and Destructor

  /// \name Setter and Getter
  /// @{

  void set_omega_H(double value)
  {
    omega_H = value;
  }

  double get_omega_H()
  {
    return omega_H;
  }

  void set_omega_P(double value)
  {
    omega_P = value;
  }

  double get_omega_P()
  {
    return omega_P;
  }

  void set_edgelength_TH(double value)
  {
    edgelength_TH = value;
  }

  double get_edgelength_TH()
  {
    return edgelength_TH;
  }

  void set_delta_area(double value)
  {
    delta_area = value;
  }

  double get_delta_area()
  {
    return delta_area;
  }

  void set_is_medially_centered(bool value)
  {
    is_medially_centered = value;
  }

  bool get_is_medially_centered()
  {
    return is_medially_centered;
  }

  void set_max_iterations(int value)
  {
    max_iterations = value;
  }

  int get_max_iterations()
  {
    return max_iterations;
  }

  /**
   * Get the pointer to the contracted mesh.
   */
  HalfedgeGraph* get_halfedge_graph()
  {
    return polyhedron;
  }

  /**
   * Get the pointer to the copy of source mesh.
   */
  HalfedgeGraph* get_mesh()
  {
    return mesh;
  }

  /// \cond SKIP_FROM_MANUAL

  void set_alpha_TH(double value)
  {
    alpha_TH = value;
  }

  double get_alpha_TH()
  {
    return alpha_TH;
  }

  void set_zero_TH(double value)
  {
    zero_TH = value;
  }

  double get_zero_TH()
  {
    return zero_TH;
  }

  /**
   * @brief If set to true, the copy of the source mesh will be deleted when
   *        the destructor is called.
   */
  void set_own_polyhedron(bool value)
  {
    owns_polyhedron = value;
  }

  bool get_own_polyhedron()
  {
    return owns_polyhedron;
  }

  /// \endcond

  /// @cond CGAL_DOCUMENT_INTERNAL

  /**
   * Get the positions of fixed(degenerate) points.
   *
   * @param fixed_points
   *        return the positions of fixed points
   */
  void get_fixed_points(std::vector<Point>& fixed_points)
  {
    fixed_points.clear();
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end())
      {
        vertex_descriptor vd = *vb;
        fixed_points.push_back(boost::get(hg_point_pmap, vd));
      }
    }
  }

  /**
   * Get the positions of non-fixed(non-degenerate) points.
   *
   * @param non_fixed_points
   *        return the positions of non-fixed points
   */
  void get_non_fixed_points(std::vector<Point>& non_fixed_points)
  {
    non_fixed_points.clear();
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      if (is_vertex_fixed_map.find(id) == is_vertex_fixed_map.end())
      {
          vertex_descriptor vd = *vb;
          non_fixed_points.push_back(boost::get(hg_point_pmap, vd));
      }
    }
  }

  /**
   * Get the Voronoi pole for the surface mesh.
   *
   * @param max_poles
   *        for each mesh vertex, record its correspondent Voronoi pole position
   */
  void get_poles(std::vector<Point>& max_poles)
  {
    max_poles.resize(boost::num_vertices(*polyhedron));
    vertex_iterator vb, ve;
    int cnt = 0;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int vid = boost::get(vertex_id_pmap, v);
      max_poles[cnt++] = cell_dual[poles[vid]];
    }
  }

  /// @endcond

  /// @} Setter and Getter


  /// \name Public Algorithm API
  /// @{

  /**
   * Extract the skeleton curve for the mesh. The algorithm repeatedly
   * contract the mesh until convergence, then turn the contracted mesh
   * to a curve skeleton.
   *
   * @param g
   *        a boost::graph containing the connectivity of the skeleton
   * @param points
   *        the embedding of the skeletal points
   */
  void extract_skeleton(Graph& g, GraphPointPMap& points)
  {
    run_to_converge();
    convert_to_skeleton(g, points);
  }

  /**
   * Extract the skeleton curve for the mesh. The algorithm repeatedly
   * contract the mesh until convergence, then turn the contracted mesh
   * to a curve skeleton.
   *
   * @param g
   *        a boost::graph containing the connectivity of the skeleton
   * @param points
   *        the locations of the skeletal points
   * @param skeleton_to_surface
   *        for each skeletal point, record its correspondent surface points
   */
  void extract_skeleton(Graph& g, GraphPointPMap& points,
                        GraphCorrelationPMap& skeleton_to_surface)
  {
    run_to_converge();
    convert_to_skeleton(g, points);
    get_correspondent_vertices(skeleton_to_surface);
  }

  /**
   * Contract the mesh by mean curvature flow.
   */
  void contract_geometry()
  {
    MCFSKEL_DEBUG(std::cerr << "before contract geometry";)

    update_vertex_id();

    compute_edge_weight();

    int nver = boost::num_vertices(*polyhedron);
    int nrows;
    if (is_medially_centered)
    {
      nrows = nver * 3;
      if (!are_poles_computed)
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

    // copy to mesh
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor vi = *vb;
      int id = boost::get(vertex_id_pmap, vi);
      int i = new_id[id];
      Point p(X[i], Y[i], Z[i]);
      boost::put(hg_point_pmap, vi, p);
    }

    MCFSKEL_DEBUG(std::cerr << "leave contract geometry\n";)
  }

  /**
   * Collapse short edges with length less than `edgelength_TH`.
   */
  int collapse_edges()
  {
    internal::Fixed_edge_map<HalfedgeGraph> fixed_edge_map;
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
   * Split triangles with one angle greater than `alpha_TH`.
   */
  int split_triangles()
  {
    MCFSKEL_DEBUG(std::cerr << "before split\n";)

    int num_splits = 0;
    while (true)
    {
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

  /**
   * Run a combination of `collapse_edges` and `split_triangles`.
   */
  int update_topology()
  {
    MCFSKEL_DEBUG(std::cerr << "before collapse edges\n";)

    int num_collapses = collapse_edges();
    MCFSKEL_INFO(std::cerr << "collapse " << num_collapses << " edges.\n";)

    int num_splits = split_triangles();
    MCFSKEL_INFO(std::cerr << "split " << num_splits << " edges.\n";)

    return num_collapses + num_splits;
  }

  /**
   * Fix degenerate vertices.
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

  /**
   * Run an iteration of `contract_geometry`, `update_topology` and
   * `detect_degeneracies`.
   */
  void contract()
  {
    contract_geometry();
    update_topology();
    detect_degeneracies();

    MCFSKEL_INFO(double area = internal::get_surface_area(*polyhedron, hg_point_pmap);)
    MCFSKEL_INFO(std::cout << "area " << area << "\n";)
  }

  /**
   * Run iterations of `contract_geometry()`, `update_topology()` and
   * `detect_degeneracies()` until the change of surface area during one
   * iteration is less than `delta_area` * original surface area.
   */
  void run_to_converge()
  {
    double last_area = 0;
    int num_iteration = 0;
    while (true)
    {
      MCFSKEL_INFO(std::cout << "iteration " << num_iteration + 1 << "\n";)

      contract_geometry();
      update_topology();
      detect_degeneracies();

      double area = internal::get_surface_area(*polyhedron, hg_point_pmap);
      double area_ratio = fabs(last_area - area) / original_area;

      MCFSKEL_INFO(std::cout << "area " << area << "\n";)
      MCFSKEL_INFO(std::cout << "|area - last_area| / original_area "
                             << area_ratio << "\n";)

      if (area_ratio < delta_area)
      {
        break;
      }
      last_area = area;

      num_iteration++;
      if (num_iteration >= max_iterations)
      {
        break;
      }
    }
  }

  /**
   * Convert the contracted mesh to a skeleton curve.
   * @param g
   *        a boost::graph containing the connectivity of the skeleton
   * @param points
   *        the locations of the skeletal points
   */
  void convert_to_skeleton(Graph& g, GraphPointPMap& points)
  {
    Skeleton skeleton(*polyhedron, vertex_id_pmap, edge_id_pmap, hg_point_pmap);

    skeleton.extract_skeleton(g, points, skeleton_to_surface_map);
  }

  /**
   * Get the correspondent surface points for the skeleton.
   *
   * @param skeleton_to_surface
   *        for each skeletal point, record its correspondent surface points
   */
  void get_correspondent_vertices(GraphCorrelationPMap& skeleton_to_surface)
  {
    typename std::map<vertex_desc, std::vector<int> >::iterator iter;
    for (iter = skeleton_to_surface_map.begin();
         iter != skeleton_to_surface_map.end(); ++iter)
    {
      vertex_desc i = iter->first;

      skeleton_to_surface[i] = std::vector<int>();
      for (size_t j = 0; j < skeleton_to_surface_map[i].size(); ++j)
      {
        int id = skeleton_to_surface_map[i][j];
        if (correspondence.find(id) != correspondence.end())
        {
          skeleton_to_surface[i].insert(skeleton_to_surface[i].end(),
                                        correspondence[id].begin(),
                                        correspondence[id].end());
        }

        if (id < max_id)
        {
          skeleton_to_surface[i].push_back(id);
        }
      }
    }
  }

  /// @} Public Algorithm API


private:

  // --------------------------------------------------------------------------
  // Initialization
  // --------------------------------------------------------------------------

  /// Initialize the parameters for MCF_Skeleton
  void init_args(MCF_skel_args<HalfedgeGraph> Skeleton_args)
  {
    omega_H = Skeleton_args.omega_H;
    omega_P = Skeleton_args.omega_P;
    edgelength_TH = Skeleton_args.edgelength_TH;
    delta_area = Skeleton_args.delta_area;
    max_iterations = Skeleton_args.max_iterations;
    is_medially_centered = Skeleton_args.is_medially_centered;
    weight_calculator = Weight_calculator();
    alpha_TH = 110;
    zero_TH = 1e-7;
  }

  /// Initialize some global data structures such as vertex id.
  void init()
  {
    are_poles_computed = false;

    vertex_iterator vb, ve;

    alpha_TH *= (M_PI / 180.0);
    double area = internal::get_surface_area(*polyhedron, hg_point_pmap);
    original_area = area;

    // initialize index maps
    vertex_id_count = 0;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      boost::put(vertex_id_pmap, *vb, vertex_id_count++);
    }
    vertex_id_count = 0;
    for (boost::tie(vb, ve) = boost::vertices(*mesh); vb != ve; ++vb)
    {
      boost::put(vertex_id_pmap, *vb, vertex_id_count++);
    }

    max_id = vertex_id_count;

    edge_iterator eb, ee;
    int idx = 0;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, idx++);
    }

    is_vertex_fixed_map.clear();
    correspondence.clear();

    if (is_medially_centered)
    {
      compute_voronoi_pole();
    }
  }

  // --------------------------------------------------------------------------
  // Contraction
  // --------------------------------------------------------------------------

  /// Compute cotangent weights of all edges.
  void compute_edge_weight()
  {
    edge_weight.clear();
    edge_weight.reserve(boost::num_edges(*polyhedron));
    edge_iterator eb, ee;
    for(boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      edge_weight.push_back(this->weight_calculator(*eb, *polyhedron));
    }
  }

  /// Assemble the left hand side.
  void assemble_LHS(typename SparseLinearAlgebraTraits_d::Matrix& A)
  {
    MCFSKEL_DEBUG(std::cerr << "start LHS\n";)

    int nver = boost::num_vertices(*polyhedron);

    Point_inside_polyhedron_3<HalfedgeGraph, Kernel> test_inside(*polyhedron);

    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);

      int i = new_id[id];
      // if the vertex is fixed
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end())
      {
        A.set_coef(i + nver, i, 1.0 / zero_TH, true);
      }
      else
      {
        A.set_coef(i + nver, i, omega_H, true);
        if (is_medially_centered)
        {
          if (id < max_id)
          {
            if (test_inside(cell_dual[poles[id]]) == CGAL::ON_BOUNDED_SIDE)
            {
              A.set_coef(i + nver * 2, i, omega_P, true);
            }
          }
        }
      }
    }

    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      int i = new_id[id];
      double L = 1.0;
      // if the vertex is fixed
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end())
      {
        L = 0;
      }
      double diagonal = 0;
      in_edge_iterator e, e_end;
      for (boost::tie(e, e_end) = boost::in_edges(*vb, *polyhedron); e != e_end; ++e)
      {
        vertex_descriptor vj = boost::source(*e, *polyhedron);
        double wij = edge_weight[boost::get(edge_id_pmap, *e)] * 2.0;
        int jd = boost::get(vertex_id_pmap, vj);
        int j = new_id[jd];
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

    Point_inside_polyhedron_3<HalfedgeGraph, Kernel> test_inside(*polyhedron);

    // assemble right columns of linear system
    int nver = boost::num_vertices(*polyhedron);
    vertex_iterator vb, ve;
    for (int i = 0; i < nver; ++i)
    {
      Bx[i] = 0;
      By[i] = 0;
      Bz[i] = 0;
    }

    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor vi = *vb;
      int id = boost::get(vertex_id_pmap, vi);
      int i = new_id[id];

      double oh, op = 0.0;
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end())
      {
        oh = 1.0 / zero_TH;
      }
      else
      {
        oh = omega_H;
        if (is_medially_centered)
        {
          if (id < max_id)
          {
            if (test_inside(cell_dual[poles[id]]) == CGAL::ON_BOUNDED_SIDE)
            {
              op = omega_P;
            }
          }
        }
      }
      Bx[i + nver] = boost::get(hg_point_pmap, vi).x() * oh;
      By[i + nver] = boost::get(hg_point_pmap, vi).y() * oh;
      Bz[i + nver] = boost::get(hg_point_pmap, vi).z() * oh;
      if (is_medially_centered)
      {
        double x = to_double(cell_dual[poles[id]].x());
        double y = to_double(cell_dual[poles[id]].y());
        double z = to_double(cell_dual[poles[id]].z());
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
    new_id.clear();
    int cnt = 0;
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      new_id[id] = cnt++;
    }
  }

  // --------------------------------------------------------------------------
  // Edge collapse
  // --------------------------------------------------------------------------

  /// Collapse short edges using simplification package.
  int collapse_edges_simplification()
  {
    internal::Fixed_edge_map<HalfedgeGraph> fixed_edge_map;

    edge_iterator eb, ee;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      vertex_descriptor vi = boost::source(*eb, *polyhedron);
      vertex_descriptor vj = boost::target(*eb, *polyhedron);
      size_t vi_idx = boost::get(vertex_id_pmap, vi);
      size_t vj_idx = boost::get(vertex_id_pmap, vj);

      if (is_vertex_fixed_map.find(vi_idx) != is_vertex_fixed_map.end()
       && is_vertex_fixed_map.find(vj_idx) != is_vertex_fixed_map.end())
      {
        fixed_edge_map.set_is_fixed(*eb, true);
      }
    }

    int edge_id = 0;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, edge_id++);
    }

    // This is a stop predicate (defines when the algorithm terminates).
    // The simplification stops when the length of all edges is greater
    // than the minimum threshold.
    CGAL::internal::Minimum_length_predicate<HalfedgeGraph> stop(edgelength_TH);

    // midpoint placement without geometric test
    SMS::Geometric_test_skipper< SMS::Midpoint_placement<HalfedgeGraph> > placement;

    internal::Track_correspondence_visitor<HalfedgeGraph, HalfedgeGraphPointPMap> vis;
    if (is_medially_centered)
    {
      vis = internal::Track_correspondence_visitor<HalfedgeGraph, HalfedgeGraphPointPMap>
            (&hg_point_pmap, &correspondence, &poles, &cell_dual, max_id);
    }
    else
    {
      vis = internal::Track_correspondence_visitor<HalfedgeGraph, HalfedgeGraphPointPMap>
            (&hg_point_pmap, &correspondence, max_id);
    }

    int r = SMS::edge_collapse
                (*polyhedron
                ,stop
                ,CGAL::get_cost(SMS::Edge_length_cost<HalfedgeGraph>())
                      .get_placement(placement)
                      .visitor(vis)
                      .edge_is_border_map(fixed_edge_map)
                );

    return r;
  }

  /// Track correspondent original surface points during collapse.
  void track_correspondence(vertex_descriptor v0, vertex_descriptor v1,
                            vertex_descriptor v)
  {
    int id0 = v0->id();
    int id1 = v1->id();
    int vid = v->id();
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

    if (correspondence.find(to) == correspondence.end())
    {
      correspondence[to] = std::vector<int>();
    }
    // only track vertex in original mesh
    if (from < max_id)
    {
      correspondence[to].push_back(from);
    }
    std::map<int, std::vector<int> >::iterator iter = correspondence.find(from);
    if (iter != correspondence.end())
    {
      for (size_t i = 0; i < (iter->second).size(); ++i)
      {
        correspondence[to].push_back((iter->second)[i]);
      }
      (iter->second).clear();
      correspondence.erase(iter);
    }

    if (is_medially_centered)
    {
      Point pole0 = Point(to_double(cell_dual[poles[id0]].x()),
                          to_double(cell_dual[poles[id0]].y()),
                          to_double(cell_dual[poles[id0]].z()));
      Point pole1 = Point(to_double(cell_dual[poles[id1]].x()),
                          to_double(cell_dual[poles[id1]].y()),
                          to_double(cell_dual[poles[id1]].z()));
      Point p1 = boost::get(hg_point_pmap, v1);
      double dis_to_pole0 = sqrt(squared_distance(pole0, p1));
      double dis_to_pole1 = sqrt(squared_distance(pole1, p1));
      if (dis_to_pole0 < dis_to_pole1)
      {
        poles[id1] = poles[id0];
      }
      std::map<int, int>::iterator pole_iter = poles.find(id0);
      poles.erase(pole_iter);
    }
  }

  /// Collapse short edges by iteratively linear search.
  int collapse_edges_linear(internal::Fixed_edge_map<HalfedgeGraph>& fixed_edge_map)
  {
    std::vector<edge_descriptor> edges;
    edges.reserve(boost::num_edges(*polyhedron));
    edge_iterator eb, ee;

    boost::tie(eb, ee) = boost::edges(*polyhedron);
    std::copy(eb, ee, std::back_inserter(edges));

    int cnt = 0;
    for (size_t i = 0; i < edges.size(); ++i)
    {
      edge_descriptor ed = edges[i];
      if (fixed_edge_map.is_fixed(ed))
      {
        continue;
      }

      vertex_descriptor vi = boost::source(ed, *polyhedron);
      vertex_descriptor vj = boost::target(ed, *polyhedron);
      double edge_length = sqrt(squared_distance(boost::get(hg_point_pmap, vi),
                                                 boost::get(hg_point_pmap, vj)));
      if (internal::is_collapse_ok(*polyhedron, ed) && edge_length < edgelength_TH)
      {
        Point p = midpoint(
          boost::get(vertex_point, *polyhedron, boost::source(ed, *polyhedron)),
          boost::get(vertex_point, *polyhedron, boost::target(ed, *polyhedron)));

        // invalidate the edges that will be collapsed
        // since the mesh is closed, 6 halfedges will be collapsed
        fixed_edge_map.set_is_fixed(ed, true);
        fixed_edge_map.set_is_fixed(ed->opposite(), true);
        fixed_edge_map.set_is_fixed(ed->prev(), true);
        fixed_edge_map.set_is_fixed(ed->prev()->opposite(), true);
        fixed_edge_map.set_is_fixed(ed->opposite()->prev(), true);
        fixed_edge_map.set_is_fixed(ed->opposite()->prev()->opposite(), true);

        vertex_descriptor v = SMS::halfedge_collapse(ed, *polyhedron);
        boost::put(vertex_point, *polyhedron, v, p);

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
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      vertex_descriptor vi = boost::source(*eb, *polyhedron);
      vertex_descriptor vj = boost::target(*eb, *polyhedron);
      size_t vi_idx = boost::get(vertex_id_pmap, vi);
      size_t vj_idx = boost::get(vertex_id_pmap, vj);

      if (is_vertex_fixed_map.find(vi_idx) != is_vertex_fixed_map.end()
       && is_vertex_fixed_map.find(vj_idx) != is_vertex_fixed_map.end())
      {
        fixed_edge_map.set_is_fixed(*eb, true);
      }
    }
  }

  // --------------------------------------------------------------------------
  // Triangle split
  // --------------------------------------------------------------------------

  /// Compute the incident angles for all the halfedges.
  void compute_incident_angle()
  {
    halfedge_angle.clear();
    int ne = boost::num_edges(*polyhedron);
    halfedge_angle.resize(ne, 0);

    edge_iterator eb, ee;
    int idx = 0;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, idx++);
    }

    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      int e_id = boost::get(edge_id_pmap, *eb);
      edge_descriptor ed = *eb;

      if (ed->is_border())
      {
        halfedge_angle[e_id] = -1;
      }
      else
      {
        vertex_descriptor vi = boost::source(ed, *polyhedron);
        vertex_descriptor vj = boost::target(ed, *polyhedron);
        edge_descriptor ed_next = ed->next();
        vertex_descriptor vk = boost::target(ed_next, *polyhedron);
        Point pi = boost::get(hg_point_pmap, vi);
        Point pj = boost::get(hg_point_pmap, vj);
        Point pk = boost::get(hg_point_pmap, vk);

        double dis2_ij = squared_distance(pi, pj);
        double dis2_ik = squared_distance(pi, pk);
        double dis2_jk = squared_distance(pj, pk);
        double dis_ij = sqrt(dis2_ij);
        double dis_ik = sqrt(dis2_ik);
        double dis_jk = sqrt(dis2_jk);

        // A degenerate triangle will never undergo a split (but rather a collapse...)
        if (dis_ij < zero_TH || dis_ik < zero_TH || dis_jk < zero_TH)
        {
          halfedge_angle[e_id] = -1;
        }
        else
        {
          halfedge_angle[e_id] =
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
    Point ps = boost::get(hg_point_pmap, vs);
    Point pt = boost::get(hg_point_pmap, vt);
    Point pk = boost::get(hg_point_pmap, vk);
    CGAL::internal::Vector vec_st = CGAL::internal::Vector(ps, pt);
    CGAL::internal::Vector vec_sk = CGAL::internal::Vector(ps, pk);

    vec_st.normalize();
    double t = vec_st.dot(vec_sk);
    Point st = Point(vec_st[0] * t, vec_st[1] * t, vec_st[2] * t);
    Point pn = Point(ps[0] + st[0], ps[1] + st[1], ps[2] + st[2]);

    // project the pole
    if (is_medially_centered)
    {
      int sid = boost::get(vertex_id_pmap, vs);
      int tid = boost::get(vertex_id_pmap, vt);
      Point pole_s = cell_dual[poles[sid]];
      Point pole_t = cell_dual[poles[tid]];
      Vector pole_st = pole_t - pole_s;
      Vector p_projector = pole_st / sqrt(pole_st.squared_length());
      Point pole_n = pole_s + p_projector * t;
      poles[vertex_id_count] = cell_dual.size();
      cell_dual.push_back(pole_n);
    }
    return pn;
  }

  /// Split triangles with an angle greater than `alpha_TH`.
  int split_flat_triangle()
  {
    int ne = boost::num_edges(*polyhedron);
    compute_incident_angle();

    int cnt = 0;
    edge_iterator eb, ee;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      edge_descriptor ei = *eb;
      edge_descriptor ej = ei->opposite();
      int ei_id = boost::get(edge_id_pmap, ei);
      int ej_id = boost::get(edge_id_pmap, ej);
      if (ei_id < 0 || ei_id >= ne
       || ej_id < 0 || ej_id >= ne)
      {
        continue;
      }

      vertex_descriptor vs = boost::source(ei, *polyhedron);
      vertex_descriptor vt = boost::target(ei, *polyhedron);

      double angle_i = halfedge_angle[ei_id];
      double angle_j = halfedge_angle[ej_id];
      if (angle_i < alpha_TH || angle_j < alpha_TH)
      {
        continue;
      }

      edge_descriptor ek;
      if (angle_i > angle_j)
      {
        ek = ei->next();
      }
      else
      {
        ek = ej->next();
      }
      vertex_descriptor vk = boost::target(ek, *polyhedron);
      Point pn = project_vertex(vs, vt, vk);
      edge_descriptor en = internal::mesh_split(*polyhedron, hg_point_pmap, ei, pn);
      // set id for new vertex
      boost::put(vertex_id_pmap, en->vertex(), vertex_id_count++);
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
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int idx = boost::get(vertex_id_pmap, v);

      if (is_vertex_fixed_map.find(idx) == is_vertex_fixed_map.end())
      {
        bool willbefixed = internal::is_vertex_degenerate(*polyhedron, hg_point_pmap,
                                                          v, edgelength_TH);
        if (willbefixed)
        {
          is_vertex_fixed_map[idx] = willbefixed;
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
    double elength_fixed = edgelength_TH;
    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int idx = boost::get(vertex_id_pmap, v);
      if (is_vertex_fixed_map.find(idx) == is_vertex_fixed_map.end())
      {
        bool willbefixed = false;
        int bad_counter = 0;

        in_edge_iterator eb, ee;
        for (boost::tie(eb, ee) = boost::in_edges(v, *polyhedron); eb != ee; ++eb)
        {
          edge_descriptor edge = *eb;
          vertex_descriptor v0 = boost::source(edge, *polyhedron);
          vertex_descriptor v1 = boost::target(edge, *polyhedron);
          double length = sqrt(squared_distance(boost::get(hg_point_pmap, v0),
                                                boost::get(hg_point_pmap, v1)));
          if (length < elength_fixed)
          {
            if (!internal::is_collapse_ok(*polyhedron, edge))
            {
              bad_counter++;
            }
          }
        }
        willbefixed = (bad_counter >= 2);
        if (willbefixed)
        {
          is_vertex_fixed_map[idx] = willbefixed;
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
    cell_dual.clear();
    point_to_pole.clear();
    point_to_pole.resize(boost::num_vertices(*polyhedron));

    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int vid = boost::get(vertex_id_pmap, v);
      Exact_point tp((boost::get(hg_point_pmap, v)).x(),
                     (boost::get(hg_point_pmap, v)).y(),
                     (boost::get(hg_point_pmap, v)).z());
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
      cell_dual.push_back(pt);
      // each cell has 4 incident vertices
      for (int i = 0; i < 4; ++i)
      {
        TriVertex_handle vt = cell->vertex(i);
        int id = vt->info();
        point_to_pole[id].push_back(cell_id);
      }
      cell_id++;
    }

    poles.clear();
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
        Point cell_point = cell_dual[pole_id];
        Vector vt = cell_point - surface_point;
        Vector n = normals[i];

        double t = vt * n;

        // choose the one with maximum distance along the normal
        if (t < 0 && t < max_neg_t)
        {
          max_neg_i = pole_id;
          max_neg_t = t;
        }
      }

      poles[i] = max_neg_i;
    }
    are_poles_computed = true;
  }

  /// Compute an approximate vertex normal for all vertices.
  void compute_vertex_normal()
  {
    normals.resize(boost::num_vertices(*polyhedron));

    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int vid = boost::get(vertex_id_pmap, v);
      normals[vid] = internal::get_vertex_normal<typename HalfedgeGraph::Vertex,Kernel>(*v);
    }
  }

};

/// \ingroup PkgMeanCurvatureSkeleton3
/// @brief Extract a medially centered curve skeleton for the mesh.
///
/// @pre the polyhedron is a watertight triangular mesh
///
/// @tparam HalfedgeGraph
///         a model of `HalfedgeGraph`
/// @tparam Graph
///         a model of boost::adjacency_list
///         data structure for skeleton curve
/// @tparam VertexIndexMap
///         a model of `ReadWritePropertyMap`</a>
///         with MCF_Skeleton::vertex_descriptor as key and
///         `unsigned int` as value type
/// @tparam EdgeIndexMap
///         a model of `ReadWritePropertyMap`</a>
///         with MCF_Skeleton::edge_descriptor as key and
///         `unsigned int` as value type
/// @tparam GraphCorrelationPMap
///         a model of `ReadWritePropertyMap`</a>
///         with Graph::vertex_descriptor as key and
///         `std::vector<int>` as value type
/// @tparam HalfedgeGraphPointPMap
///         a model of `ReadWritePropertyMap`</a>
///         with HalfedgeGraph::vertex_descriptor as key and
///         MCF_Skeleton::Point as value type
/// @tparam GraphPointPMap
///         a model of `ReadWritePropertyMap`</a>
///         with Graph::vertex_descriptor as key and
///         MCF_Skeleton::Point as value type
/// @tparam SparseLinearAlgebraTraits_d
///         a model of `SparseLinearAlgebraTraitsWithPreFactor_d`
///
/// @param P
///        triangulated surface mesh used to extract skeleton
/// @param Vertex_index_map
///        property map for associating an id to each vertex
/// @param Edge_index_map
///        property map for associating an id to each edge
/// @param Skeleton_args
///        parameters for MCF_Skeleton algorithm
/// @param g
///        a boost::graph containing the connectivity of the skeleton
/// @param points
///        the locations of the skeletal points
/// @param skeleton_to_surface
///        for each skeletal point, record its correspondent surface points
template <class HalfedgeGraph,
          class Graph,
          class VertexIndexMap,
          class EdgeIndexMap,
          class GraphCorrelationPMap,
          class GraphPointPMap,
          class HalfedgeGraphPointPMap,
          class SparseLinearAlgebraTraits_d>
void extract_skeleton(HalfedgeGraph& P,
                      VertexIndexMap Vertex_index_map,
                      EdgeIndexMap Edge_index_map,
                      MCF_skel_args<HalfedgeGraph> Skeleton_args,
                      Graph& g, GraphPointPMap& points,
                      GraphCorrelationPMap& skeleton_to_surface)
{
  typedef CGAL::MCF_Skeleton<HalfedgeGraph, Graph, VertexIndexMap, EdgeIndexMap,
  GraphCorrelationPMap, GraphPointPMap, HalfedgeGraphPointPMap, SparseLinearAlgebraTraits_d> MCFSKEL;

  MCFSKEL mcs(P, Vertex_index_map, Edge_index_map, Skeleton_args);

  mcs.run_to_converge();
  mcs.convert_to_skeleton(g, points);

  mcs.get_correspondent_vertices(skeleton_to_surface);
}

template <class HalfedgeGraph,
          class Graph,
          class VertexIndexMap,
          class EdgeIndexMap,
          class GraphCorrelationPMap,
          class GraphPointPMap,
          class HalfedgeGraphPointPMap>
void extract_skeleton(HalfedgeGraph& P,
                      VertexIndexMap Vertex_index_map,
                      EdgeIndexMap Edge_index_map,
                      MCF_skel_args<HalfedgeGraph> Skeleton_args,
                      Graph& g, GraphPointPMap& points,
                      GraphCorrelationPMap& skeleton_to_surface)
{
  typedef CGAL::MCF_default_solver<double>::type Sparse_linear_solver;

  extract_skeleton<HalfedgeGraph, Graph, VertexIndexMap, EdgeIndexMap,
                   GraphCorrelationPMap, GraphPointPMap, HalfedgeGraphPointPMap,
                   Sparse_linear_solver>
      (P, Vertex_index_map, Edge_index_map, Skeleton_args, g, points, skeleton_to_surface);
}

template <class HalfedgeGraph,
          class Graph,
          class VertexIndexMap,
          class EdgeIndexMap,
          class GraphCorrelationPMap,
          class GraphPointPMap>
void extract_skeleton(HalfedgeGraph& P,
                      VertexIndexMap Vertex_index_map,
                      EdgeIndexMap Edge_index_map,
                      MCF_skel_args<HalfedgeGraph> Skeleton_args,
                      Graph& g, GraphPointPMap& points,
                      GraphCorrelationPMap& skeleton_to_surface)
{
  typedef typename boost::property_map<HalfedgeGraph, CGAL::vertex_point_t>::type HalfedgeGraphPointPMap;
  typedef CGAL::MCF_default_solver<double>::type Sparse_linear_solver;

  extract_skeleton<HalfedgeGraph, Graph, VertexIndexMap, EdgeIndexMap,
                   GraphCorrelationPMap, GraphPointPMap, HalfedgeGraphPointPMap,
                   Sparse_linear_solver>
      (P, Vertex_index_map, Edge_index_map, Skeleton_args, g, points, skeleton_to_surface);
}

} //namespace CGAL

#endif // CGAL_MEAN_CURVATURE_SKELETON_H
