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
 * @brief The class `Mean_curvature_skeleton` containing the API to extract
 * curve skeleton for a closed triangular mesh.
 */

#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/graph/copy.hpp>

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

#include <queue>

namespace SMS = CGAL::Surface_mesh_simplification;

namespace CGAL {

/// @cond CGAL_DOCUMENT_INTERNAL

/// \ingroup PkgMeanCurvatureSkeleton3
///@brief Edge collapse algorithm tag
enum Collapse_algorithm_tag
{
  SIMPLIFICATION, /**< algorithm from simplification package */
  LINEAR          /**< iterative linear search */
};

/// \ingroup PkgMeanCurvatureSkeleton3
///@brief Degeneracy detection algorithm tag
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
/// @tparam HalfedgeGraph
///         a model of `HalfedgeGraph`
/// @tparam SparseLinearAlgebraTraits_d
///         a model of `SparseLinearAlgebraTraitsWithPreFactor_d`
/// @tparam VertexIndexMap
///         a model of `ReadWritePropertyMap`</a>
///         with Mean_curvature_skeleton::vertex_descriptor as key and
///         `unsigned int` as value type
/// @tparam EdgeIndexMap
///         a model of `ReadWritePropertyMap`</a>
///         with Mean_curvature_skeleton::edge_descriptor as key and
///         `unsigned int` as value type
/// @tparam Graph
///         a model of boost::adjacency_list
///         data structure for skeleton curve
/// @cond CGAL_DOCUMENT_INTERNAL
/// @tparam Collapse_algorithm_tag
///         tag for selecting the edge collapse algorithm
/// @tparam Degeneracy_algorithm_tag
///         tag for selecting the degeneracy detection algorithm
/// @endcond
#ifdef DOXYGEN_RUNNING
template <class HalfedgeGraph,
          class SparseLinearAlgebraTraits_d,
          class VertexIndexMap,
          class EdgeIndexMap,
          class Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> >
#else
template <class HalfedgeGraph,
          class SparseLinearAlgebraTraits_d,
          class VertexIndexMap,
          class EdgeIndexMap,
          class Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>,
          Collapse_algorithm_tag Collapse_tag = LINEAR,
          Degeneracy_algorithm_tag Degeneracy_tag = EULER>
#endif
class Mean_curvature_skeleton
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
  VertexIndexMap, EdgeIndexMap>                                                Skeleton;

  // Repeat Graph types
  typedef typename boost::graph_traits<Graph>::vertex_descriptor               vertex_desc;
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

  /** source triangulated surface mesh for skeletonization */
  HalfedgeGraph& polyhedron;
  /** storing indices of all vertices */
  VertexIndexMap vertex_id_pmap;
  /** storing indices of all edges */
  EdgeIndexMap edge_id_pmap;

  /** controls the velocity of movement and approximation quality */
  double omega_H;
  /** controls the smoothness of the medial approximation */
  double omega_P;
  /** edges with length less than `edgelength_TH` will be collapsed */
  double edgelength_TH;
  /** triangles with angle greater than `alpha_TH` will be split */
  double alpha_TH;
  /** value very close to zero */
  double zero_TH;
  /** run_to_converge will stop if the change of area in one iteration
   *  is less than `area_TH` */
  double area_TH;
  /** surface area of original mesh */
  double original_area;
  /** maximum number of iterations */
  int iteration_TH;

  /** cotangent weight calculator */
  Weight_calculator weight_calculator;
  /** storing the weights for edges */
  std::vector<double> edge_weight;
  /** the sparse solver */
  SparseLinearAlgebraTraits_d m_solver;

  /** assign a unique id to a new vertex */
  int vertex_id_count;
  /** the maximum id for original surface. vertices with ids
   *  greater than `max_id` are created during split,
   *  thus will not be considered in correspondence tracking */
  int max_id;
  /** used when assembling the matrix */
  std::map<int, int> new_id;

  /** store the id of fixed vertices */
  std::map<size_t, bool> is_vertex_fixed_map;

  /** the incident angle for a halfedge */
  std::vector<double> halfedge_angle;

  /** record the correspondence between final surface
   *  and original surface points */
  std::map<int, std::vector<int> > correspondence;
  /** record the correspondence between skeletal points
   *  and original surface points */
  std::map<vertex_desc, std::vector<int> > skeleton_to_surface;

  /** should the skeleton be medially centered? */
  bool is_medially_centered;
  /** record the corresponding pole of a point */
  std::map<int, int> poles;
  /** the normal of surface points */
  std::vector<Vector> normals;
  /** the dual of a cell in Triangulation(a Voronoi point) */
  std::vector<Point> cell_dual;

// Public methods
public:

  /// \name Constructor and Destructor
  /// @{

  /// \cond SKIP_FROM_MANUAL
  Mean_curvature_skeleton(HalfedgeGraph& P,
                          VertexIndexMap Vertex_index_map,
                          EdgeIndexMap Edge_index_map,
                          double omega_H,
                          double edgelength_TH,
                          double area_TH = 1e-5,
                          int iteration_TH = 500
                          )
    :polyhedron(P),
     vertex_id_pmap(Vertex_index_map),
     edge_id_pmap(Edge_index_map),
     omega_H(omega_H),
     edgelength_TH(edgelength_TH),
     alpha_TH(110),
     zero_TH(1e-7),
     area_TH(area_TH),
     iteration_TH(iteration_TH),
     weight_calculator(Weight_calculator()),
     is_medially_centered(false)
  {
    init();
  }
  /// \endcond

  /**
   * The constructor of a Mean_curvature_skeleton object.
   *
   * @pre the polyhedron is a watertight triangular mesh
   * @param P
   *        triangulated surface mesh used to extract skeleton
   * @param Vertex_index_map
   *        property map for associating an id to each vertex
   * @param Edge_index_map
   *        property map for associating an id to each edge
   * @param omega_H
   *        controls the velocity of movement and approximation quality
   * @param omega_P
   *        controls the smoothness of the medial approximation
   * @param edgelength_TH
   *        edges with length less than `edgelength_TH` will be collapsed
   * @param is_medially_centered
   *        should the skeleton be medially centered?
   * @param area_TH
   *        run_to_converge will stop if the change of area in one iteration
   *        is less than `area_TH`
   * @param iteration_TH
   *        the maximum number of iterations to run
   */
  Mean_curvature_skeleton(HalfedgeGraph& P,
                          VertexIndexMap Vertex_index_map,
                          EdgeIndexMap Edge_index_map,
                          double omega_H,
                          double omega_P,
                          double edgelength_TH,
                          bool is_medially_centered,
                          double area_TH = 1e-5,
                          int iteration_TH = 500
                          )
    :polyhedron(P),
     vertex_id_pmap(Vertex_index_map),
     edge_id_pmap(Edge_index_map),
     omega_H(omega_H),
     omega_P(omega_P),
     edgelength_TH(edgelength_TH),
     alpha_TH(110),
     zero_TH(1e-7),
     area_TH(area_TH),
     iteration_TH(iteration_TH),
     weight_calculator(Weight_calculator()),
     is_medially_centered(is_medially_centered)
  {
    init();
  }

  /// @} Constructor and Destructor


  /// @cond CGAL_DOCUMENT_INTERNAL

  /// \name Setter and Getter
  /// @{

  void set_omega_H(double value)
  {
    omega_H = value;
  }

  void set_omega_P(double value)
  {
    omega_P = value;
  }

  void set_edgelength_TH(double value)
  {
    edgelength_TH = value;
  }

  void set_alpha_TH(double value)
  {
    alpha_TH = value;
  }

  void set_zero_TH(double value)
  {
    zero_TH = value;
  }

  void set_medially_centered(bool value)
  {
    is_medially_centered = value;
  }

  void set_iteration_TH(int value)
  {
    iteration_TH = value;
  }

  HalfedgeGraph& get_halfedge_graph()
  {
    return polyhedron;
  }

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
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end())
      {
        vertex_descriptor vd = *vb;
        fixed_points.push_back(vd->point());
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
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      if (is_vertex_fixed_map.find(id) == is_vertex_fixed_map.end())
      {
          vertex_descriptor vd = *vb;
          non_fixed_points.push_back(vd->point());
      }
    }
  }

  /**
   * Get the correspondent surface points for the skeleton.
   *
   * @param corr
   *        for each skeletal point, record its correspondent surface points
   */
  void get_correspondent_vertices(std::map<vertex_desc, std::vector<int> >& corr)
  {
    corr = skeleton_to_surface;
  }

  /**
   * Get the Voronoi pole for the surface mesh.
   *
   * @param max_poles
   *        for each mesh vertex, record its correspondent Voronoi pole position
   */
  void get_poles(std::vector<Point>& max_poles)
  {
    max_poles.resize(boost::num_vertices(polyhedron));
    vertex_iterator vb, ve;
    int cnt = 0;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int vid = boost::get(vertex_id_pmap, v);
      max_poles[cnt++] = cell_dual[poles[vid]];
    }
  }

  /// @} Setter and Getter

  /// \endcond


  /// \name Public Algorithm API
  /// @{

  /**
   * Extract the skeleton curve for the mesh. The algorithm repeatedly
   * contract the mesh until convergence, then turn the contracted mesh
   * to a curve skeleton.
   *
   * @param g
   *        a boost::graph containing the connectivity of the skeleotn
   * @param points
   *        the locations of the skeletal points
   */
  void extract_skeleton(Graph& g, std::map<vertex_desc, Point>& points)
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
   *        a boost::graph containing the connectivity of the skeleotn
   * @param points
   *        the locations of the skeletal points
   * @param corr
   *        for each skeletal point, record its correspondent surface points
   */
  void extract_skeleton(Graph& g, std::map<vertex_desc, Point>& points,
                        std::map<vertex_desc, std::vector<int> >& corr)
  {
    run_to_converge();
    convert_to_skeleton(g, points);

    corr = skeleton_to_surface;
  }

  /// @cond CGAL_DOCUMENT_INTERNAL

  /**
   * Contract the mesh by mean curvature flow.
   */
  void contract_geometry()
  {
    MCFSKEL_DEBUG(std::cerr << "before contract geometry";)

    update_vertex_id();

    compute_edge_weight();

    int nver = boost::num_vertices(polyhedron);
    int nrows;
    if (is_medially_centered)
    {
      nrows = nver * 3;
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
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor vi = *vb;
      int id = boost::get(vertex_id_pmap, vi);
      int i = new_id[id];
      Point p(X[i], Y[i], Z[i]);
      vi->point() = p;
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

    MCFSKEL_INFO(double area = internal::get_surface_area(polyhedron);)
    MCFSKEL_INFO(std::cout << "area " << area << "\n";)
  }

  /**
   * Run iterations of `contract_geometry()`, `update_topology()` and
   * `detect_degeneracies()` until the change of surface area during one
   * iteration is less than `area_TH` * original surface area.
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

      double area = internal::get_surface_area(polyhedron);
      double area_ratio = fabs(last_area - area) / original_area;

      MCFSKEL_INFO(std::cout << "area " << area << "\n";)
      MCFSKEL_INFO(std::cout << "|area - last_area| / original_area "
                             << area_ratio << "\n";)

      if (area_ratio < area_TH)
      {
        break;
      }
      last_area = area;

      num_iteration++;
      if (num_iteration >= iteration_TH)
      {
        break;
      }
    }
  }

  /**
   * Convert the contracted mesh to a skeleton curve.
   */
  void convert_to_skeleton(Graph& g, std::map<vertex_desc, Point>& points)
  {
    Skeleton skeleton(polyhedron);
    std::map<vertex_desc, std::vector<int> > record;
    skeleton.extract_skeleton(g, points, record);

    skeleton_to_surface.clear();
    typename std::map<vertex_desc, std::vector<int> >::iterator iter;
    for (iter = record.begin(); iter != record.end(); ++iter)
    {
      vertex_desc i = iter->first;
      for (size_t j = 0; j < record[i].size(); ++j)
      {
        int id = record[i][j];
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
    int cnt = 0;
    for (iter = skeleton_to_surface.begin(); iter != skeleton_to_surface.end();
         ++iter)
    {
      vertex_desc i = iter->first;
      cnt += skeleton_to_surface[i].size();
    }

    MCFSKEL_INFO(std::cout << "tracked " << cnt << " vertices\n";)
  }

  /// \endcond

  /// @} Public Algorithm API


private:

  // --------------------------------------------------------------------------
  // Initialization
  // --------------------------------------------------------------------------

  /// Initialize some global data structures such as vertex id.
  void init()
  {
    vertex_iterator vb, ve;

    alpha_TH *= (M_PI / 180.0);
    double area = internal::get_surface_area(polyhedron);
    original_area = area;

    // initialize index maps
    vertex_id_count = 0;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      boost::put(vertex_id_pmap, *vb, vertex_id_count++);
    }
    max_id = vertex_id_count;

    edge_iterator eb, ee;
    int idx = 0;
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
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
    edge_weight.reserve(boost::num_edges(polyhedron));
    edge_iterator eb, ee;
    for(boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      edge_weight.push_back(this->weight_calculator(*eb, polyhedron));
    }
  }

  /// Assemble the left hand side.
  void assemble_LHS(typename SparseLinearAlgebraTraits_d::Matrix& A)
  {
    MCFSKEL_DEBUG(std::cerr << "start LHS\n";)

    int nver = boost::num_vertices(polyhedron);

    Point_inside_polyhedron_3<HalfedgeGraph, Kernel> test_inside(polyhedron);

    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
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

    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
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
      for (boost::tie(e, e_end) = boost::in_edges(*vb, polyhedron); e != e_end; ++e)
      {
        vertex_descriptor vj = boost::source(*e, polyhedron);
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

    Point_inside_polyhedron_3<HalfedgeGraph, Kernel> test_inside(polyhedron);

    // assemble right columns of linear system
    int nver = boost::num_vertices(polyhedron);
    vertex_iterator vb, ve;
    for (int i = 0; i < nver; ++i)
    {
      Bx[i] = 0;
      By[i] = 0;
      Bz[i] = 0;
    }

    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
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
      Bx[i + nver] = vi->point().x() * oh;
      By[i + nver] = vi->point().y() * oh;
      Bz[i + nver] = vi->point().z() * oh;
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
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
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
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      vertex_descriptor vi = boost::source(*eb, polyhedron);
      vertex_descriptor vj = boost::target(*eb, polyhedron);
      size_t vi_idx = boost::get(vertex_id_pmap, vi);
      size_t vj_idx = boost::get(vertex_id_pmap, vj);

      if (is_vertex_fixed_map.find(vi_idx) != is_vertex_fixed_map.end()
       && is_vertex_fixed_map.find(vj_idx) != is_vertex_fixed_map.end())
      {
        fixed_edge_map.set_is_fixed(*eb, true);
      }
    }

    int edge_id = 0;
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, edge_id++);
    }

    // This is a stop predicate (defines when the algorithm terminates).
    // The simplification stops when the length of all edges is greater
    // than the minimum threshold.
    CGAL::internal::Minimum_length_predicate<HalfedgeGraph> stop(edgelength_TH);

    // midpoint placement without geometric test
    SMS::Geometric_test_skipper< SMS::Midpoint_placement<HalfedgeGraph> > placement;

    internal::Track_correspondence_visitor<HalfedgeGraph> vis;
    if (is_medially_centered)
    {
      vis = internal::Track_correspondence_visitor<HalfedgeGraph>(&correspondence,
                                                   &poles, &cell_dual, max_id);
    }
    else
    {
      vis = internal::Track_correspondence_visitor<HalfedgeGraph>(&correspondence, max_id);
    }

    int r = SMS::edge_collapse
                (polyhedron
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
      Point p1 = v1->point();
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
    edges.reserve(boost::num_edges(polyhedron));
    edge_iterator eb, ee;

    boost::tie(eb, ee) = boost::edges(polyhedron);
    std::copy(eb, ee, std::back_inserter(edges));

    int cnt = 0;
    for (size_t i = 0; i < edges.size(); ++i)
    {
      edge_descriptor ed = edges[i];
      if (fixed_edge_map.is_fixed(ed))
      {
        continue;
      }

      vertex_descriptor vi = boost::source(ed, polyhedron);
      vertex_descriptor vj = boost::target(ed, polyhedron);
      double edge_length = sqrt(squared_distance(vi->point(), vj->point()));
      if (internal::is_collapse_ok(polyhedron, ed) && edge_length < edgelength_TH)
      {
        Point p = midpoint(
          boost::get(vertex_point, polyhedron, boost::source(ed, polyhedron)),
          boost::get(vertex_point, polyhedron, boost::target(ed, polyhedron)));

        // invalidate the edges that will be collapsed
        // since the mesh is closed, 6 halfedges will be collapsed
        fixed_edge_map.set_is_fixed(ed, true);
        fixed_edge_map.set_is_fixed(ed->opposite(), true);
        fixed_edge_map.set_is_fixed(ed->prev(), true);
        fixed_edge_map.set_is_fixed(ed->prev()->opposite(), true);
        fixed_edge_map.set_is_fixed(ed->opposite()->prev(), true);
        fixed_edge_map.set_is_fixed(ed->opposite()->prev()->opposite(), true);

        vertex_descriptor v = SMS::halfedge_collapse(ed, polyhedron);
        boost::put(vertex_point, polyhedron, v, p);

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
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      vertex_descriptor vi = boost::source(*eb, polyhedron);
      vertex_descriptor vj = boost::target(*eb, polyhedron);
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
    int ne = boost::num_edges(polyhedron);
    halfedge_angle.resize(ne, 0);

    edge_iterator eb, ee;
    int idx = 0;
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, idx++);
    }

    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      int e_id = boost::get(edge_id_pmap, *eb);
      edge_descriptor ed = *eb;

      if (ed->is_border())
      {
        halfedge_angle[e_id] = -1;
      }
      else
      {
        vertex_descriptor vi = boost::source(ed, polyhedron);
        vertex_descriptor vj = boost::target(ed, polyhedron);
        edge_descriptor ed_next = ed->next();
        vertex_descriptor vk = boost::target(ed_next, polyhedron);
        Point pi = vi->point();
        Point pj = vj->point();
        Point pk = vk->point();

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
    Point ps = vs->point();
    Point pt = vt->point();
    Point pk = vk->point();
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
    int ne = boost::num_edges(polyhedron);
    compute_incident_angle();

    int cnt = 0;
    edge_iterator eb, ee;
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
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

      vertex_descriptor vs = boost::source(ei, polyhedron);
      vertex_descriptor vt = boost::target(ei, polyhedron);

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
      vertex_descriptor vk = boost::target(ek, polyhedron);
      Point pn = project_vertex(vs, vt, vk);
      edge_descriptor en = internal::mesh_split(polyhedron, ei, pn);
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
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int idx = boost::get(vertex_id_pmap, v);

      if (is_vertex_fixed_map.find(idx) == is_vertex_fixed_map.end())
      {
        bool willbefixed = internal::is_vertex_degenerate(polyhedron, v, edgelength_TH);
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
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int idx = boost::get(vertex_id_pmap, v);
      if (is_vertex_fixed_map.find(idx) == is_vertex_fixed_map.end())
      {
        bool willbefixed = false;
        int bad_counter = 0;

        in_edge_iterator eb, ee;
        for (boost::tie(eb, ee) = boost::in_edges(v, polyhedron); eb != ee; ++eb)
        {
          edge_descriptor edge = *eb;
          vertex_descriptor v0 = boost::source(edge, polyhedron);
          vertex_descriptor v1 = boost::target(edge, polyhedron);
          double length = sqrt(squared_distance(v0->point(), v1->point()));
          if (length < elength_fixed)
          {
            if (!internal::is_collapse_ok(polyhedron, edge))
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
    point_to_pole.resize(boost::num_vertices(polyhedron));

    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int vid = boost::get(vertex_id_pmap, v);
      Exact_point tp((v->point()).x(), (v->point()).y(), (v->point()).z());
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
  }

  /// Compute an approximate vertex normal for all vertices.
  void compute_vertex_normal()
  {
    normals.resize(boost::num_vertices(polyhedron));

    vertex_iterator vb, ve;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor v = *vb;
      int vid = boost::get(vertex_id_pmap, v);
      normals[vid] = internal::get_vertex_normal<typename HalfedgeGraph::Vertex,Kernel>(*v);
    }
  }

};

/// \ingroup PkgMeanCurvatureSkeleton3
/// @brief Extract the medially centered curve skeleton for the mesh.
///
/// @pre the polyhedron is a watertight triangular mesh
///
/// @tparam HalfedgeGraph
///         a model of `HalfedgeGraph`
/// @tparam SparseLinearAlgebraTraits_d
///         a model of `SparseLinearAlgebraTraitsWithPreFactor_d`
/// @tparam VertexIndexMap
///         a model of `ReadWritePropertyMap`</a>
///         with Mean_curvature_skeleton::vertex_descriptor as key and
///         `unsigned int` as value type
/// @tparam EdgeIndexMap
///         a model of `ReadWritePropertyMap`</a>
///         with Mean_curvature_skeleton::edge_descriptor as key and
///         `unsigned int` as value type
/// @tparam Graph
///         a model of boost::adjacency_list
///         data structure for skeleton curve
///
/// @param g
///        a boost::graph containing the connectivity of the skeleotn
/// @param points
///        the locations of the skeletal points
/// @param corr
///        for each skeletal point, record its correspondent surface points
/// @param P
///        triangulated surface mesh used to extract skeleton
/// @param Vertex_index_map
///        property map for associating an id to each vertex
/// @param Edge_index_map
///        property map for associating an id to each edge
/// @param omega_H
///        controls the velocity of movement and approximation quality
/// @param omega_P
///        controls the smoothness of the medial approximation
/// @param edgelength_TH
///        edges with length less than `edgelength_TH` will be collapsed
/// @param is_medially_centered
///        should the skeleton be medially centered?
/// @param area_TH
///        run_to_converge will stop if the change of area in one iteration
///        is less than `area_TH`
/// @param iteration_TH
///        the maximum number of iterations to run
template <class HalfedgeGraph,
          class SparseLinearAlgebraTraits_d,
          class VertexIndexMap,
          class EdgeIndexMap,
          class Graph>
void extract_medial_skeleton(Graph& g,
                            std::map<typename boost::graph_traits<Graph>::vertex_descriptor,
                            typename HalfedgeGraph::Traits::Point_3>& points,
                            std::map<typename boost::graph_traits<Graph>::vertex_descriptor,
                            std::vector<int> >& corr,
                            HalfedgeGraph& P,
                            VertexIndexMap Vertex_index_map,
                            EdgeIndexMap Edge_index_map,
                            double omega_H,
                            double omega_P,
                            double edgelength_TH,
                            double area_TH = 1e-5,
                            int iteration_TH = 500)
{
  typedef Mean_curvature_skeleton<HalfedgeGraph, SparseLinearAlgebraTraits_d,
  VertexIndexMap, EdgeIndexMap, Graph> MCFSKEL;

  MCFSKEL mcs(P, Vertex_index_map, Edge_index_map,
  omega_H, omega_P, edgelength_TH, true, area_TH, iteration_TH);

  mcs.run_to_converge();
  mcs.convert_to_skeleton(g, points);

  mcs.get_correspondent_vertices(corr);
}

/// \ingroup PkgMeanCurvatureSkeleton3
/// @brief Extract the medially centered curve skeleton for the mesh.
///
/// @pre the polyhedron is a watertight triangular mesh
///
/// @tparam HalfedgeGraph
///         a model of `HalfedgeGraph`
/// @tparam SparseLinearAlgebraTraits_d
///         a model of `SparseLinearAlgebraTraitsWithPreFactor_d`
/// @tparam VertexIndexMap
///         a model of `ReadWritePropertyMap`</a>
///         with Mean_curvature_skeleton::vertex_descriptor as key and
///         `unsigned int` as value type
/// @tparam EdgeIndexMap
///         a model of `ReadWritePropertyMap`</a>
///         with Mean_curvature_skeleton::edge_descriptor as key and
///         `unsigned int` as value type
/// @tparam Graph
///         a model of boost::adjacency_list
///         data structure for skeleton curve
///
/// @param g
///        a boost::graph containing the connectivity of the skeleotn
/// @param points
///        the locations of the skeletal points
/// @param corr
///        for each skeletal point, record its correspondent surface points
/// @param P
///        triangulated surface mesh used to extract skeleton
/// @param Vertex_index_map
///        property map for associating an id to each vertex
/// @param Edge_index_map
///        property map for associating an id to each edge
/// @param omega_H
///        controls the velocity of movement and approximation quality
/// @param omega_P
///        controls the smoothness of the medial approximation
/// @param edgelength_TH
///        edges with length less than `edgelength_TH` will be collapsed
/// @param is_medially_centered
///        should the skeleton be medially centered?
/// @param area_TH
///        run_to_converge will stop if the change of area in one iteration
///        is less than `area_TH`
/// @param iteration_TH
///        the maximum number of iterations to run
template <class HalfedgeGraph,
          class SparseLinearAlgebraTraits_d,
          class VertexIndexMap,
          class EdgeIndexMap,
          class Graph>
void extract_medial_skeleton(Graph& g,
                            std::map<typename boost::graph_traits<Graph>::vertex_descriptor,
                            typename HalfedgeGraph::Traits::Point_3>& points,
                            HalfedgeGraph& P,
                            VertexIndexMap Vertex_index_map,
                            EdgeIndexMap Edge_index_map,
                            double omega_H,
                            double omega_P,
                            double edgelength_TH,
                            double area_TH = 1e-5,
                            int iteration_TH = 500)
{
  typedef Mean_curvature_skeleton<HalfedgeGraph, SparseLinearAlgebraTraits_d,
  VertexIndexMap, EdgeIndexMap, Graph> MCFSKEL;

  MCFSKEL mcs(P, Vertex_index_map, Edge_index_map,
  omega_H, omega_P, edgelength_TH, true, area_TH, iteration_TH);

  mcs.run_to_converge();
  mcs.convert_to_skeleton(g, points);
}

/// \ingroup PkgMeanCurvatureSkeleton3
/// @brief Extract the medially centered curve skeleton for the mesh.
///
/// @pre the polyhedron is a watertight triangular mesh
///
/// @tparam HalfedgeGraph
///         a model of `HalfedgeGraph`
/// @tparam SparseLinearAlgebraTraits_d
///         a model of `SparseLinearAlgebraTraitsWithPreFactor_d`
/// @tparam VertexIndexMap
///         a model of `ReadWritePropertyMap`</a>
///         with Mean_curvature_skeleton::vertex_descriptor as key and
///         `unsigned int` as value type
/// @tparam EdgeIndexMap
///         a model of `ReadWritePropertyMap`</a>
///         with Mean_curvature_skeleton::edge_descriptor as key and
///         `unsigned int` as value type
/// @tparam Graph
///         a model of boost::adjacency_list
///         data structure for skeleton curve
///
/// @param g
///        a boost::graph containing the connectivity of the skeleotn
/// @param points
///        the locations of the skeletal points
/// @param P
///        triangulated surface mesh used to extract skeleton
/// @param Vertex_index_map
///        property map for associating an id to each vertex
/// @param Edge_index_map
///        property map for associating an id to each edge
/// @param omega_H
///        controls the velocity of movement and approximation quality
/// @param edgelength_TH
///        edges with length less than `edgelength_TH` will be collapsed
/// @param is_medially_centered
///        should the skeleton be medially centered?
/// @param area_TH
///        run_to_converge will stop if the change of area in one iteration
///        is less than `area_TH`
/// @param iteration_TH
///        the maximum number of iterations to run
template <class HalfedgeGraph,
          class SparseLinearAlgebraTraits_d,
          class VertexIndexMap,
          class EdgeIndexMap,
          class Graph>
void extract_skeleton(Graph& g,
                      std::map<typename boost::graph_traits<Graph>::vertex_descriptor,
                      typename HalfedgeGraph::Traits::Point_3>& points,
                      std::map<typename boost::graph_traits<Graph>::vertex_descriptor,
                      std::vector<int> >& corr,
                      HalfedgeGraph& P,
                      VertexIndexMap Vertex_index_map,
                      EdgeIndexMap Edge_index_map,
                      double omega_H,
                      double edgelength_TH,
                      double area_TH = 1e-5,
                      int iteration_TH = 500)
{
  typedef Mean_curvature_skeleton<HalfedgeGraph, SparseLinearAlgebraTraits_d,
  VertexIndexMap, EdgeIndexMap, Graph> MCFSKEL;

  MCFSKEL mcs(P, Vertex_index_map, Edge_index_map,
  omega_H, 0.0, edgelength_TH, false, area_TH, iteration_TH);

  mcs.run_to_converge();
  mcs.convert_to_skeleton(g, points);

  mcs.get_correspondent_vertices(corr);
}

/// \ingroup PkgMeanCurvatureSkeleton3
/// @brief Extract the medially centered curve skeleton for the mesh.
///
/// @pre the polyhedron is a watertight triangular mesh
///
/// @tparam HalfedgeGraph
///         a model of `HalfedgeGraph`
/// @tparam SparseLinearAlgebraTraits_d
///         a model of `SparseLinearAlgebraTraitsWithPreFactor_d`
/// @tparam VertexIndexMap
///         a model of `ReadWritePropertyMap`</a>
///         with Mean_curvature_skeleton::vertex_descriptor as key and
///         `unsigned int` as value type
/// @tparam EdgeIndexMap
///         a model of `ReadWritePropertyMap`</a>
///         with Mean_curvature_skeleton::edge_descriptor as key and
///         `unsigned int` as value type
/// @tparam Graph
///         a model of boost::adjacency_list
///         data structure for skeleton curve
///
/// @param g
///        a boost::graph containing the connectivity of the skeleotn
/// @param points
///        the locations of the skeletal points
/// @param P
///        triangulated surface mesh used to extract skeleton
/// @param Vertex_index_map
///        property map for associating an id to each vertex
/// @param Edge_index_map
///        property map for associating an id to each edge
/// @param omega_H
///        controls the velocity of movement and approximation quality
/// @param edgelength_TH
///        edges with length less than `edgelength_TH` will be collapsed
/// @param is_medially_centered
///        should the skeleton be medially centered?
/// @param area_TH
///        run_to_converge will stop if the change of area in one iteration
///        is less than `area_TH`
/// @param iteration_TH
///        the maximum number of iterations to run
template <class HalfedgeGraph,
          class SparseLinearAlgebraTraits_d,
          class VertexIndexMap,
          class EdgeIndexMap,
          class Graph>
void extract_skeleton(Graph& g,
                      std::map<typename boost::graph_traits<Graph>::vertex_descriptor,
                      typename HalfedgeGraph::Traits::Point_3>& points,
                      HalfedgeGraph& P,
                      VertexIndexMap Vertex_index_map,
                      EdgeIndexMap Edge_index_map,
                      double omega_H,
                      double edgelength_TH,
                      double area_TH = 1e-5,
                      int iteration_TH = 500)
{
  typedef Mean_curvature_skeleton<HalfedgeGraph, SparseLinearAlgebraTraits_d,
  VertexIndexMap, EdgeIndexMap, Graph> MCFSKEL;

  MCFSKEL mcs(P, Vertex_index_map, Edge_index_map,
  omega_H, 0.0, edgelength_TH, false, area_TH, iteration_TH);

  mcs.run_to_converge();
  mcs.convert_to_skeleton(g, points);
}

} //namespace CGAL

#endif // CGAL_MEAN_CURVATURE_SKELETON_H
