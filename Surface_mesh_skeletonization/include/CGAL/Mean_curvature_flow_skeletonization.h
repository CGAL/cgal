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

#ifndef CGAL_MEAN_CURVATURE_FLOW_SKELETONIZATION_H
#define CGAL_MEAN_CURVATURE_FLOW_SKELETONIZATION_H

#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <CGAL/Default.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/FaceGraph_to_Polyhedron_3.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/foreach.hpp>

#include <CGAL/boost/graph/iterator.h>

// Compute cotangent Laplacian
#include <CGAL/internal/Surface_mesh_skeletonization/Weights.h>

// Compute the vertex normal
#include <CGAL/internal/Surface_mesh_skeletonization/get_normal.h>

// Simplification function
#include <CGAL/boost/graph/Euler_operations.h>

// Curve skeleton data structure
#include <CGAL/internal/Surface_mesh_skeletonization/Curve_skeleton.h>

// For Voronoi diagram
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

// For debugging macro
#include <CGAL/internal/Surface_mesh_skeletonization/Debug.h>

// Some helper functions
#include <CGAL/internal/Surface_mesh_skeletonization/Utility.h>

// For detect_degenarcy
#include <CGAL/internal/Surface_mesh_skeletonization/Detect_degeneracy.h>

// Inside mesh test
#include <CGAL/Side_of_triangle_mesh.h>

// Compute bounding box
#include <CGAL/Bbox_3.h>

#include <queue>

// for default parameters
#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>  // for sparse linear system solver
#include <Eigen/SparseCholesky>
// #include <Eigen/CholmodSupport>
#endif

namespace CGAL {

namespace internal{

template < class Refs, class Point, class ID, class vertex_descriptor>
struct Skel_HDS_vertex_type : public HalfedgeDS_vertex_max_base_with_id<Refs, Point, ID>
{
  typedef HalfedgeDS_vertex_max_base_with_id<Refs, Point, ID> Base;
  Skel_HDS_vertex_type() : Base (), is_fixed(false)  {}
  Skel_HDS_vertex_type( Point const& p) : Base(p), is_fixed(false) {}
  std::vector<vertex_descriptor> vertices;
  Point pole;
  bool is_fixed;
};

template <class vertex_descriptor>
struct Skel_polyhedron_items_3: CGAL::Polyhedron_items_with_id_3 {
    template < class Refs, class Traits>
    struct Vertex_wrapper {
        typedef typename Traits::Point_3 Point;
      typedef Skel_HDS_vertex_type< Refs, Point, std::size_t, vertex_descriptor> Vertex;
    };
};

} //end of namespace internal


/// \ingroup PkgMeanCurvatureSkeleton3
/// Function object that enables to extract the mean curvature
/// flow skeleton of a triangulated surface mesh.
///
/// The algorithm used takes as input a triangulated surface mesh and iteratively contracts the surface mesh
/// following the mean curvature flow \cgalCite{tagliasacchi2012mean}. The intermediate contracted surface
/// mesh is called the <em>meso-skeleton</em>.
/// After each iteration, the meso-skeleton is locally remeshed using angle split and edge contraction.
/// The process ends when the modification of the meso-skeleton between two iterations is small.
///
/// @tparam TriangleMesh
///         a model of `FaceListGraph`
///
/// @tparam Traits
///         a model of `MeanCurvatureSkeletonizationTraits`<br>
///         <b>%Default:</b>
/// \code
///     CGAL::Kernel_traits<
///       boost::property_traits<
///          boost::property_map<TriangleMesh, CGAL::vertex_point_t>::type
///        >::value_type
///      >::Kernel
/// \endcode
///
/// @tparam VertexPointMap
///         a model of `ReadWritePropertyMap`
///         with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key and
///         `Traits::Point_3` as value type.<br>
///         <b>%Default:</b>
/// \code
///   boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type.
/// \endcode
///
/// @tparam SolverTraits_
///         a model of `NormalEquationSparseLinearAlgebraTraits_d`.<br>
///         <b>%Default:</b> If \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available
///         and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits` is provided as default parameter:
/// \code
///      CGAL::Eigen_solver_traits<
///        Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> >
///      >
/// \endcode
///
/// @cond CGAL_DOCUMENT_INTERNAL
/// @tparam Degeneracy_algorithm_tag
///         tag for selecting the degeneracy detection algorithm
/// @endcond
template <class TriangleMesh,
          class Traits_ = Default,
          class VertexPointMap_ = Default,
          class SolverTraits_ = Default>
class Mean_curvature_flow_skeletonization
{
// Public types
public:

/// \name Types
/// @{
  // Template parameters
  #ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    VertexPointMap_,
    typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type
  >::type VertexPointMap;

  typedef typename Default::Get<
    Traits_,
    typename Kernel_traits<typename boost::property_traits<typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::type>::value_type>::Kernel
  >::type Traits;
  #endif

  #ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    SolverTraits_,
  #if defined(CGAL_EIGEN3_ENABLED)
      CGAL::Eigen_solver_traits<
        Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> >
//        Eigen::CholmodDecomposition< Eigen::SparseMatrix<double> >
      >
  #else
    SolverTraits_ // no parameter provided, and Eigen is not enabled: so don't compile!
  #endif
  >::type SolverTraits;
  #endif

  /// @cond CGAL_DOCUMENT_INTERNAL
  typedef typename Traits::Point_3                                             Point;
  typedef typename Traits::Vector_3                                             Vector;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor  Input_vertex_descriptor;
  typedef CGAL::Polyhedron_3<Traits,internal::Skel_polyhedron_items_3<Input_vertex_descriptor> > mTriangleMesh;
  typedef typename boost::property_map<mTriangleMesh, CGAL::vertex_point_t>::type mVertexPointMap;
  typedef typename boost::property_map<mTriangleMesh, boost::vertex_index_t>::type VertexIndexMap;
  typedef typename boost::property_map<mTriangleMesh, boost::halfedge_index_t>::type HalfedgeIndexMap;


  struct Vmap {
    Point point;
    std::vector<Input_vertex_descriptor> vertices;
  };
  ///@endcond

  /// The graph type representing the skeleton. The vertex property
  /// `Vmap` is a struct with a member `point` of type `Traits::Point_3`
  /// and a member `vertices` of type
  /// `std::vector<boost::graph_traits<TriangleMesh>::%vertex_descriptor>`.
  /// See  <a href="http://www.boost.org/doc/libs/release/libs/graph/doc/adjacency_list.html"><tt>the boost documentation</tt></a> page for more details
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Vmap> Skeleton;


/// @}

  // Repeat mTriangleMesh types
  typedef typename boost::graph_traits<mTriangleMesh>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<mTriangleMesh>::halfedge_descriptor     halfedge_descriptor;
  typedef typename boost::graph_traits<mTriangleMesh>::vertex_iterator         vertex_iterator;
  typedef typename boost::graph_traits<mTriangleMesh>::edge_descriptor         edge_descriptor;
  typedef typename boost::graph_traits<mTriangleMesh>::edge_iterator           edge_iterator;

  // Cotangent weight calculator
  typedef typename internal::Cotangent_weight<mTriangleMesh,
  internal::Cotangent_value_minimum_zero<mTriangleMesh,
  internal::Cotangent_value_Meyer_secure<mTriangleMesh> > >                    Weight_calculator;

  typedef internal::Curve_skeleton<mTriangleMesh,
                                   VertexIndexMap,
                                   HalfedgeIndexMap,
                                   mVertexPointMap>                            Curve_skeleton;

  // Repeat Triangulation types
  typedef CGAL::Exact_predicates_exact_constructions_kernel                    Exact_kernel;
  typedef CGAL::Triangulation_vertex_base_with_info_3
                                            <vertex_descriptor, Exact_kernel>  Vb;
  typedef CGAL::Triangulation_data_structure_3<Vb>                             Tds;
  typedef CGAL::Delaunay_triangulation_3<Exact_kernel, Tds>                    Delaunay;
  typedef typename Delaunay::Point                                             Exact_point;
  typedef typename Delaunay::Cell_handle                                       Cell_handle;
  typedef typename Delaunay::Vertex_handle                                     TriVertex_handle;
  typedef typename Delaunay::Finite_cells_iterator                             Finite_cells_iterator;

// Data members
private:

  /** The meso-skeleton */
  mTriangleMesh m_tmesh;

  /** Storing indices of all vertices. */
  VertexIndexMap m_vertex_id_pmap;
  /** Storing indices of all edges. */
  HalfedgeIndexMap m_hedge_id_pmap;
  /** Storing the point for mTriangleMesh vertex_descriptor. */
  mVertexPointMap m_tmesh_point_pmap;
  /** Traits class. */
  Traits m_traits;

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
  std::size_t m_max_iterations;
  /** Should the skeleton be medially centered? */
  bool m_is_medially_centered;
  /** Are poles computed? */
  bool m_are_poles_computed;

  /** Cotangent weight calculator. */
  Weight_calculator m_weight_calculator;
  /** Storing the weights for edges. */
  std::vector<double> m_edge_weight;
  /** The sparse solver. */
  SolverTraits m_solver;

  /** Assign a unique id to a new vertex. */
  int m_vertex_id_count;
  /** The maximum id for original surface. vertices with ids
   *  greater than `m_max_id` are created during split,
   *  thus will not be considered in correspondence tracking. */
  int m_max_id;
  /** Used when assembling the matrix. */
  std::map<int, int> m_new_id;

  /** The incident angle for a halfedge. */
  std::vector<double> m_halfedge_angle;

  /** The normal of surface points. */
  std::vector<Vector> m_normals;


// Private functions and classes

struct Vertex_to_point{
  mVertexPointMap ppmap;
  Vertex_to_point(mVertexPointMap ppmap): ppmap(ppmap){}
  typedef typename boost::property_traits<mVertexPointMap>::reference result_type;
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
  return std::sqrt(diag);
}


double init_min_edge_length()
{
  vertex_iterator vb, ve;
  boost::tie(vb, ve) = vertices(m_tmesh);
  Vertex_to_point v_to_p(m_tmesh_point_pmap);
  Bbox_3 bbox = CGAL::bbox_3(boost::make_transform_iterator(vb, v_to_p),
                             boost::make_transform_iterator(ve, v_to_p));
  return 0.002 * diagonal_length(bbox);
}

std::size_t collapse_short_edges();
// Public methods
public:

  /// \name Constructor
  ///@{
  #ifdef DOXYGEN_RUNNING
  /**
   * The constructor of a skeletonization object.
   *
   * The algorithm parameters are initialized such that:
   * - `max_triangle_angle() == 110`
   * - `quality_speed_tradeoff() == 0.1`
   * - `medially_centered_speed_tradeoff() == 0.2`
   * - `area_variation_factor() == 0.0001`
   * - `max_iterations() == 500`
   * - `is_medially_centered() == true`
   * - `min_edge_length()` == 0.002 * the length of the diagonal of the bounding box of `tmesh`
   *
   * @pre `tmesh` is a triangulated surface mesh without borders and has exactly one connected component.
   * @param tmesh
   *        input triangulated surface mesh.
   * @param vertex_point_map
   *        property map which associates a point to each vertex of the graph.
   * @param traits
   *        an instance of the traits class.
   */
  Mean_curvature_flow_skeletonization(const TriangleMesh& tmesh,
                                      VertexPointMap vertex_point_map = get(CGAL::vertex_point, tmesh),
                                      Traits traits = Traits());

  #else

  Mean_curvature_flow_skeletonization(const TriangleMesh& tmesh,
                                      VertexPointMap vertex_point_map,
                                      const Traits& traits = Traits())
    : m_traits(traits)
  {
    init(tmesh, vertex_point_map);
  }

  Mean_curvature_flow_skeletonization(const TriangleMesh& tmesh,
                                      const Traits& traits = Traits())
    : m_traits(traits)
  {
    init(tmesh, get(vertex_point, tmesh));
  }
  #endif
  /// @} Constructor

  /// \name Local Remeshing Parameters
  /// @{

  /// During the local remeshing step, a triangle will be split
  /// if it has an angle larger than `max_triangle_angle()`.
  double max_triangle_angle()
  {
    return m_alpha_TH;
  }

  /// During the local remeshing step, an edge will be collapse
  /// if it is length is less than `min_edge_length()`.
  double min_edge_length()
  {
    return m_min_edge_length;
  }

 void set_max_triangle_angle(double value)
  {
    m_alpha_TH = value;
  }

  void set_min_edge_length(double value)
  {
    m_min_edge_length = value;
  }
  /// @}

  /// \name Algorithm Termination Parameters
  /// @{

  /// Maximum number of iterations performed by `contract_until_convergence()`.
  std::size_t max_iterations()
  {
    return m_max_iterations;
  }

  /// The convergence is considered to be reached if the variation of the area of
  /// the meso-skeleton after one iteration is smaller than
  /// `area_variation_factor()*original_area` where `original_area` is the area of the input
  /// triangle mesh.
  double area_variation_factor()
  {
    return m_delta_area;
  }

  void set_max_iterations(std::size_t value)
  {
    m_max_iterations = value;
  }

  void set_area_variation_factor(double value)
  {
    m_delta_area = value;
  }
  /// @}

  /// \name Vertex Motion Parameters
  /// @{

  /// \cgalAdvancedBegin
  /// Controls the velocity of movement and approximation quality:
  /// decreasing this value makes the mean curvature flow based contraction converge
  /// faster, but results in a skeleton of lower quality.
  /// This parameter corresponds to \f$ w_H \f$ in the original publication.
  /// \cgalAdvancedEnd
  double quality_speed_tradeoff()
  {
    return m_omega_H;
  }

  /// If `true`, the meso-skeleton placement will be attracted by an approximation
  /// of the medial axis of the mesh during the contraction steps, so will be the result skeleton.
  // (an additional energy is used during the contraction using the Voronoi poles of the input triangulated mesh
  // as attractors).
  bool is_medially_centered()
  {
    return m_is_medially_centered;
  }

  /// \cgalAdvancedBegin
  /// Controls the smoothness of the medial approximation:
  /// increasing this value results in a (less smooth) skeleton closer
  /// to the medial axis, as well as a lower convergence speed.
  /// It is only used if `is_medially_centered()==true`.
  /// This parameter corresponds to \f$ w_M \f$ in the original publication.
  /// \cgalAdvancedEnd
  double medially_centered_speed_tradeoff()
  {
    return m_omega_P;
  }

  void set_quality_speed_tradeoff(double value)
  {
    m_omega_H = value;
  }

  void set_is_medially_centered(bool value)
  {
    m_is_medially_centered = value;
  }

  void set_medially_centered_speed_tradeoff(double value)
  {
    m_omega_P = value;
  }

  /// \cond SKIP_FROM_MANUAL
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
    BOOST_FOREACH(vertex_descriptor vd, vertices(m_tmesh))
    {
      if (vd->is_fixed)
        fixed_points.push_back(get(m_tmesh_point_pmap, vd));
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
    BOOST_FOREACH(vertex_descriptor vd, vertices(m_tmesh))
    {
      if (!vd->is_fixed)
        non_fixed_points.push_back(get(m_tmesh_point_pmap, vd));
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
    max_poles.resize(num_vertices(m_tmesh));
    int cnt = 0;
    BOOST_FOREACH(vertex_descriptor v, vertices(m_tmesh))
    {
      max_poles[cnt++] = v->pole;
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

   * @param skeleton
   *        graph that will contain the skeleton of the input triangulated surface mesh.
   *        For each vertex descriptor `vd` of `skeleton`, the corresponding point
   *        and the set of input vertices that contracted to `vd` can be retrieved
   *        using `skeleton[vd].point` and `skeleton[vd].vertices` respectively.
   */
  void operator()(Skeleton& skeleton)
  {
    contract_until_convergence();
    convert_to_skeleton(skeleton);
  }
  /// @}

  /// \name Low Level Functions
  /// \cgalAdvancedBegin
  /// The following functions enable the user to run the mean curvature flow skeletonization algorithm step by step.
  /// \cgalAdvancedEnd
  /// @{

  /**
   * Runs one contraction step following the mean curvature flow.
   */
  void contract_geometry()
  {
    MCFSKEL_DEBUG(std::cerr << "before contract geometry";)

    update_vertex_id();

    compute_edge_weight();

  // AF: attention: num_vertices will not decrease for a Surface_mesh
    int nver = static_cast<int>(num_vertices(m_tmesh));
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
    typename SolverTraits::Matrix A(nrows, nver);
    assemble_LHS(A);

    typename SolverTraits::Vector X(nver), Bx(nrows);
    typename SolverTraits::Vector Y(nver), By(nrows);
    typename SolverTraits::Vector Z(nver), Bz(nrows);
    assemble_RHS(Bx, By, Bz);

    MCFSKEL_DEBUG(std::cerr << "before solve\n";)

    // solve "At * A * X = At * B".
    m_solver.normal_equation_factor(A);
    m_solver.normal_equation_solver(Bx, X);
    m_solver.normal_equation_solver(By, Y);
    m_solver.normal_equation_solver(Bz, Z);

    MCFSKEL_DEBUG(std::cerr << "after solve\n";)

    // copy to surface mesh
    BOOST_FOREACH(vertex_descriptor vd, vertices(m_tmesh))
    {
      int id = static_cast<int>(get(m_vertex_id_pmap, vd));
      int i = m_new_id[id];
      Point p =  m_traits.construct_point_3_object()(X[i], Y[i], Z[i]);
      put(m_tmesh_point_pmap, vd, p);
    }

    MCFSKEL_DEBUG(std::cerr << "leave contract geometry\n";)
  }

  /**
   * Collapses edges of the meso-skeleton with length less than `min_edge_length()` and returns the number of edges collapsed.
   */
  std::size_t collapse_edges()
  {
    return collapse_short_edges();
  }

  /**
   * Splits faces of the meso-skeleton having one angle greater than `max_triangle_angle()` and returns the number of faces split.
   */
  std::size_t split_faces()
  {
    MCFSKEL_DEBUG(std::cerr << "before split\n";)

    std::size_t num_splits = 0;
    while (true)
    {
      if (num_vertices(m_tmesh) <= 3)
      {
        break;
      }
      std::size_t cnt = split_flat_triangles();
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
  std::size_t remesh()
  {
    MCFSKEL_DEBUG(std::cerr << "before collapse edges\n";)

    std::size_t num_collapses = collapse_edges();
    MCFSKEL_INFO(std::cerr << "collapse " << num_collapses << " edges.\n";)

    std::size_t num_splits = split_faces();
    MCFSKEL_INFO(std::cerr << "split " << num_splits << " edges.\n";)

    return num_collapses + num_splits;
  }
#endif

  /**
   * Prevents degenerate vertices to move during the following contraction steps and returns the number of newly fixed vertices.
   */
  std::size_t detect_degeneracies()
  {
    return detect_degeneracies_in_disk();
  }

  /// Performs subsequent calls to `contract_geometry()`, `collapse_edges()`, `split_faces()` and `detect_degeneracies()`
  void contract()
  {
    contract_geometry();
    remesh();
    detect_degeneracies();

    MCFSKEL_DEBUG(print_edges();)

    MCFSKEL_INFO(double area = internal::get_surface_area(m_tmesh, m_tmesh_point_pmap);)
    MCFSKEL_INFO(std::cout << "area " << area << "\n";)
  }


  /**
   * Iteratively calls `contract()`
   * until the change of surface area of the meso-skeleton after one iteration is smaller than
   * `area_variation_factor()*original_area` where `original_area` is the area of the input triangle mesh,
   * or if the maximum number of iterations has been reached.
   */
  void contract_until_convergence()
  {
    double last_area = 0;
    std::size_t num_iteration = 0;
    while (true)
    {
      MCFSKEL_INFO(std::cout << "iteration " << num_iteration + 1 << "\n";)

      contract_geometry();
      remesh();
      detect_degeneracies();

      double area = internal::get_surface_area(m_tmesh, m_tmesh_point_pmap, m_traits);
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
   * @tparam Skeleton
   *         an instantiation of <A href="http://www.boost.org/libs/graph/doc/adjacency_list.html">`boost::adjacency_list`</a>
   *         as a data structure for the skeleton curve.
   * @param skeleton
   *        graph that will contain the skeleton of `tmesh`. It should be empty before passed to the function.
   */
  void convert_to_skeleton(Skeleton& skeleton)
  {
    skeleton.clear();
    Curve_skeleton skeletonization(m_tmesh, m_vertex_id_pmap, m_hedge_id_pmap, m_tmesh_point_pmap);
    skeletonization.extract_skeleton(skeleton);
  }

  /// @} Public Algorithm API


  /// \name Access to the Meso-Skeleton
  /// @{

  /// When using the low level API it is possible to access the intermediate
  /// results of the skeletonization process, called meso-skeleton.
  /// It is a triangulated surface mesh which is model of `FaceListGraph`.
#ifdef DOXYGEN_RUNNING
  typedef unspecified_type Meso_skeleton;
#else
  typedef mTriangleMesh Meso_skeleton;
#endif
  /// Reference to the collapsed triangulated surface mesh.
  Meso_skeleton& meso_skeleton()
  {
    return m_tmesh;
  }

  /// @}

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
    m_alpha_TH = 110 * (CGAL_PI / 180.0);
    m_zero_TH = 1e-7;
  }

  /// Initialize some global data structures such as vertex id.
  void init(const TriangleMesh& tmesh,
            VertexPointMap vpm)
  {
    // copy the input FaceGraph into a mTriangleMesh
    CGAL::FaceGraph_to_Polyhedron_3<TriangleMesh,
                                    VertexPointMap,
                                    typename mTriangleMesh::HalfedgeDS,
                                    false> modifier(tmesh, vpm);

    m_tmesh.delegate(modifier);

    // copy input vertices to keep correspondance
    typename boost::graph_traits<mTriangleMesh>::vertex_iterator vit=vertices(m_tmesh).first;
    BOOST_FOREACH(Input_vertex_descriptor vd, vertices(tmesh) )
      (*vit++)->vertices.push_back(vd);

    //init indices
    typedef typename boost::graph_traits<mTriangleMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<mTriangleMesh>::halfedge_descriptor halfedge_descriptor;
    std::size_t i=0;
    BOOST_FOREACH( vertex_descriptor vd, vertices(m_tmesh) )
      vd->id()=i++;
    i=0;
    BOOST_FOREACH( halfedge_descriptor hd, halfedges(m_tmesh) )
      hd->id()=i++;
    m_hedge_id_pmap = get(boost::halfedge_index, m_tmesh);
    m_vertex_id_pmap = get(boost::vertex_index, m_tmesh);
      //, m_hedge_id_pmap(get(boost::halfedge_index, m_tmesh))
    m_are_poles_computed = false;

    m_original_area = internal::get_surface_area(m_tmesh, m_tmesh_point_pmap, m_traits);

    m_vertex_id_count = static_cast<int>(num_vertices(m_tmesh));
    m_max_id = m_vertex_id_count;

    if (m_is_medially_centered)
      compute_voronoi_pole();

    init_args();
  }

  // --------------------------------------------------------------------------
  // Utilities
  // --------------------------------------------------------------------------
  double get_x(const Vector& v){ return m_traits.compute_x_3_object()(v); }
  double get_y(const Vector& v){ return m_traits.compute_y_3_object()(v); }
  double get_z(const Vector& v){ return m_traits.compute_z_3_object()(v); }
  double get_x(const Point& v){ return m_traits.compute_x_3_object()(v); }
  double get_y(const Point& v){ return m_traits.compute_y_3_object()(v); }
  double get_z(const Point& v){ return m_traits.compute_z_3_object()(v); }
  // --------------------------------------------------------------------------
  // Contraction
  // --------------------------------------------------------------------------

  /// Compute cotangent weights of all edges.
  void compute_edge_weight()
  {
    m_edge_weight.clear();
    m_edge_weight.reserve(2 * num_edges(m_tmesh));
    BOOST_FOREACH(halfedge_descriptor hd, halfedges(m_tmesh))
    {
      m_edge_weight.push_back(m_weight_calculator(hd, m_tmesh, m_traits));
    }
  }

  /// Assemble the left hand side.
  void assemble_LHS(typename SolverTraits::Matrix& A)
  {
    MCFSKEL_DEBUG(std::cerr << "start LHS\n";)

    std::size_t nver = num_vertices(m_tmesh);

    Side_of_triangle_mesh<mTriangleMesh, Traits> test_inside(m_tmesh);

    BOOST_FOREACH(vertex_descriptor vd, vertices(m_tmesh))
    {
      int id = static_cast<int>(get(m_vertex_id_pmap, vd));

      int i = m_new_id[id];
      // if the vertex is fixed
      if (vd->is_fixed)
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
            if (test_inside(vd->pole) == CGAL::ON_BOUNDED_SIDE)
            {
              A.set_coef(i + nver * 2, i, m_omega_P, true);
            }
          }
        }
      }
    }

    BOOST_FOREACH(vertex_descriptor vd, vertices(m_tmesh))
    {
      int id = static_cast<int>(get(m_vertex_id_pmap, vd));
      int i = m_new_id[id];
      double L = 1.0;
      // if the vertex is fixed
      if (vd->is_fixed)
      {
        L = 0;
      }
      double diagonal = 0;
      BOOST_FOREACH(edge_descriptor ed, in_edges(vd, m_tmesh))
      {
        vertex_descriptor vj = source(ed, m_tmesh);
        double wij = m_edge_weight[get(m_hedge_id_pmap, halfedge(ed, m_tmesh))] * 2.0;
        int jd = static_cast<int>(get(m_vertex_id_pmap, vj));
        int j = m_new_id[jd];
        A.set_coef(i, j, wij * L, true);
        diagonal += -wij;
      }
      A.set_coef(i, i, diagonal, true);
    }

    MCFSKEL_DEBUG(std::cerr << "end LHS\n";)
  }

  /// Assemble the right hand side.
  void assemble_RHS(typename SolverTraits::Vector& Bx,
                    typename SolverTraits::Vector& By,
                    typename SolverTraits::Vector& Bz)
  {
    MCFSKEL_DEBUG(std::cerr << "start RHS\n";)

    Side_of_triangle_mesh<mTriangleMesh, Traits> test_inside(m_tmesh);

    // assemble right columns of linear system
    int nver = static_cast<int>(num_vertices(m_tmesh));
    for (int i = 0; i < nver; ++i)
    {
      Bx[i] = 0;
      By[i] = 0;
      Bz[i] = 0;
    }

    BOOST_FOREACH(vertex_descriptor vd, vertices(m_tmesh))
    {
      int id = static_cast<int>(get(m_vertex_id_pmap, vd));
      int i = m_new_id[id];

      double oh, op = 0.0;
      if (vd->is_fixed)
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
            if (test_inside(vd->pole) == CGAL::ON_BOUNDED_SIDE)
            {
              op = m_omega_P;
            }
          }
        }
      }
      Bx[i + nver] = get_x(get(m_tmesh_point_pmap, vd)) * oh;
      By[i + nver] = get_y(get(m_tmesh_point_pmap, vd)) * oh;
      Bz[i + nver] = get_z(get(m_tmesh_point_pmap, vd)) * oh;
      if (m_is_medially_centered)
      {
        double x = get_x(vd->pole);
        double y = get_y(vd->pole);
        double z = get_z(vd->pole);
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

    BOOST_FOREACH(vertex_descriptor vd, vertices(m_tmesh))
    {
      int id = static_cast<int>(get(m_vertex_id_pmap, vd));
      m_new_id[id] = cnt++;
    }
  }

  // --------------------------------------------------------------------------
  // Edge collapse
  // --------------------------------------------------------------------------

  /// Track correspondent original surface points during collapse.
  void update_pole(vertex_descriptor v0, vertex_descriptor vkept)
  {
    if (m_is_medially_centered)
    {
      const Point& pole0 = v0->pole;
      const Point& pole1 = vkept->pole;

      Point p1 = get(m_tmesh_point_pmap, vkept);
      double dis_to_pole0 = m_traits.compute_squared_distance_3_object()(pole0, p1);
      double dis_to_pole1 = m_traits.compute_squared_distance_3_object()(pole1, p1);
      if (dis_to_pole0 < dis_to_pole1)
        vkept->pole = v0->pole;
    }
  }

  bool edge_should_be_collapsed(edge_descriptor ed)
  {
    halfedge_descriptor h = halfedge(ed, m_tmesh);

    vertex_descriptor vi = source(h, m_tmesh);
    vertex_descriptor vj = target(h, m_tmesh);

    // an edge cannot be collapsed if both vertices are degenerate.
    if (vi->is_fixed && vj->is_fixed) return false;

    double sq_edge_length = m_traits.compute_squared_distance_3_object()(
                              get(m_tmesh_point_pmap, vi),
                              get(m_tmesh_point_pmap, vj));
    return sq_edge_length < m_min_edge_length * m_min_edge_length;
  }

  // --------------------------------------------------------------------------
  // Triangle split
  // --------------------------------------------------------------------------

  /// Compute the incident angles for all the halfedges.
  void compute_incident_angle()
  {
    m_halfedge_angle.clear();
    int ne = 2 * static_cast<int>(num_edges(m_tmesh));
    m_halfedge_angle.resize(ne, 0);

    int idx = 0;
    BOOST_FOREACH(halfedge_descriptor hd, halfedges(m_tmesh))
    {
      put(m_hedge_id_pmap, hd, idx++);
    }

    BOOST_FOREACH(halfedge_descriptor hd, halfedges(m_tmesh))
    {
      int e_id = static_cast<int>(get(m_hedge_id_pmap, hd));

      if (is_border(hd, m_tmesh))
      {
        m_halfedge_angle[e_id] = -1;
      }
      else
      {
        vertex_descriptor vi = source(hd, m_tmesh);
        vertex_descriptor vj = target(hd, m_tmesh);
        halfedge_descriptor hd_next = next(hd, m_tmesh);
        vertex_descriptor vk = target(hd_next, m_tmesh);
        Point pi = get(m_tmesh_point_pmap, vi);
        Point pj = get(m_tmesh_point_pmap, vj);
        Point pk = get(m_tmesh_point_pmap, vk);

        double dis2_ij = m_traits.compute_squared_distance_3_object()(pi, pj);
        double dis2_ik = m_traits.compute_squared_distance_3_object()(pi, pk);
        double dis2_jk = m_traits.compute_squared_distance_3_object()(pj, pk);
        double dis_ij = std::sqrt(dis2_ij);
        double dis_ik = std::sqrt(dis2_ik);
        double dis_jk = std::sqrt(dis2_jk);

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

  void normalize(Vector& v)
  {
    double norm = std::sqrt(m_traits.compute_squared_length_3_object()(v));
    v = m_traits.construct_divided_vector_3_object()(v, norm);
  }

  /// Project the vertex `vk` to the line of `vs` and `vt`.
  Point project_vertex(const vertex_descriptor vs,
                       const vertex_descriptor vt,
                       const vertex_descriptor vk,
                       const vertex_descriptor vnew)
  {
    Point ps = get(m_tmesh_point_pmap, vs);
    Point pt = get(m_tmesh_point_pmap, vt);
    Point pk = get(m_tmesh_point_pmap, vk);
    Vector vec_st = m_traits.construct_vector_3_object()(ps, pt);
    Vector vec_sk = m_traits.construct_vector_3_object()(ps, pk);

    normalize(vec_st);
    double t = m_traits.compute_scalar_product_3_object()(vec_st,vec_sk);
    Point st = m_traits.construct_point_3_object()( get_x(vec_st) * t, get_y(vec_st) * t, get_z(vec_st) * t);
    Point pn = m_traits.construct_point_3_object()( get_x(ps) + get_x(st), get_y(ps) + get_y(st), get_z(ps) + get_z(st));

    // project the pole
    if (m_is_medially_centered)
    {
      const Point& pole_s = vs->pole;
      const Point& pole_t = vt->pole;
      Vector pole_st = m_traits.construct_vector_3_object()(pole_s, pole_t );
      normalize(pole_st);
      vnew->pole =  m_traits.construct_translated_point_3_object()(
                        pole_s,
                        m_traits.construct_scaled_vector_3_object()(pole_st, t)
                     );
    }
    return pn;
  }

  /// Split triangles with an angle greater than `alpha_TH`.
  std::size_t split_flat_triangles()
  {
    compute_incident_angle();

    // collect edges that should be split because
    // both opposite angle are larger than the
    // threshold
    std::vector<edge_descriptor> edges_to_split;
    BOOST_FOREACH(edge_descriptor ed, edges(m_tmesh))
    {
      halfedge_descriptor ei = halfedge(ed, m_tmesh);
      halfedge_descriptor ej = opposite(ei, m_tmesh);
      int ei_id = static_cast<int>(get(m_hedge_id_pmap, ei));
      int ej_id = static_cast<int>(get(m_hedge_id_pmap, ej));

      double angle_i = m_halfedge_angle[ei_id];
      double angle_j = m_halfedge_angle[ej_id];
      if (angle_i >= m_alpha_TH && angle_j >= m_alpha_TH)
        edges_to_split.push_back( ed );
    }

    // now split the edge
    std::size_t cnt = 0;
    BOOST_FOREACH(edge_descriptor ed, edges_to_split)
    {
      halfedge_descriptor ei = halfedge(ed, m_tmesh);
      halfedge_descriptor ej = opposite(ei, m_tmesh);
      int ei_id = static_cast<int>(get(m_hedge_id_pmap, ei));
      int ej_id = static_cast<int>(get(m_hedge_id_pmap, ej));

      vertex_descriptor vs = source(ei, m_tmesh);
      vertex_descriptor vt = target(ei, m_tmesh);

      double angle_i = m_halfedge_angle[ei_id];
      double angle_j = m_halfedge_angle[ej_id];

      halfedge_descriptor ek = next(angle_i > angle_j ? ei : ej, m_tmesh);
      vertex_descriptor vk = target(ek, m_tmesh);

      // split the edge
      halfedge_descriptor en = m_tmesh.split_edge(ei);
      // split the incident faces
      Euler::split_face(en, next(ei,m_tmesh), m_tmesh);
      if (! is_border(ej,m_tmesh))
      {
        Euler::split_face(ej, next(next(ej,m_tmesh),m_tmesh), m_tmesh);
      }

      // set id for new vertex
      put(m_vertex_id_pmap, target(en,m_tmesh), m_vertex_id_count++);
      Point pn = project_vertex(vs, vt, vk, target(en, m_tmesh));
      // set point of new vertex
      put(m_tmesh_point_pmap, target(en,m_tmesh), pn);
      target(en,m_tmesh)->vertices.clear(); // do no copy the info
      ++cnt;
    }
    return cnt;
  }

  // --------------------------------------------------------------------------
  // Degeneracy detection
  // --------------------------------------------------------------------------

  /// Test degeneracy of a vertex by counting the euler characteristic of
  /// its local neighborhood disk.
  std::size_t detect_degeneracies_in_disk()
  {
    std::size_t num_fixed = 0;
    BOOST_FOREACH(vertex_descriptor v, vertices(m_tmesh))
    {
      if (!v->is_fixed)
      {
        bool willbefixed = internal::is_vertex_degenerate(m_tmesh, m_tmesh_point_pmap,
                                                          v, m_min_edge_length, m_traits);
        if (willbefixed)
        {
          v->is_fixed=true;
          ++num_fixed;
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

    std::vector<std::pair<Exact_point, vertex_descriptor> > points;
    std::vector<std::vector<int> > point_to_pole(num_vertices(m_tmesh));

    BOOST_FOREACH(vertex_descriptor v, vertices(m_tmesh))
    {
      const Point& input_pt = get(m_tmesh_point_pmap, v);
      Exact_point tp(get_x(input_pt), get_y(input_pt), get_z(input_pt));
      points.push_back(std::make_pair(tp, v));
    }

    Delaunay T(points.begin(), points.end());

    Finite_cells_iterator cit;
    int cell_id = 0;
    std::vector<Point> cell_dual;
    cell_dual.reserve(T.number_of_cells());
    for (cit = T.finite_cells_begin(); cit != T.finite_cells_end(); ++cit)
    {
      Cell_handle cell = cit;
      Exact_point point = T.dual(cell);
      cell_dual.push_back(
        m_traits.construct_point_3_object()(
          to_double(point.x()),
          to_double(point.y()),
          to_double(point.z())
        )
      );
      // each cell has 4 incident vertices
      for (int i = 0; i < 4; ++i)
      {
        TriVertex_handle vt = cell->vertex(i);
        std::size_t id = get(m_vertex_id_pmap, vt->info());
        point_to_pole[id].push_back(cell_id);
      }
      ++cell_id;
    }

    typedef std::pair<Exact_point, vertex_descriptor> Pair_type;
    BOOST_FOREACH(const Pair_type& p, points)
    {
      std::size_t vid = get(m_vertex_id_pmap, p.second);
      Point surface_point = get(m_tmesh_point_pmap, p.second);

      double max_neg_t = 1;
      int max_neg_i = -1;

      for (size_t j = 0; j < point_to_pole[vid].size(); ++j)
      {
        int pole_id = point_to_pole[vid][j];
        Point cell_point = cell_dual[pole_id];
        Vector vt = m_traits.construct_vector_3_object()(surface_point, cell_point);
        Vector n = m_normals[vid];

        double t = m_traits.compute_scalar_product_3_object()(vt, n);

        // choose the one with maximum distance along the normal
        if (t < 0 && t < max_neg_t)
        {
          max_neg_i = pole_id;
          max_neg_t = t;
        }
      }

      p.second->pole = cell_dual[max_neg_i];
    }
    m_are_poles_computed = true;
  }

  /// Compute an approximate vertex normal for all vertices.
  void compute_vertex_normal()
  {
    m_normals.resize(num_vertices(m_tmesh));

    BOOST_FOREACH(vertex_descriptor v, vertices(m_tmesh))
    {
      int vid = static_cast<int>(get(m_vertex_id_pmap, v));
      m_normals[vid] = internal::get_vertex_normal(*v, m_traits);
    }
  }

  // --------------------------------------------------------------------------
  // Debug
  // --------------------------------------------------------------------------

  void print_edges()
  {
    std::map<halfedge_descriptor, bool> visited;

    BOOST_FOREACH(halfedge_descriptor hd, halfedges(m_tmesh))
    {
      if (!visited[hd])
      {
        vertex_descriptor vi = source(hd, m_tmesh);
        vertex_descriptor vj = target(hd, m_tmesh);
        size_t vi_idx = get(m_vertex_id_pmap, vi);
        size_t vj_idx = get(m_vertex_id_pmap, vj);
        std::cout << vi_idx << " " << vj_idx << "\n";

        visited[hd] = true;
        visited[opposite(hd,m_tmesh)] = true;
      }
    }
  }
};

template <class TriangleMesh,
          class Traits_,
          class VertexPointMap_,
          class SolverTraits_>
std::size_t Mean_curvature_flow_skeletonization<TriangleMesh, Traits_, VertexPointMap_, SolverTraits_>::collapse_short_edges()
{
  std::size_t cnt=0, prev_cnt=0;

  std::set<edge_descriptor> edges_to_collapse, non_topologically_valid_collapses;

  BOOST_FOREACH(edge_descriptor ed, edges(m_tmesh))
    if ( edge_should_be_collapsed(ed) )
      edges_to_collapse.insert(ed);

  do{
    prev_cnt=cnt;
    while(!edges_to_collapse.empty())
    {
      edge_descriptor ed = *edges_to_collapse.begin();
      edges_to_collapse.erase(edges_to_collapse.begin());

      // skip the edge is it became long enough
      if ( !edge_should_be_collapsed(ed) ) continue;

      if ( !Euler::does_satisfy_link_condition(ed,m_tmesh) )
      {
        non_topologically_valid_collapses.insert(ed);
        continue;
      }

      halfedge_descriptor h = halfedge(ed, m_tmesh);

      vertex_descriptor vi = source(h, m_tmesh);
      vertex_descriptor vj = target(h, m_tmesh);

      Point p = m_traits.construct_midpoint_3_object()(
                  get(vertex_point, m_tmesh, vi),
                  get(vertex_point, m_tmesh, vj) );

      // invalidate the edges that will be collapsed
      edges_to_collapse.erase(edge(prev(h, m_tmesh), m_tmesh));
      edges_to_collapse.erase(edge(prev(opposite(h, m_tmesh), m_tmesh), m_tmesh));

      non_topologically_valid_collapses.erase(ed);
      non_topologically_valid_collapses.erase(edge(prev(h, m_tmesh), m_tmesh));
      non_topologically_valid_collapses.erase(edge(prev(opposite(h, m_tmesh), m_tmesh), m_tmesh));

      // the mesh is closed, the target of h is always the one kept
      put(m_tmesh_point_pmap, vj, p);
      std::vector<Input_vertex_descriptor>& vec_kept = vj->vertices;
      std::vector<Input_vertex_descriptor>& vec_removed = vi->vertices;
      vec_kept.insert(vec_kept.end(), vec_removed.begin(), vec_removed.end());
      if (vi->is_fixed) vj->is_fixed=true;
      update_pole(vi, vj);

      vertex_descriptor v = Euler::collapse_edge(ed, m_tmesh);

      CGAL_assertion(vj==v);

      BOOST_FOREACH(edge_descriptor oed, out_edges(v, m_tmesh))
        if ( edge_should_be_collapsed(oed) ) edges_to_collapse.insert(oed);

      ++cnt;
    }
    if (prev_cnt==cnt) break;
    edges_to_collapse.swap(non_topologically_valid_collapses);
  } while(!edges_to_collapse.empty());

  return cnt;
}


} //namespace CGAL

#endif // CGAL_MEAN_CURVATURE_FLOW_SKELETONIZATION_H
