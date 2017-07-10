// Copyright (c) 2014 GeometryFactory
// All rights reserved.
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
// Author(s)     : Yin Xu, Andreas Fabri and Ilker O. Yaz

#ifndef CGAL_SURFACE_MESH_DEFORMATION_H
#define CGAL_SURFACE_MESH_DEFORMATION_H

#include <CGAL/license/Surface_mesh_deformation.h>


#include <CGAL/config.h>
#include <CGAL/Default.h>
#include <CGAL/tuple.h>

#include <CGAL/Polygon_mesh_processing/Weights.h>
#include <CGAL/Simple_cartesian.h>

#include <vector>
#include <list>
#include <utility>
#include <limits>
#include <boost/foreach.hpp>

/*
#define CGAL_DEFORM_MESH_USE_EXPERIMENTAL_SCALE // define it to activate optimal scale calculation,
// then you can define below to just scale, if not both rotate and scale will be activated
#define CGAL_DEFORM_MESH_JUST_EXPERIMENTAL_SCALE // to not to rotate but just scale
*/

// for default parameters
#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>  // for sparse linear system solver
#include <CGAL/Deformation_Eigen_polar_closest_rotation_traits_3.h>  // for 3x3 closest rotation computer
#endif

namespace CGAL {

/// \ingroup PkgSurfaceMeshDeformation
///@brief Deformation algorithm type
enum Deformation_algorithm_tag
{
  ORIGINAL_ARAP,   /**< use original as-rigid-as possible algorithm */
  SPOKES_AND_RIMS, /**< use spokes and rims version of as-rigid-as possible algorithm */
  SRE_ARAP         /**< use smooth rotation enhanced As-rigid-as-possible */
};

/// @cond CGAL_DOCUMENT_INTERNAL
namespace internal {
template<class TriangleMesh, Deformation_algorithm_tag deformation_algorithm_tag>
struct Types_selectors;

template<class TriangleMesh>
struct Types_selectors<TriangleMesh, CGAL::SPOKES_AND_RIMS> {
  typedef internal::Single_cotangent_weight_impl<TriangleMesh> Weight_calculator;

  struct ARAP_visitor{
    template <class VertexPointMap>
    void init(TriangleMesh, VertexPointMap){}

    void rotation_matrix_pre(
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor,
      TriangleMesh&){}

    template <class Square_matrix_3>
    void update_covariance_matrix(
      Square_matrix_3&,
      const Square_matrix_3&){}

    void set_sre_arap_alpha(double){}
  };
};

template<class TriangleMesh>
struct Types_selectors<TriangleMesh, CGAL::ORIGINAL_ARAP> {
  typedef internal::Cotangent_weight_impl<TriangleMesh> Weight_calculator;

  typedef typename Types_selectors<TriangleMesh, CGAL::SPOKES_AND_RIMS>
    ::ARAP_visitor ARAP_visitor;
};

template<class TriangleMesh>
struct Types_selectors<TriangleMesh, CGAL::SRE_ARAP> {
  typedef internal::Cotangent_weight_impl<TriangleMesh> Weight_calculator;

  class ARAP_visitor{
    double m_nb_edges_incident;
    double m_area;
    double m_alpha;

  public:
    ARAP_visitor(): m_alpha(0.02) {}

    template<class VertexPointMap>
    void init(TriangleMesh triangle_mesh, const VertexPointMap& vpmap)
    {
      // calculate area
      m_area = 0;
      typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
      BOOST_FOREACH(face_descriptor f, faces(triangle_mesh))
      {
        typename boost::graph_traits<TriangleMesh>::halfedge_descriptor
          h = halfedge(f, triangle_mesh);
        m_area += std::sqrt(CGAL::squared_area(
                        get(vpmap, source(h, triangle_mesh) ),
                        get(vpmap, target(h, triangle_mesh) ),
                        get(vpmap, target(next(h, triangle_mesh), triangle_mesh) ) ));
      }
    }

    void rotation_matrix_pre(
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor vi,
      TriangleMesh& hg)
    {
      typename boost::graph_traits<TriangleMesh>::in_edge_iterator e, e_end;
      cpp11::tie(e,e_end) = in_edges(vi, hg);
      m_nb_edges_incident=(double) std::distance(e,e_end);
    }

    template <class Square_matrix_3>
    void update_covariance_matrix(
      Square_matrix_3& cov,
      const Square_matrix_3& rot_mtr)
    {
      // add neighbor rotation
      cov += m_alpha * m_area * rot_mtr.transpose() / m_nb_edges_incident;
    }

    void set_sre_arap_alpha(double a){ m_alpha=a; }
  };
};

// property map that create a Simple_cartesian<double>::Point_3
// on the fly in order the deformation class to be used
// with points with minimal requirements
template <class Vertex_point_map>
struct SC_on_the_fly_pmap: public Vertex_point_map{
  typedef boost::readable_property_map_tag category;
  typedef CGAL::Simple_cartesian<double>::Point_3 value_type;
  typedef value_type reference;
  typedef typename boost::property_traits<Vertex_point_map>::key_type key_type;

  SC_on_the_fly_pmap(Vertex_point_map base):
    Vertex_point_map(base) {}

  friend value_type
  get(const SC_on_the_fly_pmap map, key_type k)
  {
    typename boost::property_traits<Vertex_point_map>::reference base=
      get(static_cast<const Vertex_point_map&>(map), k);
    return value_type(base[0], base[1], base[2]);
  }
};


}//namespace internal
/// @endcond

 ///
 /// \ingroup PkgSurfaceMeshDeformation
 /// @brief Class providing the functionalities for deforming a triangulated surface mesh
 ///
 /// @tparam TM a model of `HalfedgeGraph`
 /// @tparam VIM a model of `ReadablePropertyMap` with `Surface_mesh_deformation::vertex_descriptor` as key and `unsigned int` as value type.
 ///         The default is `boost::property_map<TM, boost::%vertex_index_t>::%type`.
 /// @tparam HIM a model of `ReadablePropertyMap` with `Surface_mesh_deformation::halfedge_descriptor` as key and `unsigned int` as value type.
 ///         The default is `boost::property_map<TM, boost::%halfedge_index_t>::%type`.
 /// @tparam TAG tag for selecting the deformation algorithm
 /// @tparam WC a model of SurfaceMeshDeformationWeights, with `WC::Triangle_mesh` being `TM`.
 ///         If `TAG` is `ORIGINAL_ARAP`, the weights must be positive to guarantee a correct energy minimization.
 ///         The default is the cotangent weighting scheme. In case `TAG` is `ORIGINAL_ARAP`, negative weights are clamped to zero.
 /// @tparam ST a model of SparseLinearAlgebraWithFactorTraits_d. If \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available
 /// and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits` is provided as default parameter.\n
  /// \code
 ///     CGAL::Eigen_solver_traits<
 ///         Eigen::SparseLU<
 ///            CGAL::Eigen_sparse_matrix<double>::EigenType,
 ///            Eigen::COLAMDOrdering<int> >  >
 /// \endcode
 /// @tparam CR a model of DeformationClosestRotationTraits_3. If \ref thirdpartyEigen "Eigen" 3.1 (or greater) is available and `CGAL_EIGEN3_ENABLED` is defined,
 /// `Deformation_Eigen_polar_closest_rotation_traits_3` is provided as default parameter.
 /// @tparam VPM a model of `ReadWritePropertyMap`</a>  with `Surface_mesh_deformation::vertex_descriptor` as key and a point as value type. The point type must be a model of `::RawPoint_3`.
 /// The default is `boost::property_map<TM, CGAL::vertex_point_t>::%type`.
template <
  class TM,
  class VIM=Default,
  class HIM=Default,
  Deformation_algorithm_tag TAG = SPOKES_AND_RIMS,
  class WC = Default,
  class ST = Default,
  class CR = Default,
  class VPM = Default
  >
class Surface_mesh_deformation
{
//Typedefs
public:

  /// \name Types
  /// @{
  // typedefed template parameters, main reason is doxygen creates autolink to typedefs but not template parameters
  ///
  typedef TM Triangle_mesh;
  typedef TM Halfedge_graph;

// Index maps
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    VIM,
    typename boost::property_map<Triangle_mesh, boost::vertex_index_t>::type
  >::type Vertex_index_map;
  typedef typename Default::Get<
    HIM,
    typename boost::property_map<Triangle_mesh, boost::halfedge_index_t>::type
  >::type Hedge_index_map;
#else
  ///
  typedef VIM Vertex_index_map;
  ///
  typedef HIM Hedge_index_map;
#endif

// weight calculator
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    WC,
    typename internal::Types_selectors<TM, TAG>::Weight_calculator
  >::type Weight_calculator;
#else
  ///
  typedef WC Weight_calculator;
#endif

// sparse linear solver
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    ST,
  #if defined(CGAL_EIGEN3_ENABLED)
      CGAL::Eigen_solver_traits<
          Eigen::SparseLU<
            CGAL::Eigen_sparse_matrix<double>::EigenType,
            Eigen::COLAMDOrdering<int> >  >
  #else
    ST // no parameter provided, and Eigen is not enabled: so don't compile!
  #endif
  >::type Sparse_linear_solver;
#else
  ///
  typedef ST Sparse_linear_solver;
#endif

// CR helper
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    CR,
  #if defined(CGAL_EIGEN3_ENABLED)
    Deformation_Eigen_polar_closest_rotation_traits_3
  #else
    CR // no parameter provided, and Eigen is not enabled: so don't compile!
  #endif
  >::type Closest_rotation_traits;
#else
  ///
  typedef CR Closest_rotation_traits;
#endif

// vertex point pmap
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    VPM,
    typename boost::property_map<Triangle_mesh, CGAL::vertex_point_t>::type
  >::type Vertex_point_map;
#else
  ///
  typedef VPM Vertex_point_map;
#endif

  /// The type for vertex descriptor
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
  /// The type for halfedge descriptor
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor halfedge_descriptor;
  /// The 3D point type, model of `::RawPoint_3`
  typedef typename boost::property_traits<Vertex_point_map>::value_type Point;
  /// A constant iterator range over the vertices of the region-of-interest.
  /// It is a model of `ConstRange` with `vertex_descriptor` as iterator value type.
  typedef std::vector<vertex_descriptor> Roi_vertex_range;
/// @}

private:
  typedef Surface_mesh_deformation<TM, VIM, HIM, TAG, WC, ST, CR> Self;
  // Repeat Triangle_mesh types
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_iterator     vertex_iterator;
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_iterator       halfedge_iterator;
  typedef typename boost::graph_traits<Triangle_mesh>::in_edge_iterator    in_edge_iterator;
  typedef typename boost::graph_traits<Triangle_mesh>::out_edge_iterator   out_edge_iterator;

  typedef typename Closest_rotation_traits::Matrix CR_matrix;
  typedef typename Closest_rotation_traits::Vector CR_vector;

// Data members.
  Triangle_mesh& m_triangle_mesh;                   /**< Source triangulated surface mesh for modeling */

  std::vector<Point> original;                        ///< original positions of roi (size: ros + boundary_of_ros)
  std::vector<Point> solution;                        ///< storing position of ros vertices during iterations (size: ros + boundary_of_ros)

  Vertex_index_map vertex_index_map;                  ///< storing indices of all vertices
  Hedge_index_map   hedge_index_map;                  ///< storing indices of all halfedges

  std::vector<vertex_descriptor> roi;                 ///< region of interest
  std::vector<vertex_descriptor> ros;                 ///< region of solution, including roi and hard constraints on boundary of roi

  std::vector<std::size_t> ros_id_map;                ///< (size: num vertices)
  std::vector<bool>        is_roi_map;                ///< (size: num vertices)
  std::vector<bool>        is_ctrl_map;               ///< (size: num vertices)

  std::vector<double> hedge_weight;                   ///< all halfedge weights
  std::vector<CR_matrix> rot_mtr;                     ///< rotation matrices of ros vertices (size: ros)

  Sparse_linear_solver m_solver;                      ///< linear sparse solver
  unsigned int m_iterations;                          ///< number of maximal iterations
  double m_tolerance;                                 ///< tolerance of convergence

  bool need_preprocess_factorization;                 ///< is there any need to compute L and factorize
  bool need_preprocess_region_of_solution;            ///< is there any need to compute region of solution

  bool last_preprocess_successful;                    ///< stores the result of last call to preprocess()

  Weight_calculator weight_calculator;

  Vertex_point_map vertex_point_map;

public:
  typename internal::Types_selectors<TM, TAG>::ARAP_visitor arap_visitor;
private:

#ifdef CGAL_DEFORM_MESH_USE_EXPERIMENTAL_SCALE
  std::vector<double> scales;
#endif

#ifndef CGAL_CFG_NO_CPP0X_DELETED_AND_DEFAULT_FUNCTIONS
public:
  Surface_mesh_deformation(const Self&) = delete; // no copy
#else
private:
  Surface_mesh_deformation(const Self&); // no copy
#endif


// Public methods
public:

  /// \cond SKIP_FROM_MANUAL
  //vertex_point_map set by default
  Surface_mesh_deformation(Triangle_mesh& triangle_mesh,
                           Vertex_index_map vertex_index_map,
                           Hedge_index_map hedge_index_map
                          )
    : m_triangle_mesh(triangle_mesh), vertex_index_map(vertex_index_map), hedge_index_map(hedge_index_map),
      ros_id_map(std::vector<std::size_t>(num_vertices(triangle_mesh), (std::numeric_limits<std::size_t>::max)() )),
      is_roi_map(std::vector<bool>(num_vertices(triangle_mesh), false)),
      is_ctrl_map(std::vector<bool>(num_vertices(triangle_mesh), false)),
      m_iterations(5), m_tolerance(1e-4),
      need_preprocess_factorization(true),
      need_preprocess_region_of_solution(true),
      last_preprocess_successful(false),
      weight_calculator(Weight_calculator()),
      vertex_point_map(get(vertex_point, triangle_mesh))
  {
    init();
  }

  //vertex_point_map and hedge_index_map set by default
  Surface_mesh_deformation(Triangle_mesh& triangle_mesh,
                           Vertex_index_map vertex_index_map
                          )
    : m_triangle_mesh(triangle_mesh), vertex_index_map(vertex_index_map),
      hedge_index_map(get(boost::halfedge_index, triangle_mesh)),
      ros_id_map(std::vector<std::size_t>(num_vertices(triangle_mesh), (std::numeric_limits<std::size_t>::max)() )),
      is_roi_map(std::vector<bool>(num_vertices(triangle_mesh), false)),
      is_ctrl_map(std::vector<bool>(num_vertices(triangle_mesh), false)),
      m_iterations(5), m_tolerance(1e-4),
      need_preprocess_factorization(true),
      need_preprocess_region_of_solution(true),
      last_preprocess_successful(false),
      weight_calculator(Weight_calculator()),
      vertex_point_map(get(vertex_point, triangle_mesh))
  {
    init();
  }
  //vertex_point_map, hedge_index_map and vertex_index_map set by default
  Surface_mesh_deformation(Triangle_mesh& triangle_mesh)
    : m_triangle_mesh(triangle_mesh),
      vertex_index_map(get(boost::vertex_index, triangle_mesh)),
      hedge_index_map(get(boost::halfedge_index, triangle_mesh)),
      ros_id_map(std::vector<std::size_t>(num_vertices(triangle_mesh), (std::numeric_limits<std::size_t>::max)() )),
      is_roi_map(std::vector<bool>(num_vertices(triangle_mesh), false)),
      is_ctrl_map(std::vector<bool>(num_vertices(triangle_mesh), false)),
      m_iterations(5), m_tolerance(1e-4),
      need_preprocess_factorization(true),
      need_preprocess_region_of_solution(true),
      last_preprocess_successful(false),
      weight_calculator(Weight_calculator()),
      vertex_point_map(get(vertex_point, triangle_mesh))
  {
    init();
  }

  // Constructor with all the parameters provided
  Surface_mesh_deformation(Triangle_mesh& triangle_mesh,
                           Vertex_index_map vertex_index_map,
                           Hedge_index_map hedge_index_map,
                           Vertex_point_map vertex_point_map,
                           Weight_calculator weight_calculator = Weight_calculator()
                          )
    : m_triangle_mesh(triangle_mesh), vertex_index_map(vertex_index_map), hedge_index_map(hedge_index_map),
    ros_id_map(std::vector<std::size_t>(num_vertices(triangle_mesh), (std::numeric_limits<std::size_t>::max)() )),
    is_roi_map(std::vector<bool>(num_vertices(triangle_mesh), false)),
    is_ctrl_map(std::vector<bool>(num_vertices(triangle_mesh), false)),
    m_iterations(5), m_tolerance(1e-4),
    need_preprocess_factorization(true),
    need_preprocess_region_of_solution(true),
    last_preprocess_successful(false),
    weight_calculator(weight_calculator),
    vertex_point_map(vertex_point_map)
  {
    init();
  }
  /// \endcond
  #if DOXYGEN_RUNNING
/// \name Construction
/// @{
    /**
   * The constructor of a deformation object
   *
   * @pre `triangle_mesh` consists of only triangular facets
   * @param triangle_mesh triangulated surface mesh to deform
   * @param vertex_index_map property map which associates an id to each vertex, from `0` to `num_vertices(triangle_mesh)-1`.
   * @param hedge_index_map property map which associates an id to each halfedge, from `0` to `2*num_edges(triangle_mesh)-1`.
   * @param vertex_point_map property map which associates a point to each vertex of the triangle mesh.
   * @param weight_calculator function object or pointer for weight calculation
   */
  Surface_mesh_deformation(Triangle_mesh& triangle_mesh,
    Vertex_index_map vertex_index_map=get(boost::vertex_index, triangle_mesh),
    Hedge_index_map hedge_index_map=get(boost::halfedge_index, triangle_mesh),
    Vertex_point_map vertex_point_map=get(boost::vertex_point, triangle_mesh),
    Weight_calculator weight_calculator = Weight_calculator()
    );
/// @}
  #endif

private:
  void init() {
    typedef internal::SC_on_the_fly_pmap<Vertex_point_map> Wrapper;
    // compute halfedge weights
    halfedge_iterator eb, ee;
    hedge_weight.reserve(2*num_edges(m_triangle_mesh));
    for(cpp11::tie(eb, ee) = halfedges(m_triangle_mesh); eb != ee; ++eb)
    {
      hedge_weight.push_back(
        this->weight_calculator(*eb, m_triangle_mesh, Wrapper(vertex_point_map)));
    }
    arap_visitor.init(m_triangle_mesh, vertex_point_map);
  }

public:

/// \name Preprocessing
/// @{
  /**
   * Erases all the vertices from the region-of-interest (control vertices included).
   */
  void clear_roi_vertices(){
    need_preprocess_both();
    // clear roi vertices
    roi.clear();
    //set to false all bits
    is_roi_map.assign(num_vertices(m_triangle_mesh), false);
    is_ctrl_map.assign(num_vertices(m_triangle_mesh), false);
  }

  /**
   * Erases all the vertices from the set of control vertices.
   */
  void clear_control_vertices(){
    need_preprocess_factorization=true;
    //set to false all bits
    is_ctrl_map.assign(num_vertices(m_triangle_mesh), false);
  }

  /**
   * Inserts a vertex in the set of control vertices. The vertex is also inserted in the region-of-interest if it is not already in it.
   * @param vd the vertex to be inserted
   * @return `true` if `vd` is not already a control vertex.
   */
  bool insert_control_vertex(vertex_descriptor vd)
  {
    if(is_control_vertex(vd)) { return false; }
    need_preprocess_factorization=true;

    insert_roi_vertex(vd); // also insert it as roi

    is_ctrl_map[id(vd)] = true;
    return true;
  }

  /**
   * Inserts a range of vertices in the set of control vertices. The vertices are also inserted in the region-of-interest if they are not already in it.
   * @tparam InputIterator input iterator type with `vertex_descriptor` as value type
   * @param begin first iterator of the range of vertices
   * @param end past-the-end iterator of the range of vertices
   */
  template<class InputIterator>
  void insert_control_vertices(InputIterator begin, InputIterator end)
  {
    for( ;begin != end; ++begin)
    {
      insert_control_vertex(*begin);
    }
  }

  /**
   * Erases a vertex from the set of control vertices.
   * @param vd the vertex to be erased
   * @return `true` if `vd` was a control vertex.
   */
  bool erase_control_vertex(vertex_descriptor vd)
  {
    if(!is_control_vertex(vd)) { return false; }

    need_preprocess_factorization=true;
    is_ctrl_map[id(vd)] = false;
    return true;
  }

  /**
   * Inserts a range of vertices in the region-of-interest
   * @tparam InputIterator input iterator with `vertex_descriptor` as value type
   * @param begin first iterator of the range of vertices
   * @param end past-the-end iterator of the range of vertices
   */
  template<class InputIterator>
  void insert_roi_vertices(InputIterator begin, InputIterator end)
  {
    for( ;begin != end; ++begin)
    {
      insert_roi_vertex(*begin);
    }
  }

  /**
   * Inserts a vertex in the region-of-interest
   * @param vd the vertex to be inserted
   * @return `true` if `vd` is not already in the region-of-interest.
   */
  bool insert_roi_vertex(vertex_descriptor vd)
  {
    if(is_roi_vertex(vd)) { return false; }
    need_preprocess_both();

    is_roi_map[id(vd)] = true;
    roi.push_back(vd);
    return true;
  }

  /**
   * Erases a vertex from the region-of-interest and the set of control vertices.
   * \note At the next call to `preprocess()`, any vertex that is no longer in the region-of-interest will be assigned to its original position
   * (that is the position of the vertex at the time of construction or after the last call to `overwrite_initial_geometry()`).
   * @param vd the vertex to be erased
   * @return `true` `vd` was a vertex from the region-of-interest.
   */
  bool erase_roi_vertex(vertex_descriptor vd)
  {
    if(!is_roi_vertex(vd)) { return false; }

    erase_control_vertex(vd); // also erase from being control

    typename std::vector<vertex_descriptor>::iterator it = std::find(roi.begin(), roi.end(), vd);
    if(it != roi.end())
    {
      is_roi_map[id(vd)] = false;
      roi.erase(it);

      need_preprocess_both();
      return true;
    }

    CGAL_assertion(false); // inconsistency between is_roi_map, and roi vector!
    return false;
  }

  /**
   * Preprocessing function that need to be called each time the region-of-interest or the set
   * of control vertices are changed before calling `deform()`.
   * If not already done, `deform()` first calls this function.
   * \cgalAdvancedBegin
   * Collects the vertices not in the region-of-interest that are adjacent to a vertex from the
   * region-of-interest (these vertices are internally considered as fixed control vertices).
   * Then assembles and factorizes the Laplacian matrix used in the function `deform()`.
   * \cgalAdvancedEnd
   * \note A modification of the set of control vertices or the region-of-interest invalidates the
   * preprocessing data.
   * @return `true` if successful.
   * A common reason for failure is that the system is rank deficient,
   * which happens for example when all the vertices are in the region-of-interest and no control vertices are set, or
   * if the weighting scheme used features too many zero and breaks the connectivity information.
   */
  bool preprocess()
  {
    region_of_solution();
    assemble_laplacian_and_factorize();
    return last_preprocess_successful; // which is set by assemble_laplacian_and_factorize()
  }
/// @} Preprocessing

/// \name Deformation
/// @{

  /**
   * Sets the target position of a control vertex.
   * @param vd the control vertex the target position is set
   * @param target_position the new target position
   */
  void set_target_position(vertex_descriptor vd, const Point& target_position)
  {
    region_of_solution(); // we require ros ids, so if there is any need to preprocess of region of solution -do it.

    if(!is_control_vertex(vd)) { return; }
    solution[ros_id(vd)] = target_position;
  }

  /**
   * Updates the target position of `vd` by applying     a translation of vector `t`.
   *
   * @tparam Vect is a 3D vector class, with `Vect(double x,double y, double z)` being a constructor from its %Cartesian coordinates
   *         and `double Vect::operator[](int i)` with i=0,1 or 2 returning its %Cartesian coordinates.
   *
   * @param vd a control vertex
   * @param t translation vector
   * \pre `is_control_vertex(vd)`
   */
  template<class Vect>
  void translate(vertex_descriptor vd, const Vect& t)
  {
    region_of_solution(); // we require ros ids, so if there is any need to preprocess of region of solution -do it.

    std::size_t v_id = ros_id(vd);
    solution[v_id] = add_to_point(solution[v_id], t);
  }

  /**
   * Equivalent to calling the overload taking only one control vertex, for each vertex in the range `[begin,end[`.
   *
   * @tparam InputIterator input iterator type with `vertex_descriptor` as value type
   * @tparam Vect is a 3D vector class, with `Vect(double x,double y, double z)` being a constructor from its %Cartesian coordinates
   *         and `double Vect::operator[](int i)` with i=0,1 or 2 returning its %Cartesian coordinates.
   *
   * @param begin first iterator of the range of vertices
   * @param end past-the-end iterator of the range of vertices
   * @param t translation vector
   */
  template<class InputIterator, class Vect>
  void translate(InputIterator begin, InputIterator end, const Vect& t)
  {
    region_of_solution(); // we require ros ids, so if there is any need to preprocess of region of solution -do it.

    for(; begin != end; ++begin) {
      std::size_t v_id = ros_id(*begin);
      solution[v_id] = add_to_point(solution[v_id], t);
    }
  }

  /**
   * Updates the target position of `vd` by applying to its last target position
   * a rotation defined by the quaternion `quat`, the center of the rotation being
   * the origin translated by `to_rotation_center` .
   *
   * @tparam Quaternion is a quaternion class with `Vect operator*(Quaternion, Vect)` returning the product of a quaternion with a vector
   * @tparam Vect is a 3D vector class, with `Vect(double x,double y, double z)` being a constructor from its %Cartesian coordinates
   *         and `double Vect::operator[](int i)` with i=0,1 or 2 returning its %Cartesian coordinates.
   *
   * @param vd a control vertex
   * @param to_rotation_center the vector to translate the origin to the center of the rotation
   * @param quat quaternion of the rotation
   * \pre `is_control_vertex(vd)`
   * \pre `quad` represents a rotation
   */
  template <typename Quaternion, typename Vect>
  void rotate(vertex_descriptor vd, const Vect& to_rotation_center, const Quaternion& quat)
  {
    CGAL_precondition( is_control_vertex(vd) );
    region_of_solution(); // we require ros ids, so if there is any need to preprocess of region of solution -do it.

    std::size_t v_id = ros_id(vd);
    Vect v = quat * sub_to_vect(solution[v_id], to_rotation_center);
    solution[v_id] = Point( to_rotation_center[0]+v[0],
                            to_rotation_center[1]+v[1],
                            to_rotation_center[2]+v[2] );
  }

  /**
   * Equivalent to calling the overload taking only one control vertex, for each vertex in the range `[begin,end[`.
   *
   * @tparam InputIterator input iterator type with `vertex_descriptor` as value type
   * @tparam Quaternion is a quaternion class with `Vect operator*(Quaternion, Vect)` returning the product of a quaternion with a vector
   * @tparam Vect is a 3D vector class, with `Vect(double x,double y, double z)` being a constructor from its %Cartesian coordinates
   *         and `double Vect::operator[](int i)` with i=0,1 or 2 returning its %Cartesian coordinates.
   *
   * @param begin first iterator of the range of vertices
   * @param end past-the-end iterator of the range of vertices
   * @param to_rotation_center the vector to translate the origin to the center of the rotation
   * @param quat quaternion of the rotation
   * \pre `quad` represents a rotation
   */
  template <typename InputIterator, typename Vect, typename Quaternion>
  void rotate(InputIterator begin, InputIterator end, const Vect& to_rotation_center, const Quaternion& quat)
  {
    region_of_solution(); // we require ros ids, so if there is any need to preprocess of region of solution -do it.

    for(; begin != end; ++begin) {
      std::size_t v_id = ros_id(*begin);
      Vect v = quat * sub_to_vect(solution[v_id], to_rotation_center);
      solution[v_id] = Point( to_rotation_center[0]+v[0],
                              to_rotation_center[1]+v[1],
                              to_rotation_center[2]+v[2] );
    }
  }

  /**
    * Returns the target position of a control vertex.
    * \param vd a control vertex
    * \pre `is_control_vertex(vd)`
    */
  const Point& target_position(vertex_descriptor vd)
  {
    region_of_solution();

    CGAL_precondition( is_control_vertex(vd) );
    return solution[ ros_id(vd) ];
  }

  /**
   * Deforms the region-of-interest according to the deformation algorithm, using the target positions of each control vertex set by using `rotate()`, `translate()`, or `set_target_position()`.
   * The points associated to each vertex of the input triangle mesh that are inside the region-of-interest are updated.
   * \note Nothing happens if `preprocess()` returns `false`.
   * @see set_iterations(unsigned int iterations), set_tolerance(double tolerance), deform(unsigned int iterations, double tolerance)
   */
  void deform()
  {
    deform(m_iterations, m_tolerance);
  }

  /**
   * Same as `deform()` but the number of iterations and the tolerance are one-time parameters.
   * @param iterations number of iterations for optimization procedure
   * @param tolerance tolerance of convergence (see explanations set_tolerance(double tolerance))
   */
  void deform(unsigned int iterations, double tolerance)
  {
    preprocess();

    if(!last_preprocess_successful) {
      CGAL_warning(false);
      return;
    }
    // Note: no energy based termination occurs at first iteration
    // because comparing energy of original model (before deformation) and deformed model (deformed_1_iteration)
    // simply does not make sense, comparison is meaningful between deformed_(i)_iteration & deformed_(i+1)_iteration

    double energy_this = 0; // initial value is not important, because we skip first iteration
    double energy_last;

    // iterations
    for ( unsigned int ite = 0; ite < iterations; ++ite)
    {
      // main steps of optimization
      update_solution();
#ifndef CGAL_DEFORM_MESH_JUST_EXPERIMENTAL_SCALE
      optimal_rotations();
#endif
#ifdef CGAL_DEFORM_MESH_USE_EXPERIMENTAL_SCALE
      optimal_scales();
#endif
      // energy based termination
      if(tolerance > 0.0 && (ite + 1) < iterations) // if tolerance <= 0 then don't compute energy
      {                                             // also no need compute energy if this iteration is the last iteration
        energy_last = energy_this;
        energy_this = energy();
        CGAL_warning(energy_this >= 0);

        if(ite != 0) // skip first iteration
        {
          double energy_dif = std::abs((energy_last - energy_this) / energy_this);
          if ( energy_dif < tolerance ) { break; }
        }
      }
    }
    // copy solution to the target surface mesh
    assign_solution();
  }

  /**
   * Resets the points associated to the vertices of the region-of-interest at their
   * initial positions at time of the functor construction or after
   * the last call to `overwrite_initial_geometry()`.
   * \note if the region-of-interest or the set of control vertices have been
   * modified since the last call to `preprocess()`, it will be called prior
   * to the reset.
   */
  void reset()
  {
    region_of_solution(); // since we are using original vector

    //restore the current positions to be the original positions
    BOOST_FOREACH(vertex_descriptor vd, roi_vertices())
    {
      put(vertex_point_map, vd, original[ros_id(vd)]);
      solution[ros_id(vd)]=original[ros_id(vd)];
    }

    // also set rotation matrix to identity
    std::fill(rot_mtr.begin(), rot_mtr.end(),
              Closest_rotation_traits().identity_matrix());
  }

  /**
   * Sets the initial positions of the vertices from the region-of-interest to the current positions. Calling this function has the same effect as creating
   * a new deformation object with the current deformed triangle mesh, keeping the region-of-interest and the set of control vertices.
   * \note if the region-of-interest or the set of control vertices have been modified since the last call to `preprocess()`,
   * it will be called prior to the overwrite.
   *
   * \cgalAdvancedBegin
   * This function might have a non-negligible effect on the result.
   * The Laplacian matrix of the free vertices and the optimal rotations
   * are computed using the original positions of the points
   * associated to the vertices. Thus, if a deformed version of the surface mesh
   * is used as reference, the surface mesh properties the algorithm
   * tries to preserve are those of an altered version (which are already
   * degraded).
   * \cgalAdvancedEnd
   */
  void overwrite_initial_geometry()
  {
    typedef internal::SC_on_the_fly_pmap<Vertex_point_map> Wrapper;
    if(roi.empty()) { return; } // no ROI to overwrite

    region_of_solution(); // the roi should be preprocessed since we are using original_position vec

    BOOST_FOREACH(vertex_descriptor vd, roi_vertices())
    {
      original[ros_id(vd)] = get(vertex_point_map, vd);
    }

    // now I need to compute weights for halfedges incident to roi vertices
    std::vector<bool> is_weight_computed(2*num_edges(m_triangle_mesh), false);
    BOOST_FOREACH(vertex_descriptor vd, roi_vertices())
    {
      in_edge_iterator e, e_end;
      for (cpp11::tie(e,e_end) = in_edges(vd, m_triangle_mesh); e != e_end; e++)
      {
        halfedge_descriptor he = halfedge(*e, m_triangle_mesh);
        std::size_t id_e = id(he);
        if(is_weight_computed[id_e]) { continue; }

        hedge_weight[id_e] = weight_calculator(he, m_triangle_mesh, Wrapper(vertex_point_map));
        is_weight_computed[id_e] = true;

        halfedge_descriptor e_opp = opposite(he, m_triangle_mesh);
        std::size_t id_e_opp = id(e_opp);

        hedge_weight[id_e_opp] = weight_calculator(e_opp, m_triangle_mesh, Wrapper(vertex_point_map));
        is_weight_computed[id_e_opp] = true;
      }
    }

    // also set rotation matrix to identity
    std::fill(rot_mtr.begin(), rot_mtr.end(), Closest_rotation_traits().identity_matrix());

    need_preprocess_both(); // now we need reprocess
  }

/// @} Deformation

/// \name Utilities
/// @{

  /**
   * Gets the default number of iterations (5) or the value passed to the function `set_iterations()`
   */
  unsigned int iterations()
  { return m_iterations; }

  /**
   * Gets the default tolerance parameter (1e-4) or the value passed to the function `set_tolerance()`
   */
  double tolerance()
  { return m_tolerance; }

  /**
   * Sets the maximum number of iterations ran by `deform()`
   */
  void set_iterations(unsigned int iterations)
  { this->m_iterations = iterations; }

   /// @brief Sets the tolerance of the convergence used in `deform()`.
   /// If `tolerance==0`, no energy based termination criteria is used (preventing to do the energy computation at each iteration step)
   ///
   /// `tolerance >` \f$|\mathrm{energy}(m_i) - \mathrm{energy}(m_{i-1})| / \mathrm{energy}(m_i)\f$ will be used as a termination criterium.
  void set_tolerance(double tolerance)
  { this->m_tolerance = tolerance; }

  /**
   * Returns the range of vertices in the region-of-interest.
   */
  const Roi_vertex_range& roi_vertices() const
  {
    return roi;
  }

  /**
   * Tests whether a vertex is inside the region-of-interest.
   * @param vd the query vertex
   * @return `true` if `vd` has been inserted to (and not erased from) the region-of-interest.
   */
  bool is_roi_vertex(vertex_descriptor vd) const
  { return is_roi_map[id(vd)]; }

  /**
   * Tests whether a vertex is a control vertex.
   * @param vd the query vertex
   * @return `true` if `vd` has been inserted to (and not erased from) the set of control vertices.
   */
  bool is_control_vertex(vertex_descriptor vd) const
  { return is_ctrl_map[id(vd)]; }

  /**
   * Provides access to the triangle mesh being deformed
   */
  const Triangle_mesh& triangle_mesh() const
  { return m_triangle_mesh; }

  const Triangle_mesh& halfedge_graph() const
  { return m_triangle_mesh; }

  /**
   * Sets the alpha coefficient that determines the weight of the bending term (rotation smoothness) for the SRE-ARAP deformation technique.
   * The range of values can be from 0 to infinity. When alpha=0, the method reverts to ARAP. When alpha is increased, neighboring rotations become similar to each other, where alpha=infinity results in a global rigid transformation of the ROI. The value of alpha depends on the surface area and shape. Since alpha is not too sensitive, it can be tweaked in most cases by multiplying it by powers of 10.
   * The default value for alpha is 0.02.
   */
  void set_sre_arap_alpha(double a)
  {
    arap_visitor.set_sre_arap_alpha(a);
  }

/// @} Utilities


private:

  /// Assigns id to one ring neighbor of vd, and also push them into push_vector
  void assign_ros_id_to_one_ring(vertex_descriptor vd,
                             std::size_t& next_id,
                             std::vector<vertex_descriptor>& push_vector)
  {
    in_edge_iterator e, e_end;
    for (cpp11::tie(e,e_end) = in_edges(vd, m_triangle_mesh); e != e_end; e++)
    {
      vertex_descriptor vt = source(*e, m_triangle_mesh);
      if(ros_id(vt) == (std::numeric_limits<std::size_t>::max)())  // neighboring vertex which is outside of roi and not visited previously (i.e. need an id)
      {
        ros_id(vt) = next_id++;
        push_vector.push_back(vt);
      }
    }
  }

  /// Find region of solution, including roi and hard constraints, which is the 1-ring vertices out roi
  /// Contains four parts:
  ///  - if there is any vertex which is no longer roi, set its position back to original position
  ///  - assign a ros id to vertices inside ros + ros-boundary
  ///  - reinitialize rotation matrices, if a vertex is previously ros, use its previous matrix, otherwise set zero
  ///  - reinitialize original, and solution,
  ///      + if a vertex is previously roi, then use its original position in old_origional, else use point().
  ///        In both case we are using "original position" of the vertex.
  ///      + same for solution (it is required to prevent jumping effects)
  void region_of_solution()
  {
    if(!need_preprocess_region_of_solution) { return; }
    need_preprocess_region_of_solution = false;

    std::vector<std::size_t>  old_ros_id_map = ros_id_map;
    std::vector<CR_matrix>    old_rot_mtr    = rot_mtr;
    std::vector<Point>        old_solution   = solution;
    std::vector<Point>        old_original   = original;

    // any vertices which are no longer ROI, should be assigned to their original position, so that:
    // IF a vertex is ROI (actually if its ros + boundary) previously (when previous region_of_solution() is called)
    // we have its original position in old_original
    // ELSE
    // we have its original position in vertex->point()
    // (current ros is actually old ros - we did not change it yet)
    for(typename std::vector<vertex_descriptor>::iterator it = ros.begin(); it != ros.end(); ++it)
    {
      if(!is_roi_vertex(*it)) {
        put(vertex_point_map, *it, old_original[ old_ros_id_map[id(*it)] ]);
      }
    }

    ////////////////////////////////////////////////////////////////
    // assign id to vertices inside: roi, boundary of roi (roi + boundary of roi = ros),
    //                               and boundary of ros
    // keep in mind that id order does not have to be compatible with ros order
    ros.clear(); // clear ros
    ros.insert(ros.end(), roi.begin(), roi.end());

    ros_id_map.assign(num_vertices(m_triangle_mesh), (std::numeric_limits<std::size_t>::max)()); // use max as not assigned mark

    for(std::size_t i = 0; i < roi.size(); i++)  // assign id to all roi vertices
    { ros_id(roi[i]) = i; }

    // now assign an id to vertices on boundary of roi
    std::size_t next_ros_index = roi.size();
    for(std::size_t i = 0; i < roi.size(); i++)
    { assign_ros_id_to_one_ring(roi[i], next_ros_index, ros); }

    std::vector<vertex_descriptor> outside_ros;
    // boundary of ros also must have ids because in CR calculation,
    // one-ring neighbor of ROS vertices are reached.
    for(std::size_t i = roi.size(); i < ros.size(); i++)
    { assign_ros_id_to_one_ring(ros[i], next_ros_index, outside_ros); }
    ////////////////////////////////////////////////////////////////

    // initialize the rotation matrices (size: ros)
    rot_mtr.resize(ros.size());
    for(std::size_t i = 0; i < rot_mtr.size(); i++)
    {
      std::size_t v_ros_id = ros_id(ros[i]);
      std::size_t v_id = id(ros[i]);

      // any vertex which is previously ROS has a rotation matrix
      // use that matrix to prevent jumping effects
      if(old_ros_id_map[v_id] != (std::numeric_limits<std::size_t>::max)()
          && old_ros_id_map[v_id] < old_rot_mtr.size()) {
          // && boundary of ros vertices also have ids so check whether it is ros
        rot_mtr[v_ros_id] = old_rot_mtr[ old_ros_id_map[v_id] ];
      }
      else {
        rot_mtr[v_ros_id] = Closest_rotation_traits().identity_matrix();
      }
    }

    // initialize solution and original (size: ros + boundary_of_ros)

    // for simplifying coding effort, I also put boundary of ros into solution and original
    // because boundary of ros vertices are reached in optimal_rotations() and energy()
    solution.resize(ros.size() + outside_ros.size());
    original.resize(ros.size() + outside_ros.size());

    for(std::size_t i = 0; i < ros.size(); i++)
    {
      std::size_t v_ros_id = ros_id(ros[i]);
      std::size_t v_id = id(ros[i]);

      if(is_roi_vertex(ros[i]) && old_ros_id_map[v_id] != (std::numeric_limits<std::size_t>::max)()) {
        // if it is currently roi and previously ros + boundary
        // (actually I just need to assign old's to new's if a vertex is currently and previously ROI
        // but no harm on assigning if its also previously ros + boundary because
        // those(old_original & old_solution) will be equal to original position)
        original[v_ros_id] = old_original[old_ros_id_map[v_id]];
        solution[v_ros_id] = old_solution[old_ros_id_map[v_id]];
      }
      else {
        solution[v_ros_id] = get(vertex_point_map, ros[i]);
        original[v_ros_id] = get(vertex_point_map, ros[i]);
      }
    }

    for(std::size_t i = 0; i < outside_ros.size(); ++i)
    {
      std::size_t v_ros_id = ros_id(outside_ros[i]);
      original[v_ros_id] = get(vertex_point_map, outside_ros[i]);
      solution[v_ros_id] = get(vertex_point_map, outside_ros[i]);
    }

#ifdef CGAL_DEFORM_MESH_USE_EXPERIMENTAL_SCALE
    scales.resize(ros.size());
    std::fill(scales.begin(), scales.end(), 1.0);
#endif
  }

  /// Assemble Laplacian matrix A of linear system A*X=B
  void assemble_laplacian_and_factorize()
  {
    if(TAG == SPOKES_AND_RIMS)
    {
      assemble_laplacian_and_factorize_spokes_and_rims();
    }
    else
    {
      assemble_laplacian_and_factorize_arap();
    }
  }
  /// Construct matrix that corresponds to left-hand side of eq:lap_ber in user manual
  /// Also constraints are integrated as eq:lap_energy_system in user manual
  void assemble_laplacian_and_factorize_arap()
  {
    if(!need_preprocess_factorization) { return; }
    need_preprocess_factorization = false;

    typename Sparse_linear_solver::Matrix A(ros.size());

    /// assign cotangent Laplacian to ros vertices
    for(std::size_t k = 0; k < ros.size(); k++)
    {
      vertex_descriptor vi = ros[k];
      std::size_t vi_id = ros_id(vi);
      if ( is_roi_vertex(vi) && !is_control_vertex(vi) )          // vertices of ( roi - ctrl )
      {
        double diagonal = 0;
        in_edge_iterator e, e_end;
        for (cpp11::tie(e,e_end) = in_edges(vi, m_triangle_mesh); e != e_end; e++)
        {
          halfedge_descriptor he = halfedge(*e, m_triangle_mesh);
          vertex_descriptor vj = source(he, m_triangle_mesh);
          double wij = hedge_weight[id(he)];  // edge(pi - pj)
          double wji = hedge_weight[id(opposite(he, m_triangle_mesh))]; // edge(pi - pj)
          double total_weight = wij + wji;

          A.set_coef(vi_id, ros_id(vj), -total_weight, true); // off-diagonal coefficient
          diagonal += total_weight;
        }
        // diagonal coefficient
        A.set_coef(vi_id, vi_id, diagonal, true);
      }
      else
      {
        A.set_coef(vi_id, vi_id, 1.0, true);
      }
    }

    // now factorize
    double D;
    last_preprocess_successful = m_solver.factor(A, D);
    CGAL_warning(last_preprocess_successful);
  }
  /// Construct matrix that corresponds to left-hand side of eq:lap_ber_rims in user manual
  /// Also constraints are integrated as eq:lap_energy_system in user manual
  void assemble_laplacian_and_factorize_spokes_and_rims()
  {
    if(!need_preprocess_factorization) { return; }
    need_preprocess_factorization = false;

    typename Sparse_linear_solver::Matrix A(ros.size());

    /// assign cotangent Laplacian to ros vertices
    for(std::size_t k = 0; k < ros.size(); k++)
    {
      vertex_descriptor vi = ros[k];
      std::size_t vi_id = ros_id(vi);
      if ( is_roi_vertex(vi) && !is_control_vertex(vi) ) // vertices of ( roi - ctrl ): free vertices
      {
        double diagonal = 0;
        out_edge_iterator e, e_end;
        for (cpp11::tie(e,e_end) = out_edges(vi, m_triangle_mesh); e != e_end; e++)
        {
          halfedge_descriptor he = halfedge(*e, m_triangle_mesh);
          double total_weight = 0;
          // an edge contribute to energy only if it is part of an incident triangle
          // (i.e it should not be a border edge)
          if(!is_border(he, m_triangle_mesh))
          {
            double wji = hedge_weight[id(he)]; // edge(pj - pi)
            total_weight += wji;
          }

          halfedge_descriptor opp = opposite(he, m_triangle_mesh);
          if(!is_border(opp, m_triangle_mesh))
          {
            double wij = hedge_weight[id(opp)]; // edge(pi - pj)
            total_weight += wij;
          }

          // place coefficient to matrix
          vertex_descriptor vj = target(he, m_triangle_mesh);
          A.set_coef(vi_id, ros_id(vj), -total_weight, true); // off-diagonal coefficient
          diagonal += total_weight;
        }
        // diagonal coefficient
        A.set_coef(vi_id, vi_id, diagonal, true);
      }
      else // constrained vertex
      {
        A.set_coef(vi_id, vi_id, 1.0, true);
      }
    }

    // now factorize
    double D;
    last_preprocess_successful = m_solver.factor(A, D);
    CGAL_warning(last_preprocess_successful);
  }

  /// Local step of iterations, computing optimal rotation matrices
  void optimal_rotations()
  {
    if(TAG == SPOKES_AND_RIMS)
    {
      optimal_rotations_spokes_and_rims();
    }
    else
    {
      optimal_rotations_arap();
    }
  }
  void optimal_rotations_arap()
  {
    Closest_rotation_traits cr_traits;
    CR_matrix cov = cr_traits.zero_matrix();

    // only accumulate ros vertices
    for ( std::size_t k = 0; k < ros.size(); k++ )
    {
      vertex_descriptor vi = ros[k];
      std::size_t vi_id = ros_id(vi);
      // compute covariance matrix (user manual eq:cov_matrix)
      cov = cr_traits.zero_matrix();

      in_edge_iterator e, e_end;

      arap_visitor.rotation_matrix_pre(vi, m_triangle_mesh);

      for (cpp11::tie(e,e_end) = in_edges(vi, m_triangle_mesh); e != e_end; e++)
      {
        halfedge_descriptor he=halfedge(*e, m_triangle_mesh);
        vertex_descriptor vj = source(he, m_triangle_mesh);
        std::size_t vj_id = ros_id(vj);

        const CR_vector& pij = sub_to_CR_vector(original[vi_id], original[vj_id]);
        const CR_vector& qij = sub_to_CR_vector(solution[vi_id], solution[vj_id]);
        double wij = hedge_weight[id(he)];

        cr_traits.add_scalar_t_vector_t_vector_transpose(cov, wij, pij, qij); // cov += wij * (pij * qij)

        if ( vj_id < rot_mtr.size() )
          arap_visitor.update_covariance_matrix(cov, rot_mtr[vj_id]);
      }

      cr_traits.compute_close_rotation(cov, rot_mtr[vi_id]);
    }
  }
  void optimal_rotations_spokes_and_rims()
  {
    Closest_rotation_traits cr_traits;
    CR_matrix cov =cr_traits.zero_matrix();

    // only accumulate ros vertices
    for ( std::size_t k = 0; k < ros.size(); k++ )
    {
      vertex_descriptor vi = ros[k];
      std::size_t vi_id = ros_id(vi);
      // compute covariance matrix
      cov = cr_traits.zero_matrix();

      //iterate through all triangles
      out_edge_iterator e, e_end;
      for (cpp11::tie(e,e_end) = out_edges(vi, m_triangle_mesh); e != e_end; e++)
      {
        halfedge_descriptor he = halfedge(*e, m_triangle_mesh);
        if(is_border(he, m_triangle_mesh)) { continue; } // no facet
        // iterate edges around facet
        halfedge_descriptor hedge_around_facet = he;
        do
        {
          vertex_descriptor v1 = target(hedge_around_facet, m_triangle_mesh);
          vertex_descriptor v2 = source(hedge_around_facet, m_triangle_mesh);

          std::size_t v1_id = ros_id(v1); std::size_t v2_id = ros_id(v2);

          const CR_vector& p12 = sub_to_CR_vector(original[v1_id], original[v2_id]);
          const CR_vector& q12 = sub_to_CR_vector(solution[v1_id], solution[v2_id]);
          double w12 = hedge_weight[id(hedge_around_facet)];

          cr_traits.add_scalar_t_vector_t_vector_transpose(cov, w12, p12, q12); // cov += w12 * (p12 * q12);

        } while( (hedge_around_facet = next(hedge_around_facet, m_triangle_mesh)) != he);
      }

      cr_traits.compute_close_rotation(cov, rot_mtr[vi_id]);
    }
  }

#ifdef CGAL_DEFORM_MESH_USE_EXPERIMENTAL_SCALE
  void optimal_scales()
  {
    for ( std::size_t k = 0; k < ros.size(); k++ )
    {
      vertex_descriptor vi = ros[k];
      std::size_t vi_id = ros_id(vi);
      // compute covariance matrix (user manual eq:cov_matrix)
      double eT_eR = 0, eRT_eR = 0;

      in_edge_iterator e, e_end;
      for (cpp11::tie(e,e_end) = in_edges(vi, m_triangle_mesh); e != e_end; e++)
      {
        halfedge_descriptor he = *e;
        vertex_descriptor vj = source(he, m_triangle_mesh);
        std::size_t vj_id = ros_id(vj);

        const CR_vector& pij = sub_to_CR_vector(original[vi_id], original[vj_id]);
        const CR_vector& qij = sub_to_CR_vector(solution[vi_id], solution[vj_id]);

        double wij = hedge_weight[id(he)];

        const CR_vector& pRij = rot_mtr[vi_id] * pij;
        eRT_eR += pRij[0]*pRij[0] + pRij[1]*pRij[1] + pRij[2]*pRij[2];
        eT_eR  += qij[0]*pRij[0]  + qij[1]*pRij[1]  + qij[2]*pRij[2];
      }

      scales[vi_id] = eT_eR / eRT_eR;
    }
  }
#endif

  /// Global step of iterations, updating solution
  void update_solution()
  {
    if(TAG == SPOKES_AND_RIMS)
    {
      update_solution_spokes_and_rims();
    }
    else
    {
      update_solution_arap();
    }
  }
  /// calculate right-hand side of eq:lap_ber in user manual and solve the system
  void update_solution_arap()
  {
    typename Sparse_linear_solver::Vector X(ros.size()), Bx(ros.size());
    typename Sparse_linear_solver::Vector Y(ros.size()), By(ros.size());
    typename Sparse_linear_solver::Vector Z(ros.size()), Bz(ros.size());

    Closest_rotation_traits cr_traits;

    // assemble right columns of linear system
    for ( std::size_t k = 0; k < ros.size(); k++ )
    {
      vertex_descriptor vi = ros[k];
      std::size_t vi_id = ros_id(vi);

      if ( is_roi_vertex(vi) && !is_control_vertex(vi) )
      {// free vertices
        // sum of right-hand side of eq:lap_ber in user manual
        CR_vector xyz = cr_traits.vector(0, 0, 0);

        in_edge_iterator e, e_end;
        for (cpp11::tie(e,e_end) = in_edges(vi, m_triangle_mesh); e != e_end; e++)
        {
          halfedge_descriptor he = halfedge(*e, m_triangle_mesh);
          vertex_descriptor vj = source(he, m_triangle_mesh);
          std::size_t vj_id = ros_id(vj);

          const CR_vector& pij = sub_to_CR_vector(original[vi_id], original[vj_id]);

          double wij = hedge_weight[id(he)];
          double wji = hedge_weight[id(opposite(he, m_triangle_mesh))];
#ifndef CGAL_DEFORM_MESH_USE_EXPERIMENTAL_SCALE
          cr_traits.add__scalar_t_matrix_p_scalar_t_matrix__t_vector(xyz, wij, rot_mtr[vi_id], wji, rot_mtr[vj_id], pij);
#else
        cr_traits.add__scalar_t_matrix_p_scalar_t_matrix__t_vector(xyz, wij * scales[vi_id], rot_mtr[vi_id],
          wji * scales[vj_id], rot_mtr[vj_id], pij);
#endif
          // corresponds xyz += (wij*rot_mtr[vi_id] + wji*rot_mtr[vj_id]) * pij
        }
        Bx[vi_id] = cr_traits.vector_coordinate(xyz, 0);
        By[vi_id] = cr_traits.vector_coordinate(xyz, 1);
        Bz[vi_id] = cr_traits.vector_coordinate(xyz, 2);
      }
      else
      {// constrained vertex
        Bx[vi_id] = solution[vi_id][0]; By[vi_id] = solution[vi_id][1]; Bz[vi_id] = solution[vi_id][2];
      }
    }

    // solve "A*X = B".
    bool is_all_solved = m_solver.linear_solver(Bx, X) && m_solver.linear_solver(By, Y) && m_solver.linear_solver(Bz, Z);
    if(!is_all_solved) {
      CGAL_warning(false);
      return;
    }
    // copy to solution
    for (std::size_t i = 0; i < ros.size(); i++)
    {
      std::size_t v_id = ros_id(ros[i]);
      Point p(X[v_id], Y[v_id], Z[v_id]);
      if (!is_ctrl_map[id(ros[i])])
        solution[v_id] = p;
    }
  }
  /// calculate right-hand side of eq:lap_ber_rims in user manual and solve the system
  void update_solution_spokes_and_rims()
  {
    typename Sparse_linear_solver::Vector X(ros.size()), Bx(ros.size());
    typename Sparse_linear_solver::Vector Y(ros.size()), By(ros.size());
    typename Sparse_linear_solver::Vector Z(ros.size()), Bz(ros.size());

    Closest_rotation_traits cr_traits;

    // assemble right columns of linear system
    for ( std::size_t k = 0; k < ros.size(); k++ )
    {
      vertex_descriptor vi = ros[k];
      std::size_t vi_id = ros_id(vi);

      if ( is_roi_vertex(vi) && !is_control_vertex(vi) )
      {// free vertices
        // sum of right-hand side of eq:lap_ber_rims in user manual
        CR_vector xyz = cr_traits.vector(0, 0, 0);

        out_edge_iterator e, e_end;
        for (cpp11::tie(e,e_end) = out_edges(vi, m_triangle_mesh); e != e_end; e++)
        {
          halfedge_descriptor he = halfedge(*e, m_triangle_mesh);
          vertex_descriptor vj = target(he, m_triangle_mesh);
          std::size_t vj_id = ros_id(vj);

          const CR_vector& pij = sub_to_CR_vector(original[vi_id], original[vj_id]);

          if(!is_border(he, m_triangle_mesh))
          {
            vertex_descriptor vn = target(next(he, m_triangle_mesh), m_triangle_mesh); // opp vertex of e_ij
            double wji = hedge_weight[id(he)] / 3.0;  // edge(pj - pi)
            cr_traits.add_scalar_t_matrix_sum_t_vector(xyz, wji, rot_mtr[vi_id], rot_mtr[vj_id], rot_mtr[ros_id(vn)], pij);
            // corresponds  xyz += wji*(rot_mtr[vi_id] + rot_mtr[vj_id] + rot_mtr[ros_id(vn)])*pij;
          }

          halfedge_descriptor opp = opposite(he, m_triangle_mesh);
          if(!is_border(opp, m_triangle_mesh))
          {
            vertex_descriptor vm = target(next(opp, m_triangle_mesh), m_triangle_mesh); // other opp vertex of e_ij
            double wij = hedge_weight[id(opp)] / 3.0;  // edge(pi - pj)
            cr_traits.add_scalar_t_matrix_sum_t_vector(xyz, wij, rot_mtr[vi_id], rot_mtr[vj_id], rot_mtr[ros_id(vm)], pij);
            // corresponds xyz += wij * ( rot_mtr[vi_id] + rot_mtr[vj_id] + rot_mtr[ros_id(vm)] ) * pij
          }
        }
        Bx[vi_id] = cr_traits.vector_coordinate(xyz, 0);
        By[vi_id] = cr_traits.vector_coordinate(xyz, 1);
        Bz[vi_id] = cr_traits.vector_coordinate(xyz, 2);
      }
      else
      {// constrained vertices
        Bx[vi_id] = solution[vi_id][0]; By[vi_id] = solution[vi_id][1]; Bz[vi_id] = solution[vi_id][2];
      }
    }
    // solve "A*X = B".
    bool is_all_solved = m_solver.linear_solver(Bx, X) && m_solver.linear_solver(By, Y) && m_solver.linear_solver(Bz, Z);
    if(!is_all_solved) {
      CGAL_warning(false);
      return;
    }

    // copy to solution
    for (std::size_t i = 0; i < ros.size(); i++)
    {
      std::size_t v_id = ros_id(ros[i]);
      Point p(X[v_id], Y[v_id], Z[v_id]);
      if (!is_ctrl_map[id(ros[i])])
        solution[v_id] = p;
    }
  }

  /// Assign solution to target surface mesh
  void assign_solution()
  {
    for(std::size_t i = 0; i < ros.size(); ++i){
      std::size_t v_id = ros_id(ros[i]);
      if(is_roi_vertex(ros[i]))
      {
        put(vertex_point_map, ros[i], solution[v_id]);
      }
    }
  }

  /// Compute modeling energy
  double energy() const
  {
    if(TAG == SPOKES_AND_RIMS)
    {
      return energy_spokes_and_rims();
    }
    else
    {
      return energy_arap();
      return 0;
    }
  }
  double energy_arap() const
  {
    Closest_rotation_traits cr_traits;

    double sum_of_energy = 0;
    // only accumulate ros vertices
    for( std::size_t k = 0; k < ros.size(); k++ )
    {
      vertex_descriptor vi = ros[k];
      std::size_t vi_id = ros_id(vi);

      in_edge_iterator e, e_end;
      for (cpp11::tie(e,e_end) = in_edges(vi, m_triangle_mesh); e != e_end; e++)
      {
        halfedge_descriptor he = halfedge(*e, m_triangle_mesh);
        vertex_descriptor vj = source(he, m_triangle_mesh);
        std::size_t vj_id = ros_id(vj);

        const CR_vector& pij = sub_to_CR_vector(original[vi_id], original[vj_id]);
        const CR_vector& qij = sub_to_CR_vector(solution[vi_id], solution[vj_id]);

        double wij = hedge_weight[id(he)];

        sum_of_energy += wij * cr_traits.squared_norm_vector_scalar_vector_subs(qij, rot_mtr[vi_id], pij);
        // sum_of_energy += wij * ( qij - rot_mtr[vi_id]*pij )^2
      }
    }
    return sum_of_energy;
  }
  double energy_spokes_and_rims() const
  {
    Closest_rotation_traits cr_traits;

    double sum_of_energy = 0;
    // only accumulate ros vertices
    for( std::size_t k = 0; k < ros.size(); k++ )
    {
      vertex_descriptor vi = ros[k];
      std::size_t vi_id = ros_id(vi);
      //iterate through all triangles
      out_edge_iterator e, e_end;
      for (cpp11::tie(e,e_end) = out_edges(vi, m_triangle_mesh); e != e_end; e++)
      {
        halfedge_descriptor he = halfedge(*e, m_triangle_mesh);
        if(is_border(he, m_triangle_mesh)) { continue; } // no facet
        // iterate edges around facet
        halfedge_descriptor hedge_around_facet = he;
        do
        {
          vertex_descriptor v1 = target(hedge_around_facet, m_triangle_mesh);
          vertex_descriptor v2 = source(hedge_around_facet, m_triangle_mesh);
          std::size_t v1_id = ros_id(v1); std::size_t v2_id = ros_id(v2);

          const CR_vector& p12 = sub_to_CR_vector(original[v1_id], original[v2_id]);
          const CR_vector& q12 = sub_to_CR_vector(solution[v1_id], solution[v2_id]);
          double w12 = hedge_weight[id(hedge_around_facet)];

          sum_of_energy += w12 * cr_traits.squared_norm_vector_scalar_vector_subs(q12, rot_mtr[vi_id], p12);
          // sum_of_energy += w12 * ( q12 - rot_mtr[vi_id]*p12 )^2

        } while( (hedge_around_facet = next(hedge_around_facet, m_triangle_mesh)) != he);
      }
    }
    return sum_of_energy;
  }

  void need_preprocess_both()
  {
    need_preprocess_factorization = true;
    need_preprocess_region_of_solution = true;
  }

  /// p1 - p2, return CR_vector
  CR_vector sub_to_CR_vector(const Point& p1, const Point& p2) const
  {
    return Closest_rotation_traits().vector(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
  }

  template<class Vect>
  Point add_to_point(const Point& p, const Vect& v) {
    return Point(p[0] + v[0], p[1] + v[1], p[2] + v[2]);
  }

  template<class Vect>
  Vect sub_to_vect(const Point& p, const Vect& v) {
    return Vect(p[0] - v[0], p[1] - v[1], p[2] - v[2]);
  }

  /// shorthand of get(vertex_index_map, v)
  std::size_t id(vertex_descriptor vd) const
  { return get(vertex_index_map, vd); }

  std::size_t& ros_id(vertex_descriptor vd)
  { return ros_id_map[id(vd)]; }
  std::size_t  ros_id(vertex_descriptor vd) const
  { return ros_id_map[id(vd)]; }

  /// shorthand of get(hedge_index_map, e)
  std::size_t id(halfedge_descriptor e) const
  {
    return get(hedge_index_map, e);
  }
};
} //namespace CGAL
#endif  // CGAL_SURFACE_MESH_DEFORMATION_H
