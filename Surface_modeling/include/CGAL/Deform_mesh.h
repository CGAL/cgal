// Copyright (c) 2011-2013 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id:$
//
// Author(s)     : Yin Xu, Andreas Fabri and Ilker O. Yaz

#ifndef CGAL_DEFORM_MESH_H
#define CGAL_DEFORM_MESH_H

#include <CGAL/internal/Surface_modeling/Weights.h>
#include <CGAL/Default.h>

#include <vector>
#include <list>
#include <utility>
#include <limits>

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

/// \ingroup PkgSurfaceModeling
///@brief Deformation algorithm type
enum Deformation_algorithm_tag
{ 
  ORIGINAL_ARAP,  /**< use original as-rigid-as possible algorithm */
  SPOKES_AND_RIMS /**< use spokes and rims version of as-rigid-as possible algorithm */
};

/// @cond CGAL_DOCUMENT_INTERNAL
namespace internal {
template<class HalfedgeGraph, Deformation_algorithm_tag deformation_algorithm_tag>
struct Weight_calculator_selector {
  typedef Uniform_weight<HalfedgeGraph> weight_calculator;
};

template<class HalfedgeGraph>
struct Weight_calculator_selector<HalfedgeGraph, CGAL::SPOKES_AND_RIMS> {
  typedef Single_cotangent_weight<HalfedgeGraph> weight_calculator;
};

template<class HalfedgeGraph>
struct Weight_calculator_selector<HalfedgeGraph, CGAL::ORIGINAL_ARAP> {
  typedef Cotangent_weight<HalfedgeGraph> weight_calculator;
};
}//namespace internal
/// @endcond

 ///
 /// \ingroup PkgSurfaceModeling
 /// @brief Class providing the functionalities for deforming a triangulated surface mesh
 ///
 /// @tparam HG a model of HalfedgeGraph 
 /// @tparam VIM a model of `ReadOnlyPropertyMap`</a>  with Deform_mesh::vertex_descriptor as key and `unsigned int` as value type,
 ///         containing unique indices to vertices with offset 0
 /// @tparam EIM a model of `ReadOnlyPropertyMap`</a>  with Deform_mesh::edge_descriptor as key and `unsigned int` as value type
 ///         containing unique indices to vertices with offset 0
 /// @tparam TAG tag for selecting the deformation algorithm
 /// @tparam WC a model of SurfaceModelingWeightCalculator, with `WC::Halfedge_graph` being `HG`
 /// @tparam ST a model of SparseLinearAlgebraTraitsWithPreFactor_d. If \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available 
 /// and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits` is provided as default parameter.\n
  /// \code
 ///     CGAL::Eigen_solver_traits<
 ///         Eigen::SparseLU<
 ///            CGAL::Eigen_sparse_matrix<double>::EigenType,
 ///            Eigen::COLAMDOrdering<int> >  >
 /// \endcode
 /// @tparam CR a model of DeformationClosestRotationTraits_3. If \ref thirdpartyEigen "Eigen" 3.1 (or greater) is available and `CGAL_EIGEN3_ENABLED` is defined, 
 /// `Deformation_Eigen_polar_closest_rotation_traits_3` is provided as default parameter.
 /// @tparam VPM a model of `ReadWritePropertyMap`</a>  with Deform_mesh::vertex_descriptor as key and a point from a \cgal Kernel as value type.
template <
  class HG, 
  class VIM, 
  class EIM,
  Deformation_algorithm_tag TAG = SPOKES_AND_RIMS,
  class WC = Default,
  class ST = Default, 
  class CR = Default,
  class VPM = Default
  >
class Deform_mesh
{
//Typedefs
public:

  /// \name Template parameters
  /// @{
  // typedefed template parameters, main reason is doxygen creates autolink to typedefs but not template parameters
  typedef HG Halfedge_graph; /**< model of HalfedgeGraph */  
  typedef VIM Vertex_index_map; /**< model of `ReadWritePropertyMap`  with Deform_mesh::vertex_descriptor as key and `unsigned int` as value type */
  typedef EIM Edge_index_map; /**< model of `ReadWritePropertyMap`</a>  with Deform_mesh::edge_descriptor as key and `unsigned int` as value type */

// weight calculator
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    WC,
    typename internal::Weight_calculator_selector<HG, TAG>::weight_calculator
  >::type Weight_calculator;
#else
  typedef WC Weight_calculator; /**< model of SurfaceModelingWeightCalculator */
#endif

// sparse linear solver
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    ST,
  #if defined(CGAL_EIGEN3_ENABLED)
    #if defined(CGAL_SUPERLU_ENABLED)
      CGAL::Eigen_solver_traits<Eigen::SuperLU<CGAL::Eigen_sparse_matrix<double>::EigenType> >
    #else
      CGAL::Eigen_solver_traits<
          Eigen::SparseLU<
            CGAL::Eigen_sparse_matrix<double>::EigenType,
            Eigen::COLAMDOrdering<int> >  >
    #endif
  #else
    ST // no parameter provided, and Eigen is not enabled: so don't compile!
  #endif
  >::type Sparse_linear_solver;
#else
  typedef ST Sparse_linear_solver; /**< model of SparseLinearAlgebraTraitsWithPreFactor_d */
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
  typedef CR Closest_rotation_traits; /**< model of DeformationClosestRotationTraits_3 */
#endif

// vertex point pmap
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    VPM,
    typename boost::property_map<Halfedge_graph, CGAL::vertex_point_t>::type
  >::type Vertex_point_map;
#else
  typedef VPM Vertex_point_map; /**<  a model of `ReadWritePropertyMap`</a>  with Deform_mesh::vertex_descriptor as key and `Point` as value type */
#endif
  /// @}

/// \name Public Types
/// @{
  /// The type for vertex representative objects
  typedef typename boost::graph_traits<Halfedge_graph>::vertex_descriptor	vertex_descriptor;
  /// The type for edge representative objects
  typedef typename boost::graph_traits<Halfedge_graph>::edge_descriptor		edge_descriptor;
  /// The 3D point type
  typedef typename boost::property_traits<Vertex_point_map>::value_type Point;
  /// %Iterator over vertices in the region-of-interest.
  typedef typename std::vector<vertex_descriptor>::const_iterator  Roi_vertex_const_iterator;
/// @}

private:
  typedef Deform_mesh<HG, VIM, EIM, TAG, WC, ST, CR> Self;
  // Repeat Halfedge_graph types
  typedef typename boost::graph_traits<Halfedge_graph>::vertex_iterator     vertex_iterator;
  typedef typename boost::graph_traits<Halfedge_graph>::edge_iterator       edge_iterator;
  typedef typename boost::graph_traits<Halfedge_graph>::in_edge_iterator    in_edge_iterator;
  typedef typename boost::graph_traits<Halfedge_graph>::out_edge_iterator   out_edge_iterator;

  typedef typename Closest_rotation_traits::Matrix CR_matrix;
  typedef typename Closest_rotation_traits::Vector CR_vector;

public:

// Data members.
  Halfedge_graph& m_halfedge_graph;															/**< Source triangulated surface mesh for modeling */

  std::vector<Point> original;                        ///< original positions of roi (size: ros + boundary_of_ros)
  std::vector<Point> solution;                        ///< storing position of ros vertices during iterations (size: ros + boundary_of_ros)

  Vertex_index_map vertex_index_map;                  ///< storing indices of all vertices
  Edge_index_map   edge_index_map;                    ///< storing indices of all edges

  std::vector<vertex_descriptor> roi;                 ///< region of interest
  std::vector<vertex_descriptor> ros;                 ///< region of solution, including roi and hard constraints on boundary of roi

  std::vector<std::size_t> ros_id_map;                ///< (size: num vertices)
  std::vector<bool>        is_roi_map;                ///< (size: num vertices)
  std::vector<bool>        is_hdl_map;                ///< (size: num vertices)

  std::vector<double> edge_weight;                    ///< all edge weights 
  std::vector<CR_matrix> rot_mtr;                     ///< rotation matrices of ros vertices (size: ros)

  Sparse_linear_solver m_solver;                      ///< linear sparse solver
  unsigned int iterations;                            ///< number of maximal iterations
  double tolerance;                                   ///< tolerance of convergence 

  bool need_preprocess_factorization;                 ///< is there any need to compute L and factorize
  bool need_preprocess_region_of_solution;            ///< is there any need to compute region of solution 

  bool last_preprocess_successful;                    ///< stores the result of last call to preprocess()

  Weight_calculator weight_calculator;

  Vertex_point_map vertex_point_map;

#ifdef CGAL_DEFORM_MESH_USE_EXPERIMENTAL_SCALE
  std::vector<double> scales;
#endif
private:
  Deform_mesh(const Self& s) { } 

// Public methods
public:

  /// \cond SKIP_FROM_MANUAL
  Deform_mesh(Halfedge_graph& halfedge_graph, 
              Vertex_index_map vertex_index_map, 
              Edge_index_map edge_index_map,
              unsigned int iterations = 5,
              double tolerance = 1e-4,
              Weight_calculator weight_calculator = Weight_calculator()
              )
    : m_halfedge_graph(halfedge_graph), vertex_index_map(vertex_index_map), edge_index_map(edge_index_map),
      ros_id_map(std::vector<std::size_t>(boost::num_vertices(halfedge_graph), (std::numeric_limits<std::size_t>::max)() )),
      is_roi_map(std::vector<bool>(boost::num_vertices(halfedge_graph), false)),
      is_hdl_map(std::vector<bool>(boost::num_vertices(halfedge_graph), false)),
      iterations(iterations), tolerance(tolerance),
      need_preprocess_factorization(true), 
      need_preprocess_region_of_solution(true),
      last_preprocess_successful(false),
      weight_calculator(weight_calculator),
      vertex_point_map(boost::get(vertex_point, halfedge_graph))
  {
    init();
  }
  /// \endcond
/// \name Construction
/// @{
    /**
   * The constructor of a deformation object
   *
   * @pre the halfedge_graph consists of only triangular facets
   * @param halfedge_graph triangulated surface mesh used to deform
   * @param vertex_index_map property map for associating an id to each vertex
   * @param edge_index_map property map for associating an id to each edge
   * @param vertex_point_map property map used to access the points associated to each vertex of the graph.
   *        It is default to `boost::get(vertex_point, halfedge_graph)` and can be omitted.
   * @param iterations see `set_iterations()` for more details
   * @param tolerance  see `set_tolerance()` for more details
   * @param weight_calculator function object or pointer for weight calculation
   */
  Deform_mesh(Halfedge_graph& halfedge_graph, 
    Vertex_index_map vertex_index_map, 
    Edge_index_map edge_index_map,
    Vertex_point_map vertex_point_map,
    unsigned int iterations = 5,
    double tolerance = 1e-4,
    Weight_calculator weight_calculator = Weight_calculator()
    )
    : m_halfedge_graph(halfedge_graph), vertex_index_map(vertex_index_map), edge_index_map(edge_index_map),
    ros_id_map(std::vector<std::size_t>(boost::num_vertices(halfedge_graph), (std::numeric_limits<std::size_t>::max)() )),
    is_roi_map(std::vector<bool>(boost::num_vertices(halfedge_graph), false)),
    is_hdl_map(std::vector<bool>(boost::num_vertices(halfedge_graph), false)),
    iterations(iterations), tolerance(tolerance),
    need_preprocess_factorization(true), 
    need_preprocess_region_of_solution(true),
    last_preprocess_successful(false),
    weight_calculator(weight_calculator),
    vertex_point_map(vertex_point_map)
  {
    init();
  }
/// @}

private:
  void init() {
    // compute edge weights
    edge_iterator eb, ee;
    edge_weight.reserve(boost::num_edges(m_halfedge_graph));
    for(boost::tie(eb, ee) = boost::edges(m_halfedge_graph); eb != ee; ++eb)
    {
      edge_weight.push_back(
        this->weight_calculator(*eb, m_halfedge_graph, vertex_point_map));
    }
  }

public:

/// \name Preprocessing
/// @{
  /**
   * Puts `*this` in the same state as after the creation (except iterations and tolerance).
   */
  void reset()
  {
    need_preprocess_both();
    // clear vertices
    roi.clear();
    is_roi_map.assign(boost::num_vertices(m_halfedge_graph), false);
    is_hdl_map.assign(boost::num_vertices(m_halfedge_graph), false);
  }
  
  /**
   * Inserts a vertex in the set of control vertices. The vertex is also inserted in the region-of-interest if it is not already in it.
   * @param vd the vertex to be inserted
   * @return `true` if the insertion is successful
   */
  bool insert_control_vertex(vertex_descriptor vd)
  {
    if(is_control(vd)) { return false; }
    need_preprocess_both();

    insert_roi_vertex(vd); // also insert it as roi

    is_hdl_map[id(vd)] = true;
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
   * Erases a vertex from control vertices.
   * @param vd the vertex to be erased
   * @return `true` if the removal is successful
   */
  bool erase_control(vertex_descriptor vd)
  {
    if(!is_control(vd)) { return false; }
    
    need_preprocess_both();
    is_hdl_map[id(vd)] = false;
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
   * @return `true` if the insertion is successful
   */
  bool insert_roi_vertex(vertex_descriptor vd)   
  {
    if(is_roi(vd)) { return false; }
    need_preprocess_both();

    is_roi_map[id(vd)] = true;
    roi.push_back(vd);
    return true;
  }

  /**
   * Erases a vertex from the region-of-interest. The vertex is also removed from control vertices if possible.
   * \note The next call to `preprocess()`, any vertex which is no longer in the region-of-interest will be assigned to its original position 
   * (that is position of the vertex at the time of construction or after the last call to `overwrite_original_positions()`).
   * @param vd the vertex to be erased
   * @return `true` if the removal is successful
   */
  bool erase_roi(vertex_descriptor vd)   
  {
    if(!is_roi(vd)) { return false; }  
    
    erase_control(vd); // also erase from being control

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
   * Returns the range of vertices in the region-of-interest.
   */
  std::pair<Roi_vertex_const_iterator, Roi_vertex_const_iterator> roi_vertices() const
  {
    return std::make_pair(roi.begin(), roi.end());
  }

  /**
   * Assembles and factorizes the Laplacian matrix used in the function `deform()`.
   * \note A modification of the set of control vertices or the region-of-interest invalidates the
   * preprocessing data.
   * @return `true` if the Laplacian matrix factorization is successful.
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
   * @param vd the control vertex to be assigned target position
   * @param target_position the new target position
   */
  void assign(vertex_descriptor vd, const Point& target_position)
  {
    region_of_solution(); // we require ros ids, so if there is any need to preprocess of region of solution -do it.

    if(!is_control(vd)) { return; }
    solution[ros_id(vd)] = target_position;
  }

  /**
   * Sets the target position of each control vertex in the range `[begin,end[` by applying a translation of vertex `t` to its original position
   * (that is its positions at the time of the functor construction or after the last call to `overwrite_original_positions()`).
   * \note A call to this function cancels the last call to `rotate()`, `translate()`, or `assign()`.
   *
   * @tparam InputIterator input iterator type with `vertex_descriptor` as value type
   * @tparam Vect is a 3D vector class, `Vect::operator[](int i)` with i=0,1 or 2 returns its coordinates
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
      solution[v_id] = add_to_point(original[v_id], t);
    }
  }

  /**
   * Sets the target position of each control vertex in the range `[begin,end[` by applying
   * a rotation around `rotation_center` defined by the quaternion `quat`, followed by a
   * translation by vector `t` to its original position (that is its positions at the time
   * of the functor construction or after the last call to `overwrite_original_positions()`).
   * \note A call to this function cancels the last call to `rotate()`, `translate()`, or `assign()`.
   *
   * @tparam InputIterator input iterator type with `vertex_descriptor` as value type
   * @tparam Quaternion is a quaternion class with `Vect operator*(Quaternion, Vect)` being defined and returns the product of a quaternion with a vector
   * @tparam Vect is a 3D vector class, `Vect(double x,double y, double z)` being a constructor available and `Vect::operator[](int i)` with i=0,1 or 2 returns its coordinates
   *
   * @param begin first iterator of the range of vertices
   * @param end past-the-end iterator of the range of vertices
   * @param rotation_center center of rotation
   * @param quat rotation holder quaternion
   * @param t post translation vector
   */
  template <typename InputIterator, typename Quaternion, typename Vect>
  void rotate(InputIterator begin, InputIterator end, const Point& rotation_center, const Quaternion& quat, const Vect& t)
  {
    region_of_solution(); // we require ros ids, so if there is any need to preprocess of region of solution -do it.

    for(; begin != end; ++begin) {
      std::size_t v_id = ros_id(*begin);
      Vect v = quat * sub_to_vector<Vect>(original[v_id], rotation_center);
      const Point& rotated = add_to_point(rotation_center, v);
      solution[v_id] = Point(rotated[0] + t[0], rotated[1] + t[1], rotated[2] + t[2]);
    }
  }

  /**
   * Deforms the region-of-interest according to the deformation algorithm, using the target positions of each control vertex set by using `rotate()`, `translate()`, or `assign()`.
   * The points associated to each vertex of the input graph that are inside the region-of-interest are updated. The initial guess for solving the
   * deformation problem is using the points associated to the input graph before calling the function.
   * \note Nothing happens if `preprocess()` returns `false`.
   * @see set_iterations(unsigned int iterations), set_tolerance(double tolerance), deform(unsigned int iterations, double tolerance)
   */
  void deform()
  {
    deform(iterations, tolerance);
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
    // copy solution to target mesh
    assign_solution();
  }
/// @} Deformation

/// \name Utilities
/// @{

  /**
   * Getter of `set_iterations()`
   */
  unsigned int get_iterations()
  { return iterations; }
  
  /**
   * Getter of `set_tolerance()`
   */
  double get_tolerance()
  { return tolerance; }

  /**
   * Sets the number of iterations used in `deform()`
   */
  void set_iterations(unsigned int iterations)
  { this->iterations = iterations; }
  
   /// @brief Sets the tolerance of convergence used in `deform()`.
   /// Set to zero if energy based termination is not required, which also eliminates energy calculation effort in each iteration. 
   ///
   /// `tolerance >` \f$|\mathrm{energy}(m_i) - \mathrm{energy}(m_{i-1})| / \mathrm{energy}(m_i)\f$ will be used as a termination criterium.
  void set_tolerance(double tolerance)
  { this->tolerance = tolerance; }

  /**
   * Queries whether a vertex is inside the region-of-interest.
   * @param vd the query vertex
   * @return `true` if `vd` has been added (and not removed) to the region-of-interest.
   */
  bool is_roi(vertex_descriptor vd) const
  { return is_roi_map[id(vd)]; }

  /**
   * Queries whether a vertex is a control vertex.
   * @param vd the query vertex
   * @return `true` if `vd` has been added (and not removed) to the set of control vertices.
   */
  bool is_control(vertex_descriptor vd) const
  { return is_hdl_map[id(vd)]; }

  /**
   * Provides access to the halfedge graph being deformed
   * @return the halfedge graph
   */
  const Halfedge_graph& halfedge_graph() const
  { return m_halfedge_graph; }

  /**
   * Sets the original positions to be the current positions for vertices inside region-of-interest. Calling this function has the same effect as creating
   * a new deformation object with the current deformed halfedge-graph, keeping the region-of-interest and control vertices.
   * \note if the region-of-interest or control vertices have been modified since the last call to `preprocess()`,
   * it will be called prior to the overwrite.
   */
  void overwrite_original_positions()
  {
    if(roi.empty()) { return; } // no ROI to overwrite

    region_of_solution(); // the roi should be preprocessed since we are using original_position vec

    Roi_vertex_const_iterator rb, re;
    for(boost::tie(rb, re) = roi_vertices(); rb != re; ++rb)
    {
      original[ros_id(*rb)] = get(vertex_point_map, (*rb));
    }

    // now I need to compute weights for edges incident to roi vertices
    std::vector<bool> is_weight_computed(boost::num_edges(m_halfedge_graph), false);
    for(boost::tie(rb, re) = roi_vertices(); rb != re; ++rb)
    {
      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(*rb, m_halfedge_graph); e != e_end; e++)
      {
        std::size_t id_e = id(*e);
        if(is_weight_computed[id_e]) { continue; }

        edge_weight[id_e] = weight_calculator(*e, m_halfedge_graph, vertex_point_map);
        is_weight_computed[id_e] = true;

        edge_descriptor e_opp = CGAL::opposite_edge(*e, m_halfedge_graph);
        std::size_t id_e_opp = id(e_opp);

        edge_weight[id_e_opp] = weight_calculator(e_opp, m_halfedge_graph, vertex_point_map);
        is_weight_computed[id_e_opp] = true;
      }
    }    

    // also set rotation matrix to identity
    std::fill(rot_mtr.begin(), rot_mtr.end(), Closest_rotation_traits().identity_matrix());

    need_preprocess_both(); // now we need reprocess
  }

/// @} Utilities


private:

  /// Assigns id to one ring neighbor of vd, and also push them into push_vector
  void assign_ros_id_to_one_ring(vertex_descriptor vd, 
                             std::size_t& next_id, 
                             std::vector<vertex_descriptor>& push_vector)
  {
    in_edge_iterator e, e_end;
    for (boost::tie(e,e_end) = boost::in_edges(vd, m_halfedge_graph); e != e_end; e++)
    {
      vertex_descriptor vt = boost::source(*e, m_halfedge_graph);    
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
    std::vector<CR_matrix>   old_rot_mtr    = rot_mtr;
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
      if(!is_roi(*it)) {
        put(vertex_point_map, *it, old_original[ old_ros_id_map[id(*it)] ]);
      }
    }

    ////////////////////////////////////////////////////////////////
    // assign id to vertices inside: roi, boundary of roi (roi + boundary of roi = ros),
    //                               and boundary of ros
    // keep in mind that id order does not have to be compatible with ros order
    ros.clear(); // clear ros    
    ros.insert(ros.end(), roi.begin(), roi.end()); 

    ros_id_map.assign(boost::num_vertices(m_halfedge_graph), (std::numeric_limits<std::size_t>::max)()); // use max as not assigned mark

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

      if(is_roi(ros[i]) && old_ros_id_map[v_id] != (std::numeric_limits<std::size_t>::max)()) { 
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
      if ( is_roi(vi) && !is_control(vi) )          // vertices of ( roi - hdl )
      {
        double diagonal = 0;
        in_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::in_edges(vi, m_halfedge_graph); e != e_end; e++)
        {
          vertex_descriptor vj = boost::source(*e, m_halfedge_graph);
          double wij = edge_weight[id(*e)];  // edge(pi - pj)
          double wji = edge_weight[id(CGAL::opposite_edge(*e, m_halfedge_graph))]; // edge(pi - pj)
          double total_weight = wij + wji;

          A.set_coef(vi_id, ros_id(vj), -total_weight, true);	// off-diagonal coefficient
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
    last_preprocess_successful = m_solver.pre_factor(A, D);
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
      if ( is_roi(vi) && !is_control(vi) ) // vertices of ( roi - hdl ): free vertices
      {
        double diagonal = 0;
        out_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::out_edges(vi, m_halfedge_graph); e != e_end; e++)
        {
          double total_weight = 0;
          // an edge contribute to energy only if it is part of an incident triangle 
          // (i.e it should not be a border edge)
          if(!boost::get(CGAL::edge_is_border, m_halfedge_graph, *e)) 
          {
            double wji = edge_weight[id(*e)]; // edge(pj - pi)
            total_weight += wji; 
          }

          edge_descriptor opp = CGAL::opposite_edge(*e, m_halfedge_graph);
          if(!boost::get(CGAL::edge_is_border, m_halfedge_graph, opp))
          {
            double wij = edge_weight[id(opp)]; // edge(pi - pj)
            total_weight += wij;
          }

          // place coefficient to matrix
          vertex_descriptor vj = boost::target(*e, m_halfedge_graph);
          A.set_coef(vi_id, ros_id(vj), -total_weight, true);	// off-diagonal coefficient
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
    last_preprocess_successful = m_solver.pre_factor(A, D);
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
      for (boost::tie(e,e_end) = boost::in_edges(vi, m_halfedge_graph); e != e_end; e++)
      {
        vertex_descriptor vj = boost::source(*e, m_halfedge_graph);
        std::size_t vj_id = ros_id(vj);

        const CR_vector& pij = sub_to_CR_vector(original[vi_id], original[vj_id]);
        const CR_vector& qij = sub_to_CR_vector(solution[vi_id], solution[vj_id]);
        double wij = edge_weight[id(*e)];

        cr_traits.scalar_vector_vector_transpose_mult(cov, wij, pij, qij); // cov += wij * (pij * qij)
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
      for (boost::tie(e,e_end) = boost::out_edges(vi, m_halfedge_graph); e != e_end; e++)
      {
        if(boost::get(CGAL::edge_is_border, m_halfedge_graph, *e)) { continue; } // no facet 
        // iterate edges around facet
        edge_descriptor edge_around_facet = *e;
        do
        {
          vertex_descriptor v1 = boost::target(edge_around_facet, m_halfedge_graph);
          vertex_descriptor v2 = boost::source(edge_around_facet, m_halfedge_graph);

          std::size_t v1_id = ros_id(v1); std::size_t v2_id = ros_id(v2);
        
          const CR_vector& p12 = sub_to_CR_vector(original[v1_id], original[v2_id]);
          const CR_vector& q12 = sub_to_CR_vector(solution[v1_id], solution[v2_id]);
          double w12 = edge_weight[id(edge_around_facet)];

          cr_traits.scalar_vector_vector_transpose_mult(cov, w12, p12, q12); // cov += w12 * (p12 * q12);

        } while( (edge_around_facet = CGAL::next_edge(edge_around_facet, m_halfedge_graph)) != *e);
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
      for (boost::tie(e,e_end) = boost::in_edges(vi, m_halfedge_graph); e != e_end; e++)
      {
        vertex_descriptor vj = boost::source(*e, m_halfedge_graph);
        std::size_t vj_id = ros_id(vj);

        const CR_vector& pij = sub_to_CR_vector(original[vi_id], original[vj_id]);
        const CR_vector& qij = sub_to_CR_vector(solution[vi_id], solution[vj_id]);
        
        double wij = edge_weight[id(*e)];

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

      if ( is_roi(vi) && !is_control(vi) ) 
      {// free vertices
        // sum of right-hand side of eq:lap_ber in user manual
        CR_vector xyz = cr_traits.vector(0, 0, 0);

        in_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::in_edges(vi, m_halfedge_graph); e != e_end; e++)
        {
          vertex_descriptor vj = boost::source(*e, m_halfedge_graph);
          std::size_t vj_id = ros_id(vj); 
          
          const CR_vector& pij = sub_to_CR_vector(original[vi_id], original[vj_id]);

          double wij = edge_weight[id(*e)];
          double wji = edge_weight[id(CGAL::opposite_edge(*e, m_halfedge_graph))];
#ifndef CGAL_DEFORM_MESH_USE_EXPERIMENTAL_SCALE
          cr_traits.scalar_matrix_scalar_matrix_vector_mult(xyz, wij, rot_mtr[vi_id], wji, rot_mtr[vj_id], pij);
#else
        cr_traits.scalar_matrix_scalar_matrix_vector_mult(xyz, wij * scales[vi_id], rot_mtr[vi_id], 
          wji * scales[vj_id], rot_mtr[vj_id], pij);
#endif
          // corresponds xyz += (wij*rot_mtr[vi_id] + wji*rot_mtr[vj_id]) * pij
        }
        Bx[vi_id] = cr_traits.vector_coeff(xyz, 0); 
        By[vi_id] = cr_traits.vector_coeff(xyz, 1); 
        Bz[vi_id] = cr_traits.vector_coeff(xyz, 2); 
      }
      else 
      {// constrained vertex
        Bx[vi_id] = solution[vi_id].x(); By[vi_id] = solution[vi_id].y(); Bz[vi_id] = solution[vi_id].z();
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

      if ( is_roi(vi) && !is_control(vi) ) 
      {// free vertices
        // sum of right-hand side of eq:lap_ber_rims in user manual
        CR_vector xyz = cr_traits.vector(0, 0, 0);

        out_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::out_edges(vi, m_halfedge_graph); e != e_end; e++)
        {
          vertex_descriptor vj = boost::target(*e, m_halfedge_graph);
          std::size_t vj_id = ros_id(vj); 

          const CR_vector& pij = sub_to_CR_vector(original[vi_id], original[vj_id]);
          
          if(!boost::get(CGAL::edge_is_border, m_halfedge_graph, *e))
          {
            vertex_descriptor vn = boost::target(CGAL::next_edge(*e, m_halfedge_graph), m_halfedge_graph); // opp vertex of e_ij
            double wji = edge_weight[id(*e)] / 3.0;  // edge(pj - pi)           
            cr_traits.scalar_mult_with_matrix_sum(xyz, wji, rot_mtr[vi_id], rot_mtr[vj_id], rot_mtr[ros_id(vn)], pij);
            // corresponds  xyz += wji*(rot_mtr[vi_id] + rot_mtr[vj_id] + rot_mtr[ros_id(vn)])*pij;
          }

          edge_descriptor opp = CGAL::opposite_edge(*e, m_halfedge_graph);
          if(!boost::get(CGAL::edge_is_border, m_halfedge_graph, opp))
          {
            vertex_descriptor vm = boost::target(CGAL::next_edge(opp, m_halfedge_graph), m_halfedge_graph); // other opp vertex of e_ij
            double wij = edge_weight[id(opp)] / 3.0;  // edge(pi - pj)
            cr_traits.scalar_mult_with_matrix_sum(xyz, wij, rot_mtr[vi_id], rot_mtr[vj_id], rot_mtr[ros_id(vm)], pij);
            // corresponds xyz += wij * ( rot_mtr[vi_id] + rot_mtr[vj_id] + rot_mtr[ros_id(vm)] ) * pij
          }
        }
        Bx[vi_id] = cr_traits.vector_coeff(xyz, 0); 
        By[vi_id] = cr_traits.vector_coeff(xyz, 1); 
        Bz[vi_id] = cr_traits.vector_coeff(xyz, 2); 
      }
      else 
      {// constrained vertices
        Bx[vi_id] = solution[vi_id].x(); By[vi_id] = solution[vi_id].y(); Bz[vi_id] = solution[vi_id].z();
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
      solution[v_id] = p;
    }
  }

  /// Assign solution to target mesh
  void assign_solution()
  {
    for(std::size_t i = 0; i < ros.size(); ++i){
      std::size_t v_id = ros_id(ros[i]);
      if(is_roi(ros[i]))
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
      for (boost::tie(e,e_end) = boost::in_edges(vi, m_halfedge_graph); e != e_end; e++)
      {
        vertex_descriptor vj = boost::source(*e, m_halfedge_graph);
        std::size_t vj_id = ros_id(vj);

        const CR_vector& pij = sub_to_CR_vector(original[vi_id], original[vj_id]);
        const CR_vector& qij = sub_to_CR_vector(solution[vi_id], solution[vj_id]);

        double wij = edge_weight[id(*e)];

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
      for (boost::tie(e,e_end) = boost::out_edges(vi, m_halfedge_graph); e != e_end; e++)
      {
        if(boost::get(CGAL::edge_is_border, m_halfedge_graph, *e)) { continue; } // no facet 
        // iterate edges around facet
        edge_descriptor edge_around_facet = *e;
        do
        {
          vertex_descriptor v1 = boost::target(edge_around_facet, m_halfedge_graph);
          vertex_descriptor v2 = boost::source(edge_around_facet, m_halfedge_graph);
          std::size_t v1_id = ros_id(v1); std::size_t v2_id = ros_id(v2);

          const CR_vector& p12 = sub_to_CR_vector(original[v1_id], original[v2_id]);
          const CR_vector& q12 = sub_to_CR_vector(solution[v1_id], solution[v2_id]);
          double w12 = edge_weight[id(edge_around_facet)];
         
          sum_of_energy += w12 * cr_traits.squared_norm_vector_scalar_vector_subs(q12, rot_mtr[vi_id], p12);
          // sum_of_energy += w12 * ( q12 - rot_mtr[vi_id]*p12 )^2

        } while( (edge_around_facet = CGAL::next_edge(edge_around_facet, m_halfedge_graph)) != *e);
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
    return Closest_rotation_traits().vector(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
  }

  template<class Vect>
  Point add_to_point(const Point& p, const Vect& v) {
    return Point(
      p.x() + v[0], p.y() + v[1], p.z() + v[2]);
  }
  template<class Vect>
  Vect sub_to_vector(const Point& p1, const Point& p2) {
    return Vect(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
  }

  /// shorthand of get(vertex_index_map, v)
  std::size_t id(vertex_descriptor vd) const
  { return get(vertex_index_map, vd); }

  std::size_t& ros_id(vertex_descriptor vd)       
  { return ros_id_map[id(vd)]; }
  std::size_t  ros_id(vertex_descriptor vd) const 
  { return ros_id_map[id(vd)]; }

  /// shorthand of get(edge_index_map, e)
  std::size_t id(edge_descriptor e) const
  {
    return get(edge_index_map, e);
  }
};
} //namespace CGAL
#endif  // CGAL_DEFORM_MESH_H
