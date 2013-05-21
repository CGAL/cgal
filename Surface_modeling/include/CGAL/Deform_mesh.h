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

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <CGAL/Default.h>

#include <vector>
#include <list>
#include <utility>
#include <limits>

// for default parameters
#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>  // for sparse linear system solver

  #include <CGAL/Deformation_Eigen_polar_closest_rotation_traits_3.h>  // for 3x3 closest rotation computer

  #if defined(CGAL_SUPERLU_ENABLED)
    #include <Eigen/SuperLUSupport>
  #else
    #include <Eigen/SparseLU>
  #endif

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
template<class Polyhedron, Deformation_algorithm_tag deformation_algorithm_tag>
struct Weight_calculator_selector {
  typedef Uniform_weight<Polyhedron> weight_calculator;
};

template<class Polyhedron>
struct Weight_calculator_selector<Polyhedron, CGAL::SPOKES_AND_RIMS> {
  typedef Single_cotangent_weight<Polyhedron> weight_calculator;
};

template<class Polyhedron>
struct Weight_calculator_selector<Polyhedron, CGAL::ORIGINAL_ARAP> {
  typedef Cotangent_weight<Polyhedron> weight_calculator;
};
}//namespace internal
/// @endcond

 ///
 /// \ingroup PkgSurfaceModeling
 /// @brief Class providing the functionalities for deforming a triangulated surface mesh
 ///
 /// @tparam P a model of HalfedgeGraph 
 /// @tparam VIM a model of `ReadWritePropertyMap`</a>  with Deform_mesh::vertex_descriptor as key and `unsigned int` as value type
 /// @tparam EIM a model of `ReadWritePropertyMap`</a>  with Deform_mesh::edge_descriptor as key and `unsigned int` as value type
 /// @tparam TAG tag for selecting the deformation algorithm
 /// @tparam WC a model of SurfaceModelingWeightCalculator, with `WC::Polyhedron` being `P`
 /// @tparam ST a model of SparseLinearAlgebraTraitsWithPreFactor_d. If \ref thirdpartyEigen "Eigen" 3.1 (or greater) is available 
 /// and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits` is provided as default parameter.\n
 /// If `CGAL_SUPERLU_ENABLED` is defined, the overload is equal to:
 /// \code CGAL::Eigen_solver_traits<Eigen::SuperLU<CGAL::Eigen_sparse_matrix<double>::EigenType> > \endcode
 /// else it is equal to:
 /// \code
 ///     CGAL::Eigen_solver_traits<
 ///         Eigen::SparseLU<
 ///            CGAL::Eigen_sparse_matrix<double, Eigen::ColMajor>::EigenType,
 ///            Eigen::COLAMDOrdering<int> >  >
 /// \endcode
 /// @tparam CR a model of DeformationClosestRotationTraits_3. If \ref thirdpartyEigen "Eigen" 3.1 (or greater) is available and `CGAL_EIGEN3_ENABLED` is defined, 
 /// `Deformation_Eigen_polar_closest_rotation_traits_3` is provided as default parameter.
template <
  class P, 
  class VIM, 
  class EIM,
  Deformation_algorithm_tag TAG = SPOKES_AND_RIMS,
  class WC = Default,
  class ST = Default, 
  class CR = Default  
  >
class Deform_mesh
{
//Typedefs
public:

  /// \name Template parameter types
  /// @{
  // typedefed template parameters, main reason is doxygen creates autolink to typedefs but not template parameters
  typedef P Polyhedron; /**< model of HalfedgeGraph */  
  typedef VIM Vertex_index_map; /**< model of `ReadWritePropertyMap`  with Deform_mesh::vertex_descriptor as key and `unsigned int` as value type */
  typedef EIM Edge_index_map; /**< model of `ReadWritePropertyMap`</a>  with Deform_mesh::edge_descriptor as key and `unsigned int` as value type */

// weight calculator
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    WC,
    typename internal::Weight_calculator_selector<P, TAG>::weight_calculator
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
            CGAL::Eigen_sparse_matrix<double, Eigen::ColMajor>::EigenType,
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
  >::type CR_helper;
#else
  typedef CR CR_helper; /**< model of DeformationClosestRotationTraits_3 */
#endif

  /// @}

  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor	vertex_descriptor; /**< The type for vertex representative objects */
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor		edge_descriptor;   /**< The type for edge representative objects */

  typedef typename Polyhedron::Traits::Vector_3  Vector; /**<The type for Vector_3 from Polyhedron traits */
  typedef typename Polyhedron::Traits::Point_3   Point;  /**<The type for Point_3 from Polyhedron traits */

private:
  typedef Deform_mesh<P, VIM, EIM, TAG, WC, ST, CR> Self;
  // Repeat Polyhedron types
  typedef typename boost::graph_traits<Polyhedron>::vertex_iterator     vertex_iterator;
  typedef typename boost::graph_traits<Polyhedron>::edge_iterator       edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::in_edge_iterator    in_edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::out_edge_iterator   out_edge_iterator;

  // Handle container types
  typedef std::list<vertex_descriptor>  Handle_container;
  typedef std::list<Handle_container>   Handle_group_container;

  typedef typename CR_helper::Matrix CR_matrix;
  typedef typename CR_helper::Vector CR_vector;

public:
  /** The type used as the representative of a group of handles*/
  typedef typename Handle_group_container::iterator                Handle_group;
  /** Const version of Handle_group*/
  typedef typename Handle_group_container::const_iterator          Const_handle_group;
  /** Iterator over the groups of handles. The type can be implicitly 
      converted to Deform_mesh::Handle_group or Deform_mesh::Const_handle_group */
  typedef typename Handle_group_container::iterator                Handle_group_iterator;
   /** Const version of Handle_group_iterator */
  typedef typename Handle_group_container::const_iterator          Handle_group_const_iterator;

  /** Iterator over vertex descriptors in a group of handles. Its value type is `vertex_descriptor` */
  typedef typename Handle_container::iterator                      Handle_iterator;
   /** Const version of Handle_iterator*/
  typedef typename Handle_container::const_iterator                Handle_const_iterator;

  /** Iterator over vertex descriptors in the region-of-interest. Its value type is `vertex_descriptor` */
  typedef typename std::vector<vertex_descriptor>::iterator        Roi_iterator;
   /** Const version of Roi_iterator*/
  typedef typename std::vector<vertex_descriptor>::const_iterator  Roi_const_iterator;

// Data members.
private:
  Polyhedron& polyhedron;															/**< Source triangulated surface mesh for modeling */

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
  Handle_group_container handle_group_list;           ///< user specified handles

  Weight_calculator weight_calculator;
private:
  Deform_mesh(const Self& s) { } 

// Public methods
public:
/// \name Preprocess Section
/// @{
  /**
   * The constructor of a deformation object
   *
   * @pre the polyhedron consists of only triangular facets
   * @param polyhedron triangulated surface mesh used to deform
   * @param vertex_index_map property map for associating an id to each vertex
   * @param edge_index_map property map for associating an id to each edge
   * @param iterations see `set_iterations()` for more details
   * @param tolerance  see `set_tolerance()` for more details
   * @param weight_calculator function object or pointer for weight calculation
   */
  Deform_mesh(Polyhedron& polyhedron, 
              Vertex_index_map vertex_index_map, 
              Edge_index_map edge_index_map,
              unsigned int iterations = 5,
              double tolerance = 1e-4,
              Weight_calculator weight_calculator = Weight_calculator())
    : polyhedron(polyhedron), vertex_index_map(vertex_index_map), edge_index_map(edge_index_map),
      ros_id_map(std::vector<std::size_t>(boost::num_vertices(polyhedron), (std::numeric_limits<std::size_t>::max)() )),
      is_roi_map(std::vector<bool>(boost::num_vertices(polyhedron), false)),
      is_hdl_map(std::vector<bool>(boost::num_vertices(polyhedron), false)),
      iterations(iterations), tolerance(tolerance),
      need_preprocess_factorization(true), 
      need_preprocess_region_of_solution(true),
      last_preprocess_successful(false),
      weight_calculator(weight_calculator)
  {
    // assign id to each vertex and edge
    vertex_iterator vb, ve;
    std::size_t id = 0;
    for(boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb, ++id)
    {
      put(vertex_index_map, *vb, id);
    }

    edge_iterator eb, ee;
    id = 0;
    for(boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb, ++id)
    {
      put(edge_index_map, *eb, id);
    }

    // compute edge weights
    edge_weight.reserve(boost::num_edges(polyhedron));
    for(boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      edge_weight.push_back(this->weight_calculator(*eb, polyhedron));
    }
  }

  /**
   * Puts the object in the same state as after the creation (except iterations and tolerance).
   */
  void reset()
  {
    need_preprocess_both();
    // clear vertices
    roi.clear();
    handle_group_list.clear();
    is_roi_map.assign(boost::num_vertices(polyhedron), false);
    is_hdl_map.assign(boost::num_vertices(polyhedron), false);
  }

  /**
   * Creates a new empty group of handles.
   * `insert_handle(Handle_group handle_group, vertex_descriptor vd)` or `insert_handle(Handle_group handle_group, InputIterator begin, InputIterator end)`
   * must be used to insert handles in a group.
   * After inserting vertices, use `translate()` or `rotate()` for applying a transformation on all the vertices inside a group.
   * @return a representative of the group of handles created (it is valid until `erase_handle(Handle_group handle_group)` is called)
   */
  Handle_group create_handle_group()
  {
    // no need to need_preprocess = true;
    handle_group_list.push_back(Handle_container());
    return --handle_group_list.end();
  }
  
  /**
   * Inserts a vertex into a group of handles. The vertex is also inserted in the region-of-interest if it is not already in it.
   * @param handle_group the group where the vertex is inserted in
   * @param vd the vertex to be inserted
   * @return `true` if the insertion is successful
   */
  bool insert_handle(Handle_group handle_group, vertex_descriptor vd)
  {
    if(is_handle(vd)) { return false; }
    need_preprocess_both();

    insert_roi(vd); // also insert it as roi

    is_handle(vd) = true;    
    handle_group->push_back(vd);
    return true;
  }

  /**
   * Inserts a range of vertices in a group of handles. The vertices are also inserted in the region-of-interest if they are not already in it.
   * @tparam InputIterator input iterator type with `vertex_descriptor` as value type
   * @param handle_group the group where the vertex is inserted in
   * @param begin first iterator of the range of vertices
   * @param end past-the-end iterator of the range of vertices
   */
  template<class InputIterator>
  void insert_handle(Handle_group handle_group, InputIterator begin, InputIterator end)
  {
    for( ;begin != end; ++begin)
    {
      insert_handle(handle_group, *begin);
    }
  }

  /**
   * Erases a group of handles. Its representative becomes invalid.
   * @param handle_group the group to be erased
   */
  void erase_handle(Handle_group handle_group)
  {
    need_preprocess_both();
    for(typename Handle_container::iterator it = handle_group->begin(); it != handle_group->end(); ++it) 
    {
      is_handle(*it) = false;
    }
    handle_group_list.erase(handle_group);
  }

  /**
   * Erases a vertex from a group of handles.
   * \note The group of handles is not erased even if it becomes empty.
   * @param handle_group the group where the vertex is erased from
   * @param vd the vertex to be erased
   * @return `true` if the removal is successful
   */
  bool erase_handle(Handle_group handle_group, vertex_descriptor vd)
  {
    if(!is_handle(vd)) { return false; }
    
    typename Handle_container::iterator it = std::find(handle_group->begin(), handle_group->end(), vd);
    if(it != handle_group->end())
    {
      is_handle(*it) = false;
      handle_group->erase(it);   
      need_preprocess_both();
      return true;
      // Although the handle group might get empty, we do not delete it from handle_group
    }
    return false; // OK (vd is handle but placed in other handle group than handle_group argument)
  }

  /**
   * Erases a vertex by searching it through all groups of handles.
   * \note The group of handles is not erased even if it becomes empty.
   * @param vd the vertex to be erased
   * @return `true` if the removal is successful
   */
  bool erase_handle(vertex_descriptor vd)
  {
    if(!is_handle(vd)) { return false; }

    for(Handle_group it = handle_group_list.begin(); it != handle_group_list.end(); ++it)
    {
      if(erase_handle(it, vd)) { return true; }
    }

    CGAL_assertion(false);// inconsistency between is_handle_map, and handle_group_list
    return false;
  }

  /** 
   * Provides access to all the groups of handles.
   * @return the range of iterators over all groups of handles
   */
  std::pair<Handle_group_iterator, Handle_group_iterator> handle_groups()
  {
    return std::make_pair(handle_group_list.begin(), handle_group_list.end());
  }

  /** 
   * Const version
   */
  std::pair<Handle_group_const_iterator, Handle_group_const_iterator> handle_groups() const
  {
    return std::make_pair(handle_group_list.begin(), handle_group_list.end());
  }

  /** 
   * Provides access to all the handles inside a group.
   * @param handle_group the group containing handles
   * @return the range of iterators over all handles inside a group
   * 
   */
  std::pair<Handle_iterator, Handle_iterator> handles(Handle_group handle_group)
  {
    return std::make_pair(handle_group->begin(), handle_group->end());
  }

  /** 
   * Const version
   */
  std::pair<Handle_const_iterator, Handle_const_iterator> handles(Const_handle_group handle_group) const
  {
    return std::make_pair(handle_group->begin(), handle_group->end());
  }

  /**
   * Inserts a range of vertices in the region-of-interest
   * @tparam InputIterator input iterator with `vertex_descriptor` as value type
   * @param begin first iterator of the range of vertices
   * @param end past-the-end iterator of the range of vertices
   */
  template<class InputIterator>
  void insert_roi(InputIterator begin, InputIterator end)
  {
    for( ;begin != end; ++begin)
    {
      insert_roi(*begin);
    }
  }

  /**
   * Inserts a vertex in the region-of-interest
   * @param vd the vertex to be inserted
   * @return `true` if the insertion is successful
   */
  bool insert_roi(vertex_descriptor vd)   
  {
    if(is_roi(vd)) { return false; }
    need_preprocess_both();

    is_roi(vd) = true;
    roi.push_back(vd);
    return true;
  }

  /**
   * Erases a vertex from the region-of-interest. The vertex is also removed from any group of handles.
   * \note The next call to `preprocess()`, any vertex which is no longer in the region-of-interest will be assigned to its original position 
   * (that is position of the vertex at the time of construction or after the last call to `overwrite_original_positions()`).
   * @param vd the vertex to be erased
   * @return `true` if the removal is successful
   */
  bool erase_roi(vertex_descriptor vd)   
  {
    if(!is_roi(vd)) { return false; }  
    
    erase_handle(vd); // also erase from being handle

    typename std::vector<vertex_descriptor>::iterator it = std::find(roi.begin(), roi.end(), vd);
    if(it != roi.end())
    {
      is_roi(vd) = false;
      roi.erase(it);

      need_preprocess_both();
      return true;
    }
    
    CGAL_assertion(false); // inconsistency between is_roi_map, and roi vector!
    return false;
  }

  /** 
   * Provides access to the vertices in the region-of-interest. 
   * \note Deleting a vertex from the region-of-interest invalidates iterators. 
   * @return an iterator range
   */
  std::pair<Roi_iterator, Roi_iterator> roi_vertices()
  {
    return std::make_pair(roi.begin(), roi.end());
  }

  /** 
   * Const version
   */
  std::pair<Roi_const_iterator, Roi_const_iterator> roi_vertices() const
  {
    return std::make_pair(roi.begin(), roi.end());
  }

  /**
   * Triggers the necessary precomputation work before beginning deformation.
   * \note Calling this function is optional.
   * \note The insertion / removal of a vertex in a group of handles or in the region-of-interest invalidates the
   * preprocessing data.
   * @return `true` if Laplacian matrix factorization is successful.
   * A common reason for failure is that the system is rank deficient, 
   * which happens if there is no path between a free vertex and a handle vertex (i.e. both fixed and user-inserted).
   */
  bool preprocess()
  {
    region_of_solution();
    assemble_laplacian_and_factorize();
    return last_preprocess_successful; // which is set by assemble_laplacian_and_factorize()
  }
/// @} Preprocess Section

/// \name Deform Section
/// @{  
  /**
   * Sets the transformation to apply to all the vertices in a group of handles to be a translation by vector `t`.
   * \note This transformation is applied on the original positions of the vertices 
   * (that is positions of vertices at the time of construction or after the last call to `overwrite_original_positions()`). 
   * \note A call to this function cancels the last call to `rotate()`, `translate()`, or `assign()`.
   * @param handle_group the representative of the group of handles to be translated
   * @param t translation vector 
   */
  void translate(Const_handle_group handle_group, const Vector& t)
  {
    region_of_solution(); // we require ros ids, so if there is any need to preprocess of region of solution -do it.

    for(typename Handle_container::const_iterator it = handle_group->begin();
      it != handle_group->end(); ++it)
    {
      std::size_t v_id = ros_id(*it);
      solution[v_id] = original[v_id] + t;
    }
  }

  /**
   * Sets the transformation to apply to all the vertices in a group of handles to be a rotation around `rotation_center`
   * defined by the quaternion `quat`, followed by a translation by vector `t`.
   * \note This transformation is applied on the original positions of the vertices 
   * (that is positions of vertices at the time of construction or after the last call to `overwrite_original_positions()`).  
   * \note A call to this function cancels the last call to `rotate()`, `translate()`, or `assign()`.
   * @tparam Quaternion is a quaternion class with `Vect operator*(Quaternion, Vect)` being defined and returns the product of a quaternion with a vector
   * @tparam Vect is a 3D vector class, `Vect(double x,double y, double z)` being a constructor available and `Vect::operator[](int i)` with i=0,1 or 2 returns its coordinates
   * @param handle_group the representative of the group of handles to be rotated and translated
   * @param rotation_center center of rotation
   * @param quat rotation holder quaternion
   * @param t post translation vector
   */
  template <typename Quaternion, typename Vect>
  void rotate(Const_handle_group handle_group, const Point& rotation_center, const Quaternion& quat, const Vect& t)
  {
    region_of_solution(); // we require ros ids, so if there is any need to preprocess of region of solution -do it.

    for(typename Handle_container::const_iterator it = handle_group->begin();
      it != handle_group->end(); ++it)
    {
      std::size_t v_id = ros_id(*it);

      Point p = CGAL::ORIGIN + ( original[v_id] - rotation_center);
      Vect v = quat * Vect(p.x(),p.y(),p.z());
      p = Point(v[0], v[1], v[2]) + (rotation_center - CGAL::ORIGIN);
      p = p + Vector(t[0],t[1],t[2]);

      solution[v_id] = p;
    }
  }

  /**
   * Assigns the target position of a handle vertex 
   * @param vd the handle vertex to be assigned target position
   * @param target_position the new target position
   */
  void assign(vertex_descriptor vd, const Point& target_position)
  {
    region_of_solution(); // we require ros ids, so if there is any need to preprocess of region of solution -do it.

    if(!is_handle(vd)) { return; }
    solution[ros_id(vd)] = target_position;
  }

  /**
   * Deforms the region-of-interest according to the deformation algorithm, applying for each group of handles the transformation provided by `rotate()` or `translate()`
   * to their original positions, or using target positions provided by `assign()`. 
   * The coordinates of the vertices of the input graph that are inside the region-of-interest are updated. The initial guess for solving the
   * deformation problem is using the coordinates of the input graph before calling the function.
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
      optimal_rotations(); 

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
/// @} Deform Section

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
   * @return `true` if the vertex is inside the region-of-interest
   */
  bool is_roi(vertex_descriptor vd) const
  { return is_roi_map[id(vd)]; }

  /**
   * Queries whether a vertex is a handle.
   * @param vd the query vertex
   * @return `true` if the vertex is inside any group of handles
   */
  bool is_handle(vertex_descriptor vd) const
  { return is_hdl_map[id(vd)]; }

  /**
   * Provides access to halfedge graph being deformed
   * @return the halfedge graph
   */
  const Polyhedron& halfedge_graph() const
  { return polyhedron; }
  
  /**
   * Sets the original positions to be the current positions for vertices inside region-of-interest. Calling this function has the same effect as creating
   * a new deformation object with the current deformed polyhedron, keeping the region-of-interest and the groups of handles.
   * \note if the region-of-interest or any group of handles have been modified since the last call to `preprocess()`,
   * it will be called prior to the overwrite.
   */
  void overwrite_original_positions()
  {
    if(roi.empty()) { return; } // no ROI to overwrite

    region_of_solution(); // the roi should be preprocessed since we are using original_position vec

    Roi_iterator rb, re;
    for(boost::tie(rb, re) = roi_vertices(); rb != re; ++rb)
    {
      original[ros_id(*rb)] = (*rb)->point();
    }

    // now I need to compute weights for edges incident to roi vertices
    std::vector<bool> is_weight_computed(boost::num_edges(polyhedron), false);
    for(boost::tie(rb, re) = roi_vertices(); rb != re; ++rb)
    {
      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(*rb, polyhedron); e != e_end; e++)
      {
        std::size_t id_e = id(*e);
        if(is_weight_computed[id_e]) { continue; }

        edge_weight[id_e] = weight_calculator(*e, polyhedron);
        is_weight_computed[id_e] = true;

        edge_descriptor e_opp = CGAL::opposite_edge(*e, polyhedron);
        std::size_t id_e_opp = id(e_opp);

        edge_weight[id_e_opp] = weight_calculator(e_opp, polyhedron);
        is_weight_computed[id_e_opp] = true;
      }
    }    

    // also set rotation matrix to identity
    std::fill(rot_mtr.begin(), rot_mtr.end(), CR_helper().identity_matrix());

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
    for (boost::tie(e,e_end) = boost::in_edges(vd, polyhedron); e != e_end; e++)
    {
      vertex_descriptor vt = boost::source(*e, polyhedron);    
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
        (*it)->point() = old_original[ old_ros_id_map[id(*it)] ];
      }
    }

    ////////////////////////////////////////////////////////////////
    // assign id to vertices inside: roi, boundary of roi (roi + boundary of roi = ros),
    //                               and boundary of ros
    // keep in mind that id order does not have to be compatible with ros order
    ros.clear(); // clear ros    
    ros.insert(ros.end(), roi.begin(), roi.end()); 

    ros_id_map.assign(boost::num_vertices(polyhedron), (std::numeric_limits<std::size_t>::max)()); // use max as not assigned mark

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
        rot_mtr[v_ros_id] = CR_helper().identity_matrix();
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
        solution[v_ros_id] = ros[i]->point();
        original[v_ros_id] = ros[i]->point();
      }         
    }

    for(std::size_t i = 0; i < outside_ros.size(); ++i)
    {
      std::size_t v_ros_id = ros_id(outside_ros[i]);
      original[v_ros_id] = outside_ros[i]->point();
      solution[v_ros_id] = outside_ros[i]->point();
    }
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
      if ( is_roi(vi) && !is_handle(vi) )          // vertices of ( roi - hdl )
      {
        double diagonal = 0;
        in_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
        {
          vertex_descriptor vj = boost::source(*e, polyhedron);
          double wij = edge_weight[id(*e)];  // edge(pi - pj)
          double wji = edge_weight[id(CGAL::opposite_edge(*e, polyhedron))]; // edge(pi - pj)
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
      if ( is_roi(vi) && !is_handle(vi) ) // vertices of ( roi - hdl ): free vertices
      {
        double diagonal = 0;
        out_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::out_edges(vi, polyhedron); e != e_end; e++)
        {
          double total_weight = 0;
          // an edge contribute to energy only if it is part of an incident triangle 
          // (i.e it should not be a border edge)
          if(!boost::get(CGAL::edge_is_border, polyhedron, *e)) 
          {
            double wji = edge_weight[id(*e)]; // edge(pj - pi)
            total_weight += wji; 
          }

          edge_descriptor opp = CGAL::opposite_edge(*e, polyhedron);
          if(!boost::get(CGAL::edge_is_border, polyhedron, opp))
          {
            double wij = edge_weight[id(opp)]; // edge(pi - pj)
            total_weight += wij;
          }

          // place coefficient to matrix
          vertex_descriptor vj = boost::target(*e, polyhedron);
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
    CR_helper cr_helper;
    CR_matrix cov = cr_helper.zero_matrix();

    // only accumulate ros vertices
    for ( std::size_t k = 0; k < ros.size(); k++ )
    {
      vertex_descriptor vi = ros[k];
      std::size_t vi_id = ros_id(vi);
      // compute covariance matrix (user manual eq:cov_matrix)
      cov = cr_helper.zero_matrix();

      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
      {
        vertex_descriptor vj = boost::source(*e, polyhedron);
        std::size_t vj_id = ros_id(vj);

        const CR_vector& pij = sub_to_CR_vector(original[vi_id], original[vj_id]);
        const CR_vector& qij = sub_to_CR_vector(solution[vi_id], solution[vj_id]);
        double wij = edge_weight[id(*e)];

        cr_helper.scalar_vector_vector_transpose_mult(cov, wij, pij, qij); // cov += wij * (pij * qij)
      }

      cr_helper.compute_close_rotation(cov, rot_mtr[vi_id]);
    }
  }

  void optimal_rotations_spokes_and_rims()
  {    
    CR_helper cr_helper;
    CR_matrix cov =cr_helper.zero_matrix();   

    // only accumulate ros vertices
    for ( std::size_t k = 0; k < ros.size(); k++ )
    {
      vertex_descriptor vi = ros[k];
      std::size_t vi_id = ros_id(vi);
      // compute covariance matrix
      cov = cr_helper.zero_matrix();

      //iterate through all triangles 
      out_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::out_edges(vi, polyhedron); e != e_end; e++)
      {
        if(boost::get(CGAL::edge_is_border, polyhedron, *e)) { continue; } // no facet 
        // iterate edges around facet
        edge_descriptor edge_around_facet = *e;
        do
        {
          vertex_descriptor v1 = boost::target(edge_around_facet, polyhedron);
          vertex_descriptor v2 = boost::source(edge_around_facet, polyhedron);

          std::size_t v1_id = ros_id(v1); std::size_t v2_id = ros_id(v2);
        
          const CR_vector& p12 = sub_to_CR_vector(original[v1_id], original[v2_id]);
          const CR_vector& q12 = sub_to_CR_vector(solution[v1_id], solution[v2_id]);
          double w12 = edge_weight[id(edge_around_facet)];

          cr_helper.scalar_vector_vector_transpose_mult(cov, w12, p12, q12); // cov += w12 * (p12 * q12);

        } while( (edge_around_facet = CGAL::next_edge(edge_around_facet, polyhedron)) != *e);
      }

      cr_helper.compute_close_rotation(cov, rot_mtr[vi_id]);
    }
  }

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

    CR_helper cr_helper;

    // assemble right columns of linear system
    for ( std::size_t k = 0; k < ros.size(); k++ )
    {
      vertex_descriptor vi = ros[k];
      std::size_t vi_id = ros_id(vi);

      if ( is_roi(vi) && !is_handle(vi) ) 
      {// free vertices
        // sum of right-hand side of eq:lap_ber in user manual
        CR_vector xyz = cr_helper.vector(0, 0, 0);

        in_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
        {
          vertex_descriptor vj = boost::source(*e, polyhedron);
          std::size_t vj_id = ros_id(vj); 
          
          const CR_vector& pij = sub_to_CR_vector(original[vi_id], original[vj_id]);

          double wij = edge_weight[id(*e)];
          double wji = edge_weight[id(CGAL::opposite_edge(*e, polyhedron))];

          cr_helper.scalar_matrix_scalar_matrix_vector_mult(xyz, wij, rot_mtr[vi_id], wji, rot_mtr[vj_id], pij);
          // corresponds xyz += (wij*rot_mtr[vi_id] + wji*rot_mtr[vj_id]) * pij
        }
        Bx[vi_id] = cr_helper.vector_coeff(xyz, 0); 
        By[vi_id] = cr_helper.vector_coeff(xyz, 1); 
        Bz[vi_id] = cr_helper.vector_coeff(xyz, 2); 
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

    CR_helper cr_helper;

    // assemble right columns of linear system
    for ( std::size_t k = 0; k < ros.size(); k++ )
    {
      vertex_descriptor vi = ros[k];
      std::size_t vi_id = ros_id(vi);

      if ( is_roi(vi) && !is_handle(vi) ) 
      {// free vertices
        // sum of right-hand side of eq:lap_ber_rims in user manual
        CR_vector xyz = cr_helper.vector(0, 0, 0);

        out_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::out_edges(vi, polyhedron); e != e_end; e++)
        {
          vertex_descriptor vj = boost::target(*e, polyhedron);
          std::size_t vj_id = ros_id(vj); 

          const CR_vector& pij = sub_to_CR_vector(original[vi_id], original[vj_id]);
          
          if(!boost::get(CGAL::edge_is_border, polyhedron, *e))
          {
            vertex_descriptor vn = boost::target(CGAL::next_edge(*e, polyhedron), polyhedron); // opp vertex of e_ij
            double wji = edge_weight[id(*e)] / 3.0;  // edge(pj - pi)           
            cr_helper.scalar_mult_with_matrix_sum(xyz, wji, rot_mtr[vi_id], rot_mtr[vj_id], rot_mtr[ros_id(vn)], pij);
            // corresponds  xyz += wji*(rot_mtr[vi_id] + rot_mtr[vj_id] + rot_mtr[ros_id(vn)])*pij;
          }

          edge_descriptor opp = CGAL::opposite_edge(*e, polyhedron);
          if(!boost::get(CGAL::edge_is_border, polyhedron, opp))
          {
            vertex_descriptor vm = boost::target(CGAL::next_edge(opp, polyhedron), polyhedron); // other opp vertex of e_ij
            double wij = edge_weight[id(opp)] / 3.0;  // edge(pi - pj)
            cr_helper.scalar_mult_with_matrix_sum(xyz, wij, rot_mtr[vi_id], rot_mtr[vj_id], rot_mtr[ros_id(vm)], pij);
            // corresponds xyz += wij * ( rot_mtr[vi_id] + rot_mtr[vj_id] + rot_mtr[ros_id(vm)] ) * pij
          }
        }
        Bx[vi_id] = cr_helper.vector_coeff(xyz, 0); 
        By[vi_id] = cr_helper.vector_coeff(xyz, 1); 
        Bz[vi_id] = cr_helper.vector_coeff(xyz, 2); 
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
        ros[i]->point() = solution[v_id];
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
    CR_helper cr_helper;

    double sum_of_energy = 0;    
    // only accumulate ros vertices
    for( std::size_t k = 0; k < ros.size(); k++ )
    {
      vertex_descriptor vi = ros[k];
      std::size_t vi_id = ros_id(vi);

      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
      {
        vertex_descriptor vj = boost::source(*e, polyhedron);
        std::size_t vj_id = ros_id(vj);

        const CR_vector& pij = sub_to_CR_vector(original[vi_id], original[vj_id]);
        const CR_vector& qij = sub_to_CR_vector(solution[vi_id], solution[vj_id]);

        double wij = edge_weight[id(*e)];

        sum_of_energy += wij * cr_helper.squared_norm_vector_scalar_vector_subs(qij, rot_mtr[vi_id], pij);
        // sum_of_energy += wij * ( qij - rot_mtr[vi_id]*pij )^2
      }
    }
    return sum_of_energy;
  }
  double energy_spokes_and_rims() const
  {
    CR_helper cr_helper;

    double sum_of_energy = 0;
    // only accumulate ros vertices
    for( std::size_t k = 0; k < ros.size(); k++ )
    {
      vertex_descriptor vi = ros[k];
      std::size_t vi_id = ros_id(vi);
      //iterate through all triangles 
      out_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::out_edges(vi, polyhedron); e != e_end; e++)
      {
        if(boost::get(CGAL::edge_is_border, polyhedron, *e)) { continue; } // no facet 
        // iterate edges around facet
        edge_descriptor edge_around_facet = *e;
        do
        {
          vertex_descriptor v1 = boost::target(edge_around_facet, polyhedron);
          vertex_descriptor v2 = boost::source(edge_around_facet, polyhedron);
          std::size_t v1_id = ros_id(v1); std::size_t v2_id = ros_id(v2);

          const CR_vector& p12 = sub_to_CR_vector(original[v1_id], original[v2_id]);
          const CR_vector& q12 = sub_to_CR_vector(solution[v1_id], solution[v2_id]);
          double w12 = edge_weight[id(edge_around_facet)];
         
          sum_of_energy += w12 * cr_helper.squared_norm_vector_scalar_vector_subs(q12, rot_mtr[vi_id], p12);
          // sum_of_energy += w12 * ( q12 - rot_mtr[vi_id]*p12 )^2

        } while( (edge_around_facet = CGAL::next_edge(edge_around_facet, polyhedron)) != *e);
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
    return CR_helper().vector(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
  }

  /// shorthand of get(vertex_index_map, v)
  std::size_t id(vertex_descriptor vd) const
  { return get(vertex_index_map, vd); }

  std::size_t& ros_id(vertex_descriptor vd)       
  { return ros_id_map[id(vd)]; }
  std::size_t  ros_id(vertex_descriptor vd) const 
  { return ros_id_map[id(vd)]; }

  std::vector<bool>::reference is_handle(vertex_descriptor vd)
  { return is_hdl_map[id(vd)]; }

  std::vector<bool>::reference is_roi(vertex_descriptor vd)
  { return is_roi_map[id(vd)]; }

  /// shorthand of get(edge_index_map, e)
  std::size_t id(edge_descriptor e) const
  {
    return get(edge_index_map, e);
  }

  #ifdef CGAL_DEFORM_EXPERIMENTAL      // Experimental stuff, needs further testing

  double norm_1(const Eigen::Matrix3d& X)
  {
    double sum = 0;
    for ( int i = 0; i < 3; i++ )
    {
      for ( int j = 0; j < 3; j++ )
      {
        sum += abs(X(i,j));
      }
    }
    return sum;
  }

  double norm_inf(const Eigen::Matrix3d& X)
  {
    double max_abs = abs(X(0,0));
    for ( int i = 0; i < 3; i++ )
    {
      for ( int j = 0; j < 3; j++ )
      {
        double new_abs = abs(X(i,j));
        if ( new_abs > max_abs )
        {
          max_abs = new_abs;
        }
      }
    }
    return max_abs;
  }

  // polar decomposition using Newton's method, with warm start, stable but slow
  // not used, need to be investigated later
  void polar_newton(const Eigen::Matrix3d& A, Eigen::Matrix3d &U, double tole)
  {
    Eigen::Matrix3d X = A;
    Eigen::Matrix3d Y;
    double alpha, beta, gamma;
    do 
    {
      Y = X.inverse();
      alpha = sqrt( norm_1(X) * norm_inf(X) );
      beta = sqrt( norm_1(Y) * norm_inf(Y) );
      gamma = sqrt(beta/alpha);
      X = 0.5*( gamma*X + Y.transpose()/gamma );

    } while ( abs(gamma-1) > tole );

    U = X;
  }
  
  // polar decomposition using Eigen, 5 times faster than SVD
  template<typename Mat>
  void polar_eigen(const Mat& A, Mat& R, bool& SVD)
  {
    typedef typename Mat::Scalar Scalar;
    typedef Eigen::Matrix<typename Mat::Scalar,3,1> Vec;

    const Scalar th = std::sqrt(Eigen::NumTraits<Scalar>::dummy_precision());

    Eigen::SelfAdjointEigenSolver<Mat> eig;
    feclearexcept(FE_UNDERFLOW);
    eig.computeDirect(A.transpose()*A);
    if(fetestexcept(FE_UNDERFLOW) || eig.eigenvalues()(0)/eig.eigenvalues()(2)<th)
    {
      // The computation of the eigenvalues might have diverged.
      // Fallback to an accurate SVD based decompositon method.
      Eigen::JacobiSVD<Mat> svd;
      svd.compute(A, Eigen::ComputeFullU | Eigen::ComputeFullV );
      const Mat& u = svd.matrixU(); const Mat& v = svd.matrixV();
      R = u*v.transpose();
      SVD = true;
      return;
    }

    Vec S = eig.eigenvalues().cwiseSqrt();
    R = A  * eig.eigenvectors() * S.asDiagonal().inverse()
      * eig.eigenvectors().transpose();
    SVD = false;

    if(std::abs(R.squaredNorm()-3.) > th)
    {
      // The computation of the eigenvalues might have diverged.
      // Fallback to an accurate SVD based decomposition method.
      Eigen::JacobiSVD<Mat> svd;
      svd.compute(A, Eigen::ComputeFullU | Eigen::ComputeFullV );
      const Mat& u = svd.matrixU(); const Mat& v = svd.matrixV();
      R = u*v.transpose();
      SVD = true;
      return;
    }
  }

  // Local step of iterations, computing optimal rotation matrices using Polar decomposition
  void optimal_rotations_polar()
  {
    Eigen::Matrix3d u, v;           // orthogonal matrices 
    Eigen::Vector3d w;              // singular values
    Eigen::Matrix3d cov;            // covariance matrix
    Eigen::Matrix3d r;
    Eigen::JacobiSVD<Eigen::Matrix3d> svd;      // SVD solver, for non-positive covariance matrices
    int num_svd = 0;
    bool SVD = false;

    // only accumulate ros vertices
    for ( std::size_t i = 0; i < ros.size(); i++ )
    {
      vertex_descriptor vi = ros[i];
      // compute covariance matrix
      cov.setZero();

      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
      {
        vertex_descriptor vj = boost::source(*e, polyhedron);
        Vector pij = original[get(vertex_index_map, vi)] - original[get(vertex_index_map, vj)];
        Vector qij = solution[get(vertex_index_map, vi)] - solution[get(vertex_index_map, vj)];
        double wij = edge_weight[get(edge_index_map, *e)];
        for (int j = 0; j < 3; j++)
        {
          for (int k = 0; k < 3; k++)
          {
            cov(j, k) += wij*pij[j]*qij[k]; 
          }
        }
      }

      // svd decomposition
      if (cov.determinant() > 0)
      {
        polar_eigen<Eigen::Matrix3d> (cov, r, SVD);
        //polar_newton(cov, r, 1e-4);   
        if(SVD)
          num_svd++;
        r.transposeInPlace();     // the optimal rotation matrix should be transpose of decomposition result
      }
      else
      {
        svd.compute( cov, Eigen::ComputeFullU | Eigen::ComputeFullV );
        u = svd.matrixU(); v = svd.matrixV(); w = svd.singularValues();
        r = v*u.transpose();
        num_svd++;
      }
      
      // checking negative determinant of covariance matrix
      if ( r.determinant() < 0 )    // back to SVD method
      {
        if (cov.determinant() > 0)
        {
          svd.compute( cov, Eigen::ComputeFullU | Eigen::ComputeFullV );
          u = svd.matrixU(); v = svd.matrixV(); w = svd.singularValues();
          num_svd++;
        }
        for (int j = 0; j < 3; j++)
        {
          int j0 = j;
          int j1 = (j+1)%3;
          int j2 = (j1+1)%3;
          if ( w[j0] <= w[j1] && w[j0] <= w[j2] )    // smallest singular value as j0
          {
            u(0, j0) = - u(0, j0);
            u(1, j0) = - u(1, j0);
            u(2, j0) = - u(2, j0);
            break;
          }
        }

        // re-extract rotation matrix
        r = v*u.transpose();
      }

      rot_mtr[i] = r;
    }

    double svd_percent = (double)(num_svd)/ros.size();
    CGAL_TRACE_STREAM << svd_percent*100 << "% percentage SVD decompositions;";
    CGAL_TRACE_STREAM << num_svd << " SVD decompositions\n";

  }

#endif
};
} //namespace CGAL
#endif  // CGAL_DEFORM_MESH_H
