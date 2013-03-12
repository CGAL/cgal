// Copyright (c) 2011 GeometryFactory
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
// Author(s)     : Yin Xu, Andreas Fabri

#ifndef CGAL_DEFORM_MESH_H
#define CGAL_DEFORM_MESH_H

#include <CGAL/internal/Surface_modeling/Weights.h>
#include <CGAL/internal/Surface_modeling/Spokes_and_rims_iterator.h>

#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <CGAL/FPU_extension.h>

#include <Eigen/Eigen>
#include <Eigen/SVD>

#include <set>
#include <vector>
#include <list>

//#define CGAL_DEFORM_SPOKES_AND_RIMS // this one is my understanding from spokes and rims
                                      // it uses spoke and rim edges in energy eq

#define CGAL_DEFORM_SPOKES_AND_RIMS_2 // this one is from implementation of Olga (again my understanding)
                                      // it uses all edges in facets around a vertex
namespace CGAL {

/// \ingroup PkgSurfaceModeling
/**
 * @brief Class providing the functionalities for deforming a triangulated surface mesh
 *
 * @pre @a polyhedron.is_pure_triangle()
 * @tparam Polyhedron a model of HalfedgeGraph
 * @tparam SparseLinearAlgebraTraitsWithPreFactor_d sparse linear solver for square sparse linear systems
 * @tparam VertexIndexMap a <a href="http://www.boost.org/doc/libs/release/libs/property_map/doc/ReadWritePropertyMap.html">`ReadWritePropertyMap`</a>  with vertex_descriptor as key and `unsigned int` as value type
 * @tparam EdgeIndexMap a <a href="http://www.boost.org/doc/libs/release/libs/property_map/doc/ReadWritePropertyMap.html">`ReadWritePropertyMap`</a>  with edge_descriptor as key and `unsigned int` as value type
 * @tparam WeightCalculator how to document this (should I provide a concept, like in SegmentationGeomTraits ?) */
 /// \code
 /// // a simple model to WeightCalculator concept, which provides uniform weights
 /// template<class Polyhedron>
 /// class Uniform_weight
 /// {
 /// public:
 ///   typedef typename boost::graph_traits<Polyhedron>::edge_descriptor edge_descriptor;
 ///
 ///   Uniform_weight(Polyhedron& /*polyhedron*/) { } 
 ///
 ///   double operator()(edge_descriptor e)
 ///   { return 1.0; }
 /// };
 /// \endcode
 
template <
  class Polyhedron, 
  class SparseLinearAlgebraTraits_d, 
  class VertexIndexMap, 
  class EdgeIndexMap,
  class WeightCalculator = internal::Single_cotangent_weight<Polyhedron >
  >
class Deform_mesh
{
//Typedefs
public:
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor	vertex_descriptor; /**< The type for vertex representative objects */
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor		edge_descriptor;   /**< The type for edge representative objects */

  typedef typename Polyhedron::Traits::Vector_3           Vector; /**<The type for Vector_3 from Polyhedron traits */
  typedef typename Polyhedron::Traits::Point_3            Point;  /**<The type for Point_3 from Polyhedron traits */

private:
  // Geometric types              
  typedef typename Polyhedron::Traits         Kernel;

  // Repeat Polyhedron types
  typedef typename boost::graph_traits<Polyhedron>::vertex_iterator		  vertex_iterator;
  typedef typename boost::graph_traits<Polyhedron>::edge_iterator		    edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::in_edge_iterator		in_edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::out_edge_iterator		out_edge_iterator;
  
  typedef internal::Spokes_and_rims_iterator<Polyhedron> Rims_iterator;

  // Handle container types
  typedef std::vector<vertex_descriptor>            Handle_container;
  typedef std::list<Handle_container>               Handle_group_container;
public:
  /** The type for returned handle group representative from insert_handle(vertex_descriptor vd), insert_handle(InputIterator begin, InputIterator end) */
  typedef typename Handle_group_container::iterator Handle_group;  
                                                                      
// Data members.
public:
  Polyhedron& polyhedron;															/**< Source triangulated surface mesh for modeling */

private:
  std::vector<Point> original;                        // original positions of roi (size: ros + boundary_of_ros)
  std::vector<Point> solution;                        // storing position of ros vertices during iterations (size: ros + boundary_of_ros)

  VertexIndexMap vertex_index_map;										// storing indices of ros vertices
  EdgeIndexMap   edge_index_map;										  // storing indices of ros related edges

  std::vector<vertex_descriptor> ros;									// region of solution, including roi and hard constraints on boundary of roi
  std::vector<bool> is_roi;                           // (size: ros)
  std::vector<bool> is_hdl;                           // (size: ros)

  std::vector<double> edge_weight;                    // weight of edges only those who are incident to ros 
  std::vector<Eigen::Matrix3d> rot_mtr;               // rotation matrices of ros vertices (size: ros)

  SparseLinearAlgebraTraits_d m_solver;               // linear sparse solver
  unsigned int iterations;                            // number of maximal iterations
  double tolerance;                                   // tolerance of convergence 

  bool need_preprocess;                               // is there any need to call preprocess() function
  Handle_group_container handle_groups;               // user specified handles

  WeightCalculator weight_calculator;                 // calculate weight for an edge

// Public methods
public:
  /**
   * The constructor for deformation object
   *
   * @pre @a polyhedron.is_pure_triangle()
   * @param polyhedron a triangulated surface mesh for modeling
   * @param vertex_index_map_ vertex index map for associating ids with region of interest vertices
   * @param edge_index_map_  edge index map for associating ids with region of interest edges
   * @param iterations number of iterations for each call to deform()
   * @param tolerance ...
   */
  Deform_mesh(Polyhedron& polyhedron, 
              VertexIndexMap vertex_index_map, 
              EdgeIndexMap edge_index_map,
              unsigned int iterations = 5,
              double tolerance = 0.0)
    : polyhedron(polyhedron), vertex_index_map(vertex_index_map), edge_index_map(edge_index_map),
      iterations(iterations), tolerance(tolerance), need_preprocess(true), weight_calculator(polyhedron)
  {
    CGAL_precondition(polyhedron.is_pure_triangle());   
  }

///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Vertex insertion deletion //////////////////////////////////

  /**
   * Clear the internal state of the object, after cleanup the object can be treated as if it is just constructed
   */
  void clear()
  {
    need_preprocess = true;
    // clear vertices
    //roi.clear(); 
    ros.clear();
    handle_groups.clear();
    // no need to clear vertex index map (or edge) since they are going to be reassigned 
    // (at least the useful parts will be reassigned)
  }

  /**
   * Create a new empty handle group for inserting handles.
   * insert_handle(Handle_group handle_group, vertex_descriptor vd) or insert_handle(Handle_group handle_group, InputIterator begin, InputIterator end)
   * can be used for populating a group.
   * After inserting vertices, one can use translate(Handle_group handle_group, const Vector& translation) or rotate(...)
   * to apply transformations on all vertices inside the group. 
   * @return created handle group representative (returned representative is valid until erase_handle(Handle_group handle_group) is called [or copy constructor what to do about it?])
   * @see insert_handle(vertex_descriptor vd), insert_handle(InputIterator begin, InputIterator end)
   */
  Handle_group create_handle_group()
  {
    need_preprocess = true;
    handle_groups.push_back(Handle_container());
    return --handle_groups.end();
  }

  /**
   * -- I think this function can be removed, since combination of other functions simply accomplish the same task.
   \code
    Handle_group handle_group = create_handle_group();
    insert_handle(handle_group, vd);
    // or 
    Handle_group handle_group = insert_handle(vd);
    \endcode
   * Create a new empty handle group and insert vd in it.
   * @param vd vertex to be inserted
   * @return created handle group representative
   * @see insert_handle(Handle_group handle_group, vertex_descriptor vd), 
   * insert_handle(Handle_group handle_group, InputIterator begin, InputIterator end)
   * for inserting more vertices into a handle group
   */ 
  Handle_group insert_handle(vertex_descriptor vd)
  {
    need_preprocess = true;
    Handle_group handle_group = create_handle_group();

    insert_handle(handle_group, vd);
    return handle_group;
  }
  
  /**
   * Insert a vertex into a handle group
   * @param handle_group group to be inserted into
   * @param vd vertex to be inserted
   */
  void insert_handle(Handle_group handle_group, vertex_descriptor vd)
  {
    need_preprocess = true;
    handle_group->push_back(vd);
  }

  /**
   * -- I think this function can be removed, since combination of other functions simply accomplish the same task.
   \code
    Handle_group handle_group = create_handle_group();
    insert_handle(handle_group, begin, end);
    // or 
    Handle_group handle_group = insert_handle(begin, end);
    \endcode
   * Create a new handle group and insert vertices in the range.

   * @tparam InputIterator input iterator type which points to vertex descriptors
   * @param begin iterators spesifying the range of vertices i.e. [begin, end) 
   * @param end iterators spesifying the range of vertices i.e. [begin, end) 
   * It simply corresponds to: 
   */ 
  template<class InputIterator>
  Handle_group insert_handle(InputIterator begin, InputIterator end)
  {
    need_preprocess = true;
    Handle_group handle_group = create_handle_group();

    insert_handle(handle_group, begin, end);
    return handle_group; 
  }

  /**
   * Insert vertices in the range to provided handle group
   * @tparam InputIterator input iterator type which points to vertex descriptors
   * @param handle_group group to be inserted in
   * @param begin iterators spesifying the range of vertices [begin, end) 
   * @param end iterators spesifying the range of vertices [begin, end) 
   */
  template<class InputIterator>
  void insert_handle(Handle_group handle_group, InputIterator begin, InputIterator end)
  {
    need_preprocess = true;
    for( ;begin != end; ++begin)
    {
      insert_handle(handle_group, *begin);
    }
  }

  /**
   * Erase handle group, and invalidate the representative so that it should not be used anymore.
   * @param handle_group group to be erased
   */
  void erase_handle(Handle_group handle_group)
  {
    need_preprocess = true;
    handle_groups.erase(handle_group);
  }

  /**
   * Erase a vertex from a handle group, note that handle group is not erased even if it becomes empty.
   * @param handle_group group to be erased from
   * @param vd vertex to be erased
   */
  void erase_handle(Handle_group handle_group, vertex_descriptor vd)
  {
    need_preprocess = true;
    typename Handle_container::iterator it 
      = std::find(handle_group->begin(), handle_group->end(), vd);
    if(vd != handle_group->end())
    {
      handle_group->erase(it);
      // Although the handle group might get empty, we do not delete it from handle_groups
    }
  }

  /**
   * Insert vertices in the range to region of interest
   * @tparam InputIterator input iterator type which points to vertex descriptors
   * @param begin iterators spesifying the range of vertices [begin, end) 
   * @param end iterators spesifying the range of vertices [begin, end) 
   */
  template<class InputIterator>
  void insert_roi(InputIterator begin, InputIterator end)
  {
    need_preprocess = true;
    for( ;begin != end; ++begin)
    {
      insert_roi(*begin);
    }
  }

  /**
   * Insert a vertex to region of interest
   * @param vd vertex to be inserted
   */
  void insert_roi(vertex_descriptor vd)   
  {
    need_preprocess = true;
    ros.push_back(vd);
  }

  /**
   * Erease a vertex from ROI
   * @param vd vertex to be erased
   */
  void erase_roi(vertex_descriptor vd)   
  {
    need_preprocess = true;
    typename std::vector<vertex_descriptor>::iterator it = std::find(ros.begin(), ros.end(), vd);
    if(vd != ros.end())
    {
      ros.erase(it);
    }
  }
//////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Other utilities ///////////////////////////////////////////

  /**
   * Set the number of iterations used in deform()
   */
  void set_iterations(unsigned int iterations)
  {
    this->iterations = iterations;
  }

  /**
   * Set the tolerance of convergence used in deform()
   * Set to zero if energy based termination is not required.
   */
  void set_tolerance(double tolerance)
  {
    this->tolerance = tolerance;
  }
  
  /**
   * Translate the handle group by translation,
   * in other words every handle vertex in the handle_group is translated from its original position
   * @param handle_group representative of the handle group which is subject to translation
   * @param translation translation vector 
   */
  void translate(Handle_group handle_group, const Vector& translation)
  {
    for(typename Handle_container::iterator it = handle_group->begin();
      it != handle_group->end(); ++it)
    {
      size_t v_index = get(vertex_index_map, *it);
      solution[v_index] = solution[v_index] + translation;
    }
  }

  /**
   * Rotate the handle group around rotation center by quaternion then translate it by translation
   * @tparam Quaternion quaternion type which defines a multiplication operator with Vect as quad * vector
   * @tparam Vect vector type 3 param constructable and has operator[] ...
   */
  template <typename Quaternion, typename Vect>
  void rotate(Handle_group handle_group, const Point& rotation_center, const Quaternion& quat, const Vect& translation)
  {
    for(typename Handle_container::iterator it = handle_group->begin();
      it != handle_group->end(); ++it)
    {
      size_t v_index = get(vertex_index_map, *it);

      Point p = CGAL::ORIGIN + ( original[v_index] - rotation_center);
      Vect v = quat * Vect(p.x(),p.y(),p.z());
      p = Point(v[0], v[1], v[2]) + ( rotation_center - CGAL::ORIGIN);
      p = p + Vector(translation[0],translation[1],translation[2]);

      solution[v_index] = p;
    }
  }
  
  /**
   * Assign positions in the range as target positions for the vertices in the handle group
   * @tparam InputIterator input iterator type which points to Polyhedron::Traits::Point_3
   * @param handle_group group of target vertices
   * @param begin iterators spesifying the range of positions [begin, end) 
   * @param end iterators spesifying the range of positions [begin, end) 
   */
  template<class InputIterator>
  void assign(Handle_group handle_group, InputIterator begin, InputIterator end)
  {
    for(typename Handle_container::iterator it = handle_group->begin(); 
        (it != handle_group->end()) && (begin != end); 
        ++it, ++begin)
    {
      size_t v_index = get(vertex_index_map, *it);
      solution[v_index] = *begin;
    }
  }

  /**
   * Assign the target position for the handle vertes 
   * @param vd handle vertex to be assigned target position
   * @param target_position constrained position
   */
  void assign(vertex_descriptor vd, const Point& target_position)
  {
    size_t v_index = get(vertex_index_map, vd);
    solution[v_index] = target_position;
  }
  /**
   * Reset position of deformed vertices to their original positions (i.e. positions at the time of last preprocess() call)
   */
  void undo()
  {
    for(size_t i = 0; i < ros.size(); ++i)
    {
      ros[i]->point() = original[i];
    }
  }

///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Deformation Core ///////////////////////////////////////////

  /** 
   * Necessary precomputation work before beginning deformation.
   * It needs to be called after insertion of vertices as handles or roi is done.
   * @return true if Laplacian matrix factorization is successful
   */
  bool preprocess()
  {
    need_preprocess = false;

    region_of_solution();
    compute_edge_weight(); // compute_edge_weight() has to come later then region_of_solution()

    // Assemble linear system A*X=B
    typename SparseLinearAlgebraTraits_d::Matrix A(ros.size()); // matrix is definite positive, and not necessarily symmetric
    assemble_laplacian(A);		

    // Pre-factorizing the linear system A*X=B
    double D;
    return m_solver.pre_factor(A, D);
  }

  /**
   * Deformation on roi vertices. Default iteration and tolerance values are used.
   * @see set_iterations(unsigned int iterations), set_tolerance(double tolerance), deform(unsigned int iterations, double tolerance)
   */
  void deform()
  {
    deform(iterations, tolerance);
  }

  /**
   * Deformation on roi vertices, 
   */
  void deform(unsigned int iterations, double tolerance)
  {
    CGAL_precondition(!need_preprocess || !"preprocess() need to be called before deforming!");

    // Note: no energy based termination occurs at first iteration
    // because comparing energy of original model (before deformation) and deformed model (deformed_1_iteration)
    // simply does not make sense, comparison is meaningful between deformed_(i)_iteration & deformed_(i+1)_iteration

    double energy_this = 0; // initial value is not important, because we skip first iteration
    double energy_last;

    // iterations
    for ( unsigned int ite = 0; ite < iterations; ++ite)
    {
      // main steps of optimization
      optimal_rotations_svd(); 
      update_solution();       

      // energy based termination
      if(tolerance > 0.0 && (ite + 1) < iterations) // if tolerance <= 0 then don't compute energy 
      {                                             // also no need compute energy if this iteration is the last iteration
        energy_last = energy_this;
        energy_this = energy();
        if(energy_this < 0)
        {
          // std::cout << "Negative energy" << std::endl;
        }

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

private:

  /// Compute cotangent weights of all edges 
  void compute_edge_weight()
  {
  #if defined(CGAL_DEFORM_SPOKES_AND_RIMS_2)
    compute_edge_weight_spokes_and_rims(); // nothing special
  #elif defined(CGAL_DEFORM_SPOKES_AND_RIMS)
    compute_edge_weight_spokes_and_rims();
  #else
    compute_edge_weight_arap();
  #endif
  }
  void compute_edge_weight_arap()
  {
    std::set<edge_descriptor> have_id; // edges which has assigned ids (and also weights are calculated)

    // iterate over ros vertices and calculate weights for edges which are incident to ros
    size_t next_edge_id = 0;
    for (std::size_t i = 0; i < ros.size(); i++)
    {
      vertex_descriptor vi = ros[i];
      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
      {
        typename std::set<edge_descriptor>::iterator it = have_id.find(*e);
        if(it != have_id.end()) { continue; } // we have assigned an id already, which means we also calculted the weight
        
        put(edge_index_map, *e, next_edge_id++);
        have_id.insert(*e);

        double weight = weight_calculator(*e);
        edge_weight.push_back(weight);
      }// end of edge loop
    }// end of ros loop
  }
  void compute_edge_weight_spokes_and_rims()
  {
    std::set<edge_descriptor> have_id; // edges which has assigned ids (and also weights are calculated)

    // iterate over ros vertices and calculate weights for edges which are incident to ros
    size_t next_edge_id = 0;
    for (std::size_t i = 0; i < ros.size(); i++)
    {       
      vertex_descriptor vi = ros[i];
      out_edge_iterator e_begin, e_end;
      boost::tie(e_begin, e_end) = boost::out_edges(vi, polyhedron);

      for (Rims_iterator rims_it(e_begin, polyhedron); rims_it.get_iterator() != e_end; ++rims_it )
      {
        edge_descriptor active_edge = rims_it.get_descriptor();

        typename std::set<edge_descriptor>::iterator it = have_id.find(active_edge);
        if(it == have_id.end()) // we have not assigned an id yet
        {  
          put(edge_index_map, active_edge, next_edge_id++);
          have_id.insert(active_edge);
          double weight = weight_calculator(active_edge);
          edge_weight.push_back(weight);

          edge_descriptor opp = CGAL::opposite_edge(active_edge, polyhedron);

          put(edge_index_map, opp, next_edge_id++);
          have_id.insert(opp);
          weight = weight_calculator(opp);
          edge_weight.push_back(weight);
        }
      }// end of edge loop
    }// end of ros loop
  }

  /// Assigns id to one rign neighbor of vd, and also push them into push_vector
  void assign_id_to_one_ring(vertex_descriptor vd, 
                             std::size_t& next_id, 
                             std::vector<vertex_descriptor>& push_vector,
                             std::set<vertex_descriptor>& have_id)
  {
    in_edge_iterator e, e_end;
    for (boost::tie(e,e_end) = boost::in_edges(vd, polyhedron); e != e_end; e++)
    {
      vertex_descriptor vt = boost::source(*e, polyhedron);
      typename std::set<vertex_descriptor>::iterator it = have_id.find(vt);
      if( it == have_id.end() )  // neighboring vertex which is outside of roi and not visited previously (i.e. need an id)
      {
        put(vertex_index_map, vt, next_id++);
        have_id.insert(vt);
        push_vector.push_back(vt);        
      }
    }
  }
  
  /// Find region of solution, including roi and hard constraints, which is the 1-ring vertices out roi
  void region_of_solution()
  {
    // Important: at this point ros contains the roi vertices only.
    // copy roi vertices to roi vector
    std::vector<vertex_descriptor> roi;	// we can remove this temp, but I keep to simplify things below
    roi.insert(roi.end(), ros.begin(), ros.end()); 

    ////////////////////////////////////////////////////////////////
    // assign id to vertices inside: roi, boundary of roi (roi + boundary of roi = ros),
    //                               and boundary of ros
    std::set<vertex_descriptor> have_id;         // keep vertices which are assigned an id
    
    for(std::size_t i = 0; i < roi.size(); i++)  // assign id to all roi vertices
    {
      put(vertex_index_map, roi[i], i);
    }

    have_id.insert(roi.begin(), roi.end());      // mark roi vertices since they have ids now

    // now assign an id to vertices on boundary of roi
    std::size_t next_ros_index = roi.size();
    for(std::size_t i = 0; i < roi.size(); i++)
    {
      assign_id_to_one_ring(roi[i], next_ros_index, ros, have_id);
    }

    std::vector<vertex_descriptor> outside_ros;
    // boundary of ros also must have ids because in SVD calculation,
    // one-ring neighbor of ROS vertices are reached. 
    for(std::size_t i = roi.size(); i < ros.size(); i++)
    {
      assign_id_to_one_ring(ros[i], next_ros_index, outside_ros, have_id);
    }
    ////////////////////////////////////////////////////////////////

    // initialize the rotation matrices (size: ros)
    rot_mtr.resize(ros.size());
    for(std::size_t i = 0; i < rot_mtr.size(); i++)
    {
      rot_mtr[i].setIdentity();
    }
    
    // initialize solution and original (size: ros + boundary_of_ros)

    // for simplifying coding effort, I also put boundary of ros into solution and original
    // because boundary of ros vertices are reached in optimal_rotations_svd() and energy()
    solution.resize(ros.size() + outside_ros.size());
    original.resize(ros.size() + outside_ros.size());
    
    for(std::size_t i = 0; i < ros.size(); i++)
    {
      solution[i] = ros[i]->point();
      original[i] = ros[i]->point();
    }
    for(std::size_t i = 0; i < outside_ros.size(); ++i)
    {
      original[ros.size() + i] = outside_ros[i]->point();
      solution[ros.size() + i] = outside_ros[i]->point();
    }

    // initialize flag vectors of roi, handle, ros 
    is_roi.assign(ros.size(), false);
    is_hdl.assign(ros.size(), false);
    for(std::size_t i = 0; i < roi.size(); i++)
    {
      size_t v_index = get(vertex_index_map, roi[i]);
      is_roi[v_index] = true;
    }
    for(typename Handle_group_container::iterator it_group = handle_groups.begin(); 
      it_group != handle_groups.end(); ++it_group)
    {
      for(typename Handle_container::iterator it_vertex = it_group->begin(); 
        it_vertex != it_group->end(); ++it_vertex)
      {
        size_t v_index = get(vertex_index_map, *it_vertex);
        is_hdl[v_index] = true;
      }
    }
  }

  /// Assemble Laplacian matrix A of linear system A*X=B
  void assemble_laplacian(typename SparseLinearAlgebraTraits_d::Matrix& A)
  {
  #if defined(CGAL_DEFORM_SPOKES_AND_RIMS_2)
    assemble_laplacian_spokes_and_rims_2(A);
  #elif defined(CGAL_DEFORM_SPOKES_AND_RIMS)
    assemble_laplacian_spokes_and_rims(A);
  #else
    assemble_laplacian_arap(A);
  #endif
  }
  void assemble_laplacian_arap(typename SparseLinearAlgebraTraits_d::Matrix& A)
  {
    /// assign cotangent Laplacian to ros vertices
    for(std::size_t i = 0; i < ros.size(); i++)
    {
      vertex_descriptor vi = ros[i];
      std::size_t vertex_idx_i = get(vertex_index_map, vi);
      if ( is_roi[vertex_idx_i] && !is_hdl[vertex_idx_i] )          // vertices of ( roi - hdl )
      {
        double diagonal = 0;
        in_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
        {
          vertex_descriptor vj = boost::source(*e, polyhedron);
          double wij = edge_weight[ get(edge_index_map, *e)];  // edge weights
          double wji = edge_weight[get(edge_index_map, CGAL::opposite_edge(*e, polyhedron))];
          double total_weight = wij + wji;

          std::size_t vj_index = get(vertex_index_map, vj);
          A.set_coef(i, vj_index, -total_weight, true);	// off-diagonal coefficient
          diagonal += total_weight;  
        }
        // diagonal coefficient
        A.set_coef(i, i, diagonal, true);
      }
      else
      {
        A.set_coef(i, i, 1.0, true);
      }
    }
  }
  void assemble_laplacian_spokes_and_rims(typename SparseLinearAlgebraTraits_d::Matrix& A)
  {
     /// assign cotangent Laplacian to ros vertices
    for(std::size_t i = 0; i < ros.size(); i++)
    {
      vertex_descriptor vi = ros[i];
      std::size_t vertex_idx_i = get(vertex_index_map, vi);
      if ( is_roi[vertex_idx_i] && !is_hdl[vertex_idx_i] )          // vertices of ( roi - hdl )
      {
        double diagonal = 0;
        in_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
        {
          vertex_descriptor vj = boost::source(*e, polyhedron);
          double wij = edge_weight[ get(edge_index_map, *e)];  // edge weights
          double wji = edge_weight[get(edge_index_map, CGAL::opposite_edge(*e, polyhedron))];

          double total_weight = (wij + wji);

          if (boost::get(CGAL::edge_is_border, polyhedron, *e) ||
            boost::get(CGAL::edge_is_border, polyhedron, CGAL::opposite_edge(*e, polyhedron)))
          {
              total_weight *= 1.5;
          }
          else
          {
            total_weight *= 2;
          }

          std::size_t vj_index = get(vertex_index_map, vj);
          A.set_coef(i, vj_index, -total_weight, true);	// off-diagonal coefficient
          diagonal += total_weight;  
        }
        // diagonal coefficient
        A.set_coef(i, i, diagonal, true);
      }
      else
        A.set_coef(i, i, 1.0, true);
    }
  }
  void assemble_laplacian_spokes_and_rims_2(typename SparseLinearAlgebraTraits_d::Matrix& A)
  {
    /// assign cotangent Laplacian to ros vertices
    for(std::size_t i = 0; i < ros.size(); i++)
    {
      vertex_descriptor vi = ros[i];
      std::size_t vertex_idx_i = get(vertex_index_map, vi);
      if ( is_roi[vertex_idx_i] && !is_hdl[vertex_idx_i] )          // vertices of ( roi - hdl )
      {
        double diagonal = 0;
        out_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::out_edges(vi, polyhedron); e != e_end; e++)
        {
          if(boost::get(CGAL::edge_is_border, polyhedron, *e)) { continue; } // no facet 

          vertex_descriptor vj = boost::target(*e, polyhedron);
          double wij = edge_weight[ get(edge_index_map, *e)];  // edge weights
          double total_weight = (wij * 3); 

          edge_descriptor opp = CGAL::opposite_edge(*e, polyhedron);
          if(boost::get(CGAL::edge_is_border, polyhedron, opp)) { continue; } // no facet 

          double wji = edge_weight[ get(edge_index_map, opp)];  // edge weights

          total_weight += (wji * 3);

          std::size_t vj_index = get(vertex_index_map, vj);
          A.set_coef(i, vj_index, -total_weight, true);	// off-diagonal coefficient
          diagonal += total_weight; 
        }
        // diagonal coefficient
        A.set_coef(i, i, diagonal, true);
      }
      else
      {
        A.set_coef(i, i, 1.0, true);
      }
    }
  }

  /// Local step of iterations, computing optimal rotation matrices using SVD decomposition, stable
  void optimal_rotations_svd()
  {
  #if defined(CGAL_DEFORM_SPOKES_AND_RIMS_2)
    optimal_rotations_svd_spokes_and_rims_2();
  #elif defined(CGAL_DEFORM_SPOKES_AND_RIMS)
    optimal_rotations_svd_spokes_and_rims();
  #else
    optimal_rotations_svd_arap();
  #endif
  }
  void optimal_rotations_svd_arap()
  {
    Eigen::Matrix3d u, v;           // orthogonal matrices 
    Eigen::Vector3d w;              // singular values
    Eigen::Matrix3d cov;            // covariance matrix
    Eigen::JacobiSVD<Eigen::Matrix3d> svd;       // SVD solver         
    Eigen::Matrix3d r;

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
        size_t vj_index = get(vertex_index_map, vj);

        Vector pij = original[i] - original[vj_index];
        Vector qij = solution[i] - solution[vj_index];

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
      svd.compute( cov, Eigen::ComputeFullU | Eigen::ComputeFullV );
      u = svd.matrixU(); v = svd.matrixV();

      // extract rotation matrix
      r = v*u.transpose();

      // checking negative determinant of r
      if ( r.determinant() < 0 )    // changing the sign of column corresponding to smallest singular value
      {
        u.col(2) *= -1;      // singular values are always sorted in decresing order so use column 2
        r = v*u.transpose(); // re-extract rotation matrix
      }
      
      rot_mtr[i] = r;
    }
  }
  void optimal_rotations_svd_spokes_and_rims()
  {
    Eigen::Matrix3d u, v;           // orthogonal matrices 
    Eigen::Vector3d w;              // singular values
    Eigen::Matrix3d cov;            // covariance matrix
    Eigen::JacobiSVD<Eigen::Matrix3d> svd;       // SVD solver         
    Eigen::Matrix3d r;

    // only accumulate ros vertices
    for ( std::size_t i = 0; i < ros.size(); i++ )
    {
      vertex_descriptor vi = ros[i];
      // compute covariance matrix
      cov.setZero();

      // spoke + rim edges
      out_edge_iterator e_begin, e_end;
      boost::tie(e_begin, e_end) = boost::out_edges(vi, polyhedron);

      for (Rims_iterator rims_it(e_begin, polyhedron); rims_it.get_iterator() != e_end; ++rims_it )
      { 
        edge_descriptor active_edge = rims_it.get_descriptor();

        vertex_descriptor v1 = boost::source(active_edge, polyhedron);
        vertex_descriptor v2 = boost::target(active_edge, polyhedron);
        
        size_t v1_index = get(vertex_index_map, v1);
        size_t v2_index = get(vertex_index_map, v2);
        
        Vector pij = original[v1_index] - original[v2_index];
        Vector qij = solution[v1_index] - solution[v2_index];
        
        std::size_t edge_id = get(edge_index_map, active_edge);
        double wij = edge_weight[edge_id];
        for (int j = 0; j < 3; j++)
        {
          for (int k = 0; k < 3; k++)
          {
            cov(j, k) += wij*pij[j]*qij[k]; 
          }
        }
      }
  
      // svd decomposition
      svd.compute( cov, Eigen::ComputeFullU | Eigen::ComputeFullV );
      u = svd.matrixU(); v = svd.matrixV();

      // extract rotation matrix
      r = v*u.transpose();

      // checking negative determinant of r
      if ( r.determinant() < 0 )    // changing the sign of column corresponding to smallest singular value
      {
        u.col(2) *= -1;      // singular values are always sorted in decresing order so use column 2
        r = v*u.transpose(); // re-extract rotation matrix
      }
      
      rot_mtr[i] = r;
    }
  }
  void optimal_rotations_svd_spokes_and_rims_2()
  {
    Eigen::Matrix3d u, v;           // orthogonal matrices 
    Eigen::Vector3d w;              // singular values
    Eigen::Matrix3d cov;            // covariance matrix
    Eigen::JacobiSVD<Eigen::Matrix3d> svd;       // SVD solver         
    Eigen::Matrix3d r;

    // only accumulate ros vertices
    for ( std::size_t i = 0; i < ros.size(); i++ )
    {
      vertex_descriptor vi = ros[i];
      // compute covariance matrix
      cov.setZero();

      //iterate through all triangles 
      out_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::out_edges(vi, polyhedron); e != e_end; e++)
      {
        if(boost::get(CGAL::edge_is_border, polyhedron, *e)) { continue; } // no facet 
        // iterate edges around facet
        edge_descriptor edge_around_facet = *e;
        do
        {
          vertex_descriptor v1 = boost::source(edge_around_facet, polyhedron);
          vertex_descriptor v2 = boost::target(edge_around_facet, polyhedron);

          size_t v1_index = get(vertex_index_map, v1);
          size_t v2_index = get(vertex_index_map, v2);
        
          Vector pij = original[v1_index] - original[v2_index];
          Vector qij = solution[v1_index] - solution[v2_index];

          std::size_t edge_id = get(edge_index_map, edge_around_facet);
          double wij = edge_weight[edge_id];
          for (int j = 0; j < 3; j++)
          {
            for (int k = 0; k < 3; k++)
            {
              cov(j, k) += wij*pij[j]*qij[k]; 
            }
          }
        } while( (edge_around_facet = CGAL::next_edge(edge_around_facet, polyhedron)) != *e);
      }
  
      // svd decomposition
      svd.compute( cov, Eigen::ComputeFullU | Eigen::ComputeFullV );
      u = svd.matrixU(); v = svd.matrixV();

      // extract rotation matrix
      r = v*u.transpose();

      // checking negative determinant of r
      if ( r.determinant() < 0 )    // changing the sign of column corresponding to smallest singular value
      {
        u.col(2) *= -1;      // singular values are always sorted in decresing order so use column 2
        r = v*u.transpose(); // re-extract rotation matrix
      }
      
      rot_mtr[i] = r;
    }
  }

  /// Global step of iterations, updating solution
  void update_solution()
  {
  #if defined(CGAL_DEFORM_SPOKES_AND_RIMS_2)
    update_solution_spokes_and_rims_2();
  #elif defined(CGAL_DEFORM_SPOKES_AND_RIMS)
    update_solution_spokes_and_rims();
  #else
    update_solution_arap();
  #endif
  }
  void update_solution_arap()
  {
    typename SparseLinearAlgebraTraits_d::Vector X(ros.size()), Bx(ros.size());
    typename SparseLinearAlgebraTraits_d::Vector Y(ros.size()), By(ros.size());
    typename SparseLinearAlgebraTraits_d::Vector Z(ros.size()), Bz(ros.size());

    // assemble right columns of linear system
    for ( std::size_t i = 0; i < ros.size(); i++ )
    {
      vertex_descriptor vi = ros[i];
      std::size_t vertex_idx_i = get(vertex_index_map, vi);
      if ( !is_roi[vertex_idx_i] || is_hdl[vertex_idx_i] )   // hard constraints or handle vertices
      {
        Bx[i] = solution[vertex_idx_i].x(); By[i] = solution[vertex_idx_i].y(); Bz[i] = solution[vertex_idx_i].z();
      }
      else  // ( roi - handle ) vertices
      {
        Bx[i] = 0; By[i] = 0; Bz[i] = 0;
        in_edge_iterator e, e_end;
        Point& pi = original[get(vertex_index_map, vi)];
        for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
        {
          vertex_descriptor vj = boost::source(*e, polyhedron);
          std::size_t vj_index = get(vertex_index_map, vj) ; 
          Vector pij =  pi - original[get(vertex_index_map, vj) ];
          double wij = edge_weight[get(edge_index_map, *e) ];
          double wji = edge_weight[get(edge_index_map, CGAL::opposite_edge(*e, polyhedron))];

          double x, y, z;
          x = y = z = 0.0;
          for (int j = 0; j < 3; j++)
          {
            x += ( rot_mtr[i](0, j)*wij + rot_mtr[vj_index](0, j)*wji ) * pij[j];
            y += ( rot_mtr[i](1, j)*wij + rot_mtr[vj_index](1, j)*wji ) * pij[j];
            z += ( rot_mtr[i](2, j)*wij + rot_mtr[vj_index](2, j)*wji ) * pij[j];
          }
          Bx[i] += x; By[i] += y; Bz[i] += z; 
        }
      }
    }

    // solve "A*X = B".
    m_solver.linear_solver(Bx, X); m_solver.linear_solver(By, Y); m_solver.linear_solver(Bz, Z);

    // copy to solution
    for (std::size_t i = 0; i < ros.size(); i++)
    {
      Point p(X[i], Y[i], Z[i]);
      solution[get(vertex_index_map, ros[i])] = p;
    }

  }
  void update_solution_spokes_and_rims()
  {
    typename SparseLinearAlgebraTraits_d::Vector X(ros.size()), Bx(ros.size());
    typename SparseLinearAlgebraTraits_d::Vector Y(ros.size()), By(ros.size());
    typename SparseLinearAlgebraTraits_d::Vector Z(ros.size()), Bz(ros.size());

    // assemble right columns of linear system
    for ( std::size_t i = 0; i < ros.size(); i++ )
    {
      vertex_descriptor vi = ros[i];
      std::size_t vertex_idx_i = get(vertex_index_map, vi);
      if ( !is_roi[vertex_idx_i] || is_hdl[vertex_idx_i] )   // hard constraints or handle vertices
      {
        Bx[i] = solution[vertex_idx_i].x(); By[i] = solution[vertex_idx_i].y(); Bz[i] = solution[vertex_idx_i].z();
      }
      else  // ( roi - handle ) vertices
      {
        Bx[i] = 0; By[i] = 0; Bz[i] = 0;
        in_edge_iterator e, e_end;
        Point& pi = original[get(vertex_index_map, vi)];
        for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
        {
          vertex_descriptor vj = boost::source(*e, polyhedron);
          std::size_t vj_index = get(vertex_index_map, vj); 
          Vector pij =  pi - original[get(vertex_index_map, vj)];

          double wij = edge_weight[get(edge_index_map, *e)];
          double wji = edge_weight[get(edge_index_map, CGAL::opposite_edge(*e, polyhedron))];

          bool is_e_border = boost::get(CGAL::edge_is_border, polyhedron, *e);
          bool is_opp_of_e_border = boost::get(CGAL::edge_is_border, polyhedron, CGAL::opposite_edge(*e, polyhedron));
          double x, y, z;
          x = y = z = 0.0;
          if (is_e_border || is_opp_of_e_border)
          {  
            edge_descriptor next_of_not_border_edge = is_e_border ? CGAL::next_edge(CGAL::opposite_edge(*e, polyhedron), polyhedron)
              : CGAL::next_edge(*e, polyhedron);
            std::size_t vk_index = get(vertex_index_map, boost::target(next_of_not_border_edge, polyhedron));

            double wk = is_e_border ? wji : wij;

            for (int j = 0; j < 3; j++)
            {
              x += ( rot_mtr[i](0, j)*wij + rot_mtr[vj_index](0, j)*wji +
                     rot_mtr[vk_index](0, j)*wk  ) * pij[j];

              y += ( rot_mtr[i](1, j)*wij + rot_mtr[vj_index](1, j)*wji +
                     rot_mtr[vk_index](1, j)*wk ) * pij[j];

              z += ( rot_mtr[i](2, j)*wij + rot_mtr[vj_index](2, j)*wji +
                     rot_mtr[vk_index](2, j)*wk ) * pij[j];
            }
          }
          else
          {
            edge_descriptor next_e = CGAL::next_edge(*e, polyhedron);
            vertex_descriptor next_e_target = boost::target(next_e, polyhedron);
            std::size_t vk_index = get(vertex_index_map, next_e_target) ;

            edge_descriptor next_opp_e = CGAL::next_edge(CGAL::opposite_edge(*e, polyhedron), polyhedron);
            vertex_descriptor next_opp_e_target = boost::target(next_opp_e, polyhedron);
            std::size_t vm_index = get(vertex_index_map, next_opp_e_target) ;

            for (int j = 0; j < 3; j++)
            {
              x += ( rot_mtr[i](0, j)*wij + rot_mtr[vj_index](0, j)*wji +
                     rot_mtr[vk_index](0, j)*wij + rot_mtr[vm_index](0, j)*wji) * pij[j];

              y += ( rot_mtr[i](1, j)*wij + rot_mtr[vj_index](1, j)*wji +
                     rot_mtr[vk_index](1, j)*wij + rot_mtr[vm_index](1, j)*wji) * pij[j];

              z += ( rot_mtr[i](2, j)*wij + rot_mtr[vj_index](2, j)*wji +
                     rot_mtr[vk_index](2, j)*wij + rot_mtr[vm_index](2, j)*wji) * pij[j];
            }
          }

          Bx[i] += x; By[i] += y; Bz[i] += z; 
        }
      }
    }
    // solve "A*X = B".
    m_solver.linear_solver(Bx, X); m_solver.linear_solver(By, Y); m_solver.linear_solver(Bz, Z);

    // copy to solution
    for (std::size_t i = 0; i < ros.size(); i++)
    {
      Point p(X[i], Y[i], Z[i]);
      solution[get(vertex_index_map, ros[i])] = p;
    }
  }
  void update_solution_spokes_and_rims_2()
  {
    typename SparseLinearAlgebraTraits_d::Vector X(ros.size()), Bx(ros.size());
    typename SparseLinearAlgebraTraits_d::Vector Y(ros.size()), By(ros.size());
    typename SparseLinearAlgebraTraits_d::Vector Z(ros.size()), Bz(ros.size());

    // assemble right columns of linear system
    for ( std::size_t i = 0; i < ros.size(); i++ )
    {
      vertex_descriptor vi = ros[i];
      std::size_t vertex_idx_i = get(vertex_index_map, vi);
      if ( !is_roi[vertex_idx_i] || is_hdl[vertex_idx_i] )   // hard constraints or handle vertices
      {
        Bx[i] = solution[vertex_idx_i].x(); By[i] = solution[vertex_idx_i].y(); Bz[i] = solution[vertex_idx_i].z();
      }
      else  // ( roi - handle ) vertices
      {
        Bx[i] = 0; By[i] = 0; Bz[i] = 0;
        Point& pi = original[get(vertex_index_map, vi)];

        out_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::out_edges(vi, polyhedron); e != e_end; e++)
        {
          if(boost::get(CGAL::edge_is_border, polyhedron, *e)) { continue; } // no facet 
          // iterate edges around facet
          vertex_descriptor vj = boost::target(*e, polyhedron);
          std::size_t vj_index = get(vertex_index_map, vj) ; 

          vertex_descriptor vk = boost::target(CGAL::next_edge(*e, polyhedron), polyhedron);
          std::size_t vk_index = get(vertex_index_map, vk);

          Vector pij =  pi - original[get(vertex_index_map, vj) ];

          double wij = edge_weight[get(edge_index_map, *e)];

          double x, y, z;
          x = y = z = 0.0;
          for (int j = 0; j < 3; j++)
          {
            x += ( rot_mtr[i](0, j)*wij + rot_mtr[vj_index](0, j)*wij +
                    rot_mtr[vk_index](0, j)*wij  ) * pij[j];

            y += ( rot_mtr[i](1, j)*wij + rot_mtr[vj_index](1, j)*wij +
                    rot_mtr[vk_index](1, j)*wij ) * pij[j];

            z += ( rot_mtr[i](2, j)*wij + rot_mtr[vj_index](2, j)*wij +
                    rot_mtr[vk_index](2, j)*wij ) * pij[j];
          }

          edge_descriptor opp = CGAL::opposite_edge(*e, polyhedron);
          if(boost::get(CGAL::edge_is_border, polyhedron, opp)) { continue; } // no facet

          vk = boost::target(CGAL::next_edge(opp, polyhedron), polyhedron);
          vk_index = get(vertex_index_map, vk);

          double wji = edge_weight[ get(edge_index_map, opp)];  // edge weights

          for (int j = 0; j < 3; j++)
          {
            x += ( rot_mtr[i](0, j)*wji + rot_mtr[vj_index](0, j)*wji +
                    rot_mtr[vk_index](0, j)*wji  ) * pij[j];

            y += ( rot_mtr[i](1, j)*wji + rot_mtr[vj_index](1, j)*wji +
                    rot_mtr[vk_index](1, j)*wji ) * pij[j];

            z += ( rot_mtr[i](2, j)*wji + rot_mtr[vj_index](2, j)*wji +
                    rot_mtr[vk_index](2, j)*wji ) * pij[j];
          }

          Bx[i] += x; By[i] += y; Bz[i] += z; 
        }
      }
    }
    // solve "A*X = B".
    m_solver.linear_solver(Bx, X); m_solver.linear_solver(By, Y); m_solver.linear_solver(Bz, Z);

    // copy to solution
    for (std::size_t i = 0; i < ros.size(); i++)
    {
      Point p(X[i], Y[i], Z[i]);
      solution[get(vertex_index_map, ros[i])] = p;
    }
  }


  /// Assign solution to target mesh
  void assign_solution()
  {
    for(std::size_t i = 0; i < ros.size(); ++i){
      std::size_t v_id = get(vertex_index_map, ros[i]);
      if(is_roi[v_id])
      {
        ros[i]->point() = solution[v_id];
      }
    }
  }

  /// Compute modeling energy
  double energy()
  {
  #ifdef CGAL_DEFORM_SPOKES_AND_RIMS
    return energy_spokes_and_rims();
  #else
    return energy_arap();
  #endif
  }
  double energy_arap()
  {
    double sum_of_energy = 0;
    // only accumulate ros vertices
    for( std::size_t i = 0; i < ros.size(); i++ )
    {
      vertex_descriptor vi = ros[i];
      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
      {
        vertex_descriptor vj = boost::source(*e, polyhedron);
        size_t vj_index = get(vertex_index_map, vj);

        Vector pij = original[i] - original[vj_index];
        double wij = edge_weight[get(edge_index_map, *e)];

        Vector rot_p(0, 0, 0);                 // vector rot_i*p_ij
        for (int j = 0; j < 3; j++)
        {
          double x = rot_mtr[i](0, j) * pij[j];
          double y = rot_mtr[i](1, j) * pij[j];
          double z = rot_mtr[i](2, j) * pij[j];
          Vector v(x, y, z);
          rot_p = rot_p + v;
        }
        Vector qij = solution[i] - solution[vj_index];
        sum_of_energy += wij*(qij - rot_p).squared_length();
      }
    }
    return sum_of_energy;
  }
  double energy_spokes_and_rims()
  {
    double sum_of_energy = 0;
    // only accumulate ros vertices
    for( std::size_t i = 0; i < ros.size(); i++ )
    {
      // spoke + rim edges
      out_edge_iterator e_begin, e_end;
      boost::tie(e_begin, e_end) = boost::out_edges(ros[i], polyhedron);

      for (Rims_iterator rims_it(e_begin, polyhedron); rims_it.get_iterator() != e_end; ++rims_it )
      { 
        edge_descriptor active_edge = rims_it.get_descriptor();

        vertex_descriptor v1 = boost::source(active_edge, polyhedron);
        vertex_descriptor v2 = boost::target(active_edge, polyhedron);
        
        size_t v1_index = get(vertex_index_map, v1);
        size_t v2_index = get(vertex_index_map, v2);
        
        Vector pij = original[v1_index] - original[v2_index];
        Vector qij = solution[v1_index] - solution[v2_index];
        double wij = edge_weight[get(edge_index_map, active_edge)];

        Vector rot_p(0, 0, 0);                 // vector rot_i*p_ij
        for (int j = 0; j < 3; j++)
        {
          double x = rot_mtr[i](0, j) * pij[j];
          double y = rot_mtr[i](1, j) * pij[j];
          double z = rot_mtr[i](2, j) * pij[j];
          Vector v(x, y, z);
          rot_p = rot_p + v;
        }
        sum_of_energy += wij*(qij - rot_p).squared_length();
      }
    }
    return sum_of_energy;
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
      // Fallback to an accurate SVD based decomposiiton method.
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
      // Fallback to an accurate SVD based decomposiiton method.
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

