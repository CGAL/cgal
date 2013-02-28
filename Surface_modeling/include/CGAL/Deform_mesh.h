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

#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <CGAL/FPU_extension.h>

#include <Eigen/Eigen>
#include <Eigen/SVD>

// #define CGAL_DEFORM_SPOKES_AND_RIMS

namespace CGAL {

/// \ingroup PkgSurfaceModeling
/**
 * @brief Class providing the functionalities for deforming a triangulated surface mesh
 *
 * @pre @a polyhedron.is_pure_triangle()
 * @tparam Polyhedron a model of HalfedgeGraph
 * @tparam SparseLinearAlgebraTraitsWithPreFactor_d sparse linear solver for square symmetric sparse linear systems
 * @tparam VertexIndexMap a <a href="http://www.boost.org/doc/libs/release/libs/property_map/doc/ReadWritePropertyMap.html">`ReadWritePropertyMap`</a>  with ::vertex_descriptor as key and `unsigned int` as value type
 * @tparam EdgeIndexMap a <a href="http://www.boost.org/doc/libs/release/libs/property_map/doc/ReadWritePropertyMap.html">`ReadWritePropertyMap`</a>  with ::edge_descriptor as key and `unsigned int` as value type
 * @tparam WeightCalculator how to document this (should I provide a concept, like in SegmentationGeomTraits ?)
 */
template <
  class Polyhedron, 
  class SparseLinearAlgebraTraits_d, 
  class VertexIndexMap, 
  class EdgeIndexMap,
  class WeightCalculator = internal::Cotangent_weight<Polyhedron >
  >
class Deform_mesh
{
//Typedefs
public:
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor	vertex_descriptor; /**< The type for vertex representative objects */
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor		edge_descriptor;   /**< The type for edge representative objects */

private:
  // Geometric types              
  typedef typename Polyhedron::Traits         Kernel;
  typedef typename Kernel::Vector_3           Vector;
  typedef typename Kernel::Point_3            Point;

  // Repeat Polyhedron types
  typedef typename boost::graph_traits<Polyhedron>::vertex_iterator		  vertex_iterator;
  typedef typename boost::graph_traits<Polyhedron>::edge_iterator		    edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::in_edge_iterator		in_edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::out_edge_iterator		out_edge_iterator;
  
  // Handle container types
  typedef std::vector<vertex_descriptor>            Handle_container;
  typedef std::list<Handle_container>               Handle_group_container;
public:
  /** The type for returned handle group representative from insert_handle(vertex_descriptor vd), insert_handle(InputIterator begin, InputIterator end) */
  typedef typename Handle_group_container::iterator Handle_group;  
                                                                      
// Data members.
public:
  Polyhedron& polyhedron;															/**< Source triangulated surface mesh for modeling */
  std::vector<Point> original;                        // original positions of roi

private:
  VertexIndexMap vertex_index_map;										// storing indices of ros vertices, others should be 0
  EdgeIndexMap   edge_index_map;										  // storing indices of ros related edges, others should be 0							

  std::vector<vertex_descriptor> roi;                 // region of interest, including both free and haldle vertices
  std::vector<vertex_descriptor> ros;									// region of solution, including roi and hard constraints outside roi
  std::vector<vertex_descriptor> outside_ros;         // boundary of ros, for clearing purpose

  // properties per ros vertex, indexed by vertex_index_map[vertex_descriptor] -1
  std::vector<bool> is_roi;                            
  std::vector<bool> is_hdl;                       
  std::vector<Eigen::Matrix3d> rot_mtr;               // rotation matrices of ros vertices
  std::vector<Point> solution;                        // storing position of ros vertices during iterations

  std::vector<double> edge_weight;                    // weight of edges only those who are incident to ros 

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
   * @param vertex_index_map_ zero initialized vertex index map
   * @param vertex_index_map_ zero initialized edge index map
   * @param iterations number of iterations for each call to deform()
   * @param tolerance ...
   */
  Deform_mesh(Polyhedron& polyhedron, 
              const VertexIndexMap& vertex_index_map_, 
              const EdgeIndexMap& edge_index_map_,
              unsigned int iterations = 5,
              double tolerance = 1e-4)
    : polyhedron(polyhedron), vertex_index_map(vertex_index_map_), edge_index_map(edge_index_map_),
      weight_calculator(polyhedron), need_preprocess(true), iterations(iterations), tolerance(tolerance)
  {
    CGAL_precondition(polyhedron.is_pure_triangle());
    /////////////////////////////////////////////////////////////////
    // this part should be removed since it iterates over all vertices,
    // we can achieve that by a requiring that supplied vertex_index_map, and 
    // edge_index_map should be filled by 0.

    // Q: how can it be different than looping over all vertices ?
    // the user might provide a custom pmap, such as Polyhedron_vertex_zero_default_index_map 
    // (added in demo Property_maps_for_edit_plugin.h) which use a map and returns 0 for not found keys.
    // so no actual initialization takes place.
    vertex_iterator vb, ve;
    for(boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb )
    {
      put(vertex_index_map, *vb, 0);
    }
    // this part should be removed same as above
    edge_iterator eb, ee;
    for(boost::tie(eb,ee) = boost::edges(polyhedron); eb != ee; ++eb )
    {
      put(edge_index_map, *eb, 0);
    }
    /////////////////////////////////////////////////////////////////    
  }

  /**
   * Clear the internal state of the object, after cleanup the object can be treated as if it is just constructed
   */
  void clear()
  {
    need_preprocess = true;
    //clear vertices
    roi.clear(); 
    handle_groups.clear();

    for (std::size_t i = 0; i < ros.size(); i++)
    {
      put(vertex_index_map, ros[i], 0);
    }
    for (std::size_t i = 0; i < outside_ros.size(); i++)
    {
      put(vertex_index_map, outside_ros[i], 0);
    }
    // note that cleaning/reassigning is_roi, is_hdl, ros, solution, original vectors
    // handled in region_of_solution().

    //clear edges
  #ifdef CGAL_DEFORM_SPOKES_AND_RIMS
    edge_weight.clear();
    for (std::size_t i = 0; i < ros.size(); i++)
    {
      bool rim = false;
      out_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::out_edges(ros[i], polyhedron); e != e_end;)
      {
        edge_descriptor active_edge = rim ? CGAL::next_edge(*e, polyhedron): *e;
        put(edge_index_map, active_edge, 0);
        edge_descriptor e_oppo = CGAL::opposite_edge(active_edge, polyhedron);
        put(edge_index_map, e_oppo, 0);
        if(rim) { ++e; }
        rim = !rim;
      }
    }
  #else
    edge_weight.clear();
    for (std::size_t i = 0; i < ros.size(); i++)
    {
      out_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::out_edges(ros[i], polyhedron); e != e_end; e++)
      {
        put(edge_index_map, *e, 0);
        edge_descriptor e_oppo = CGAL::opposite_edge(*e, polyhedron);
        put(edge_index_map, e_oppo, 0);
      }
    }
  #endif
  }

////////////// Handle insertion and deletion //////////////
  /**
   * Create a new empty handle group for inserting handles
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
   * Create a new empty handle group and insert vd in it.
   * @return created handle group representative
   * @see insert_handle(Handle_group handle_group, vertex_descriptor vd), 
   * insert_handle(Handle_group handle_group, InputIterator begin, InputIterator end)
   * for inserting more vertices into a handle group
   */ 
  Handle_group insert_handle(vertex_descriptor vd)
  {
    need_preprocess = true;
    handle_groups.push_back(Handle_container());
    Handle_group handle_group = --handle_groups.end();

    insert_handle(handle_group, vd);
    return handle_group;
  }
  
  /**
   * Insert vd to provided handle_group
   */
  void insert_handle(Handle_group handle_group, vertex_descriptor vd)
  {
    need_preprocess = true;
    handle_group->push_back(vd);
  }

  template<class InputIterator>
  Handle_group insert_handle(InputIterator begin, InputIterator end)
  {
    need_preprocess = true;
    handle_groups.push_back(Handle_container());
    Handle_group handle_group = --handle_groups.end();

    insert_handle(handle_group, begin, end);
    return handle_group; 
  }

  template<class InputIterator>
  void insert_handle(Handle_group handle_group, InputIterator begin, InputIterator end)
  {
    need_preprocess = true;
    for( ;begin != end; ++begin)
    {
      insert_handle(handle_group, *begin);
    }
  }

  void erase_handle(Handle_group handle_group)
  {
    need_preprocess = true;
    handle_groups.erase(handle_group);
  }

  void erase_handle(Handle_group handle_group, vertex_descriptor vd)
  {
     need_preprocess = true;
     typename Handle_container::iterator it 
       = std::find(handle_group->begin(), handle_group->end(), vd);
     if(vd != handle_group->end())
     {
       handle_group->erase(vd);
       // Although the handle group might get empty, we do not delete it from handle_groups
     }
  }

  template<class InputIterator>
  void insert_roi(InputIterator begin, InputIterator end)
  {
    need_preprocess = true;
    for( ;begin != end; ++begin)
    {
      insert_roi(*begin);
    }
  }

  void insert_roi(vertex_descriptor vd)   
  {
    need_preprocess = true;
    roi.push_back(vd);
  }
    
////////////////////////////////////////////////////////////

  /** 
   * Necessary precomputation work before beginning deformation
   * Need to be called before translate(Handle_group handle_group, const Vector& translation), deform
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
   * Set the number of iterations used in deform()
   */
  void set_iterations(unsigned int iterations)
  {
    this->iterations = iterations;
  }

  /**
   * Set the tolerance of convergence used in deform()
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
        size_t v_index = get(vertex_index_map, *it) -1;
        solution[v_index] = original[v_index] + translation;
    }
  }

#ifdef CGAL_DEFORM_ROTATION

  template <typename Quaternion, typename Vect>
  void operator()(vertex_descriptor vd, const Point& rotation_center, const Quaternion& quat, const Vect& translation)
  {
    std::size_t idx = get(vertex_index_map, vd);
    Point p = CGAL::ORIGIN + ( original[idx] - rotation_center);
    Vect v = quat * Vect(p.x(),p.y(),p.z());
    p = Point(v[0], v[1], v[2]) + ( rotation_center - CGAL::ORIGIN); 
    p = p + Vector(translation[0],translation[1],translation[2]);
   
    solution[idx] = p;
  }
#endif // CGAL_DEFORM_ROTATION

  /**
   * Deformation on roi vertices
   * Use member variables for determining iterations and tolerance
   */
  void deform()
  {
    deform(iterations, tolerance);
  }
  /**
   * Deformation on roi vertices
   */
  void deform(unsigned int iterations, double tolerance)
  {
    CGAL_precondition(!need_preprocess); // preprocess should be called first

    double energy_this = 0;
    double energy_last;
    // iterations
    CGAL_TRACE_STREAM << "iteration started...\n";
    for ( unsigned int ite = 0; ite < iterations; ite ++)
    {
      update_solution();

#ifdef CGAL_DEFORM_EXPERIMENTAL
      optimal_rotations_polar();    // polar decomposition for optimal rotations, faster than SVD but unstable 
#else
      optimal_rotations_svd();
#endif
      /* For now close energy based termination */

      // energy_last = energy_this;
      // energy_this = energy();
      //CGAL_TRACE_STREAM << ite << " iterations: energy = " << energy_this << "\n";
      //if ( abs((energy_last-energy_this)/energy_this) < tolerance )
      //{
      //  break;
      //}
    }

    CGAL_TRACE_STREAM << "iteration end!\n";

    // copy solution to target mesh
    assign_solution();
  }

  /**
   * Reset position of deformed vertices to their original positions (i.e. positions at the time of last Deform_mesh::preprocess call)
   */
  void undo()
  {
    for(size_t i = 0; i < ros.size(); ++i)
    {
      ros[i]->point() = original[i];
    }
  }

//////////////////////// private functions ////////////////////////
private:
  void compute_edge_weight()
  {
  #ifdef CGAL_DEFORM_SPOKES_AND_RIMS
    compute_edge_weight_spokes_and_rims();
  #else
    compute_edge_weight_arap();
  #endif
  }
  // compute cotangent weights of all edges 
  void compute_edge_weight_arap()
  {
    // iterate over ros vertices and calculate weights for edges which are incident to ros
    size_t next_edge_id = 1;
    for (std::size_t i = 0; i < ros.size(); i++)
    {
      vertex_descriptor vi = ros[i];
      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
      {
        std::size_t e_idx = get(edge_index_map, *e);
        if( e_idx != 0) { continue; } // we have assigned an id already, which means we also calculted the weight
        
        put(edge_index_map, *e, next_edge_id++);
        double weight = weight_calculator(*e);
        edge_weight.push_back(weight);
      }// end of edge loop
    }// end of ros loop
  }
  void compute_edge_weight_spokes_and_rims()
  {
    // iterate over ros vertices and calculate weights for edges which are incident to ros
    size_t next_edge_id = 1;
    for (std::size_t i = 0; i < ros.size(); i++)
    {       
      vertex_descriptor vi = ros[i];
      bool is_current_rim = false;
      out_edge_iterator e, e_end;
      boost::tie(e,e_end) = boost::out_edges(ros[i], polyhedron);
      edge_descriptor active_edge = *e;

      while ( e != e_end )
      {
        std::size_t e_idx = get(edge_index_map, active_edge);
        if( e_idx == 0)  // we haven't assigned an id yet, which means we haven't calculted the weight
        {
          put(edge_index_map, active_edge, next_edge_id++);
          double weight = weight_calculator(active_edge);
          edge_weight.push_back(weight);
        }

        // loop through one spoke then one rim edge
        if(!is_current_rim && !boost::get(CGAL::edge_is_border, polyhedron, *e)) // it is rim edge's turn
        {
          is_current_rim = true;
          active_edge = CGAL::next_edge(*e, polyhedron);
        }
        else // if current edge is rim OR there is no rim edge (current spoke edge is boudary)
        {    // then iterate to next spoke edge
          is_current_rim = false;
          ++e;
          active_edge = *e;
        }
      }// end of edge loop
    }// end of ros loop
  }

  // assigns id to one rign neighbor of vd, and also push them into push_vector
  void assign_id_to_one_ring(vertex_descriptor vd, std::size_t& next_id, std::vector<vertex_descriptor>& push_vector)
  {
    in_edge_iterator e, e_end;
    for (boost::tie(e,e_end) = boost::in_edges(vd, polyhedron); e != e_end; e++)
    {
      vertex_descriptor vt = boost::source(*e, polyhedron);
      std::size_t vt_index = get(vertex_index_map, vt);
      if( vt_index == 0 )  // neighboring vertices outside roi && not visited
      {
        push_vector.push_back(vt);
        put(vertex_index_map, vt, next_id++);
      }
    }
  }
  // find region of solution, including roi and hard constraints, which is the 1-ring vertices out roi
  void region_of_solution()
  {
    outside_ros.clear();
    ros.clear();
    ros.insert(ros.end(), roi.begin(), roi.end());

    // ID assign //////////////////////////////////
    // assign ids to ROI, offset is 1
    for(std::size_t i = 0; i < roi.size(); i++)
    {
      put(vertex_index_map, roi[i], i+1);
    }
    // now assign an id (in vertex_index_map) to vertices on boundary of roi
    std::size_t next_ros_index = roi.size() + 1;
    for(std::size_t i = 0; i < roi.size(); i++)
    {
      assign_id_to_one_ring(roi[i], next_ros_index, ros);
    }
    // boundary of ros also must have ids because in SVD calculation,
    // one-ring neighbor of ROS vertices are reached. 
    for(std::size_t i = roi.size(); i < ros.size(); i++)
    {
      assign_id_to_one_ring(ros[i], next_ros_index, outside_ros);
    }
    //////////////////////////////////////////////

    // initialize the rotation matrices with the same size of ROS
    rot_mtr.resize(ros.size());
    for(std::size_t i = 0; i < rot_mtr.size(); i++)
    {
      rot_mtr[i].setIdentity();
    }
    
    solution.resize(ros.size() + outside_ros.size());
    original.resize(ros.size() + outside_ros.size());
    // initialize solution
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
      is_roi[v_index-1] = true;
    }
    for(typename Handle_group_container::iterator it_group = handle_groups.begin(); 
      it_group != handle_groups.end(); ++it_group)
    {
      for(typename Handle_container::iterator it_vertex = it_group->begin(); 
        it_vertex != it_group->end(); ++it_vertex)
      {
        size_t v_index = get(vertex_index_map, *it_vertex);
        is_hdl[v_index-1] = true;
      }
    }
  }
  void assemble_laplacian(typename SparseLinearAlgebraTraits_d::Matrix& A)
  {
  #ifdef CGAL_DEFORM_SPOKES_AND_RIMS
    assemble_laplacian_spokes_and_rims(A);
  #else
    assemble_laplacian_arap(A);
  #endif
  }
  void assemble_laplacian_spokes_and_rims(typename SparseLinearAlgebraTraits_d::Matrix& A)
  {
     /// assign cotangent Laplacian to ros vertices
    for(std::size_t i = 0; i < ros.size(); i++)
      {
        vertex_descriptor vi = ros[i];
        std::size_t vertex_idx_i = get(vertex_index_map, vi) -1;
        if ( is_roi[vertex_idx_i] && !is_hdl[vertex_idx_i] )          // vertices of ( roi - hdl )
          {
            double diagonal = 0;
            in_edge_iterator e, e_end;
            for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
              {
                vertex_descriptor vj = boost::source(*e, polyhedron);
                double wij = edge_weight[ get(edge_index_map, *e) -1];  // edge weights
                double wji = edge_weight[get(edge_index_map, CGAL::opposite_edge(*e, polyhedron))-1];

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

                std::size_t vj_index = get(vertex_index_map, vj) - 1;
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
  // Assemble Laplacian matrix A of linear system A*X=B
  void assemble_laplacian_arap(typename SparseLinearAlgebraTraits_d::Matrix& A)
  {
    /// assign cotangent Laplacian to ros vertices
    for(std::size_t i = 0; i < ros.size(); i++)
      {
        vertex_descriptor vi = ros[i];
        std::size_t vertex_idx_i = get(vertex_index_map, vi) -1;
        if ( is_roi[vertex_idx_i] && !is_hdl[vertex_idx_i] )          // vertices of ( roi - hdl )
          {
            double diagonal = 0;
            in_edge_iterator e, e_end;
            for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
              {
                vertex_descriptor vj = boost::source(*e, polyhedron);
                double wij = edge_weight[ get(edge_index_map, *e) -1];  // edge weights
                double wji = edge_weight[get(edge_index_map, CGAL::opposite_edge(*e, polyhedron))-1];
                double total_weight = wij + wji;

                std::size_t vj_index = get(vertex_index_map, vj) - 1;
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

  void optimal_rotations_svd()
  {
  #ifdef CGAL_DEFORM_SPOKES_AND_RIMS
    optimal_rotations_svd_spokes_and_rims();
  #else
    optimal_rotations_svd_arap();
  #endif
  }
  // Local step of iterations, computing optimal rotation matrices using SVD decomposition, stable
  void optimal_rotations_svd_arap()
  {
    Eigen::Matrix3d u, v;           // orthogonal matrices 
    Eigen::Vector3d w;              // singular values
    Eigen::Matrix3d cov;            // covariance matrix
    Eigen::JacobiSVD<Eigen::Matrix3d> svd;       // SVD solver         
    Eigen::Matrix3d r;
    int num_neg = 0;

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

        Vector pij = original[i] - original[vj_index -1];
        Vector qij = solution[i] - solution[vj_index -1];

        double wij = edge_weight[get(edge_index_map, *e) -1];
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
        num_neg++; 
        w = svd.singularValues();
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

    CGAL_TRACE_STREAM << num_neg << " negative rotations\n";
  }
  void optimal_rotations_svd_spokes_and_rims()
  {
    Eigen::Matrix3d u, v;           // orthogonal matrices 
    Eigen::Vector3d w;              // singular values
    Eigen::Matrix3d cov;            // covariance matrix
    Eigen::JacobiSVD<Eigen::Matrix3d> svd;       // SVD solver         
    Eigen::Matrix3d r;
    int num_neg = 0;

    // only accumulate ros vertices
    for ( std::size_t i = 0; i < ros.size(); i++ )
    {
      vertex_descriptor vi = ros[i];
      // compute covariance matrix
      cov.setZero();

      // spoke + rim edges
      bool is_current_rim = false;
      out_edge_iterator e, e_end;
      boost::tie(e,e_end) = boost::out_edges(vi, polyhedron);
      edge_descriptor active_edge = *e;

      while ( e != e_end )
      {
        vertex_descriptor v1 = boost::source(active_edge, polyhedron);
        vertex_descriptor v2 = boost::target(active_edge, polyhedron);
        
        size_t v1_index = get(vertex_index_map, v1);
        size_t v2_index = get(vertex_index_map, v2);
        
        Vector pij = original[v1_index-1] - original[v2_index-1];
        Vector qij = solution[v1_index-1] - solution[v2_index-1];

        double wij = edge_weight[get(edge_index_map, active_edge)-1];
        for (int j = 0; j < 3; j++)
        {
          for (int k = 0; k < 3; k++)
          {
            cov(j, k) += wij*pij[j]*qij[k]; 
          }
        }

        // loop through one spoke then one rim edge
        if(!is_current_rim && !boost::get(CGAL::edge_is_border, polyhedron, *e)) // it is rim edge's turn
        {
          is_current_rim = true;
          active_edge = CGAL::next_edge(*e, polyhedron);
        }
        else // if current edge is rim OR there is no rim edge (current spoke edge is boudary)
        {    // then iterate to next spoke edge
          is_current_rim = false;
          ++e;
          active_edge = *e;
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
        num_neg++; 
        w = svd.singularValues();
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

    CGAL_TRACE_STREAM << num_neg << " negative rotations\n";

  }

  void update_solution()
  {
  #ifdef CGAL_DEFORM_SPOKES_AND_RIMS
    update_solution_spokes_and_rims();
  #else
    update_solution_arap();
  #endif
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
      std::size_t vertex_idx_i = get(vertex_index_map, vi)-1;
      if ( !is_roi[vertex_idx_i] || is_hdl[vertex_idx_i] )   // hard constraints or handle vertices
      {
        Bx[i] = solution[vertex_idx_i].x(); By[i] = solution[vertex_idx_i].y(); Bz[i] = solution[vertex_idx_i].z();
      }
      else  // ( roi - handle ) vertices
      {
        Bx[i] = 0; By[i] = 0; Bz[i] = 0;
        in_edge_iterator e, e_end;
        Point& pi = original[get(vertex_index_map, vi)-1];
        for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
        {
          vertex_descriptor vj = boost::source(*e, polyhedron);
          std::size_t vj_index = get(vertex_index_map, vj) -1; 
          Vector pij =  pi - original[get(vertex_index_map, vj) -1];

          double wij = edge_weight[get(edge_index_map, *e) -1];
          double wji = edge_weight[get(edge_index_map, CGAL::opposite_edge(*e, polyhedron))-1];

          bool is_e_border = boost::get(CGAL::edge_is_border, polyhedron, *e);
          bool is_opp_of_e_border = boost::get(CGAL::edge_is_border, polyhedron, CGAL::opposite_edge(*e, polyhedron));
          double x, y, z;
          x = y = z = 0.0;
          if (is_e_border || is_opp_of_e_border)
          {  
            edge_descriptor next_of_not_border_edge = is_e_border ? CGAL::next_edge(CGAL::opposite_edge(*e, polyhedron), polyhedron)
              : CGAL::next_edge(*e, polyhedron);
            std::size_t vk_index = get(vertex_index_map, boost::target(next_of_not_border_edge, polyhedron))-1;

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
            std::size_t vk_index = get(vertex_index_map, next_e_target) -1;

            edge_descriptor next_opp_e = CGAL::next_edge(CGAL::opposite_edge(*e, polyhedron), polyhedron);
            vertex_descriptor next_opp_e_target = boost::target(next_opp_e, polyhedron);
            std::size_t vm_index = get(vertex_index_map, next_opp_e_target) -1;

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
      solution[get(vertex_index_map, ros[i])-1] = p;
    }
  }
  // Global step of iterations, updating solution
  void update_solution_arap()
  {
    typename SparseLinearAlgebraTraits_d::Vector X(ros.size()), Bx(ros.size());
    typename SparseLinearAlgebraTraits_d::Vector Y(ros.size()), By(ros.size());
    typename SparseLinearAlgebraTraits_d::Vector Z(ros.size()), Bz(ros.size());

    // assemble right columns of linear system
    for ( std::size_t i = 0; i < ros.size(); i++ )
    {
      vertex_descriptor vi = ros[i];
      std::size_t vertex_idx_i = get(vertex_index_map, vi)-1;
      if ( !is_roi[vertex_idx_i] || is_hdl[vertex_idx_i] )   // hard constraints or handle vertices
      {
        Bx[i] = solution[vertex_idx_i].x(); By[i] = solution[vertex_idx_i].y(); Bz[i] = solution[vertex_idx_i].z();
      }
      else  // ( roi - handle ) vertices
      {
        Bx[i] = 0; By[i] = 0; Bz[i] = 0;
        in_edge_iterator e, e_end;
        Point& pi = original[get(vertex_index_map, vi)-1];
        for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
        {
          vertex_descriptor vj = boost::source(*e, polyhedron);
          std::size_t vj_index = get(vertex_index_map, vj) -1; 
          Vector pij =  pi - original[get(vertex_index_map, vj) -1];
          double wij = edge_weight[get(edge_index_map, *e) -1];
          double wji = edge_weight[get(edge_index_map, CGAL::opposite_edge(*e, polyhedron))-1];

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
      solution[get(vertex_index_map, ros[i])-1] = p;
    }

  }
  // Assign solution to target mesh
  void assign_solution()
  {
    for(std::size_t i = 0; i < roi.size(); ++i){
      roi[i]->point() = solution[get(vertex_index_map, roi[i])-1];
    }
  }
  // Compute modeling energy
  double energy()
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
        Point vj_original, vj_solution;
        size_t vj_index = get(vertex_index_map, vj);
        if(vj_index == 0) // outside of ROS, just take current position (since it never changes)
        {
          vj_original = vj->point();
          vj_solution = vj->point();
        }
        else
        {
           vj_original = original[vj_index -1];
           vj_solution = solution[vj_index -1];
        }
        Vector pij = original[i] - vj_original;
        double wij = edge_weight[get(edge_index_map, *e) -1];
        Vector rot_p(0, 0, 0);                 // vector rot_i*p_ij
        for (int j = 0; j < 3; j++)
        {
          double x = rot_mtr[i](0, j) * pij[j];
          double y = rot_mtr[i](1, j) * pij[j];
          double z = rot_mtr[i](2, j) * pij[j];
          Vector v(x, y, z);
          rot_p = rot_p + v;
        }
        Vector qij = solution[i] - vj_solution;
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

