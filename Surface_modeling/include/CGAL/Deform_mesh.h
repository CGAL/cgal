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

#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <CGAL/FPU_extension.h>

#include <Eigen/Eigen>
#include <Eigen/SVD>

#include <limits>

namespace CGAL {

template <class Polyhedron, class SparseLinearAlgebraTraits_d, 
          class VertexIndexMap, class EdgeIndexMap>
class Deform_mesh
{
// Public types
public:

  // Geometric types              
  typedef typename Polyhedron::Traits         Kernel;
  typedef typename Kernel::Vector_3           Vector;
  typedef typename Kernel::Point_3            Point;

  // Repeat Polyhedron types
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor		vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_iterator		  vertex_iterator;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor		  edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::edge_iterator		    edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::in_edge_iterator		in_edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::out_edge_iterator		out_edge_iterator;
  // Data members.
public:

  Polyhedron& polyhedron;															// source mesh, can not be modified

  VertexIndexMap vertex_index_map;										// storing indices of ros vertices, others should be 0
  EdgeIndexMap   edge_index_map;										  // storing indices of ros related edges, others should be 0

  std::vector<vertex_descriptor> hdl;									// user specified handles
  std::vector<vertex_descriptor> roi;                 // region of interest, including both free and haldle vertices
  std::vector<vertex_descriptor> ros;									// region of solution, including roi and hard constraints outside roi
  std::vector<vertex_descriptor> outside_ros;         // boundary of ros, for clearing purpose

  // properties per ros vertex
  std::vector<bool> is_roi;
	std::vector<bool> is_hdl; 
  std::vector<Eigen::Matrix3d> rot_mtr;               // rotation matrices of ros vertices
  std::vector<Point> solution;                        // storing position of ros vertices during iterations
  std::vector<Point> original;

  std::vector<double> edge_weight;                    // weight of edges only those who are incident to ros 

  SparseLinearAlgebraTraits_d m_solver;               // linear sparse solver
  unsigned int iterations;                            // number of maximal iterations
  double tolerance;                                   // tolerance of convergence 

  // Public methods
public:

  // The constructor gets the polyhedron that we will model
  Deform_mesh(Polyhedron& P, const VertexIndexMap& vertex_index_map_, const EdgeIndexMap& edge_index_map_)
    :polyhedron(P), vertex_index_map(vertex_index_map_), edge_index_map(edge_index_map_)
  {

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
      boost::put(vertex_index_map, *vb, 0);
    }
    // this part should be removed same as above
    edge_iterator eb, ee;
    for(boost::tie(eb,ee) = boost::edges(polyhedron); eb != ee; ++eb )
    {
      boost::put(edge_index_map, *eb, 0);
    }
    /////////////////////////////////////////////////////////////////

    // AF: What is a good number of iterations
    // YX: I will add another new threshold according to energy value.
    iterations = 5;
    tolerance = 1e-4;       
  }


  void clear()
  {
    //clear vertices
    roi.clear(); 
    hdl.clear();
    for (std::size_t i = 0; i < ros.size(); i++)
    {
      boost::put(vertex_index_map, ros[i], 0);
    }
    for (std::size_t i = 0; i < outside_ros.size(); i++)
    {
      boost::put(vertex_index_map, outside_ros[i], 0);
    }
    // note that cleaning/reassigning is_roi, is_hdl, ros, solution, original vectors
    // handled in region_of_solution().

    //clear edges

    edge_weight.clear();
    for (std::size_t i = 0; i < ros.size(); i++)
    {
      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(ros[i], polyhedron); e != e_end; e++)
      {
        boost::put(edge_index_map, *e, 0);
        edge_descriptor e_oppo = CGAL::opposite_edge(*e, polyhedron);
        boost::put(edge_index_map, e_oppo, 0);
      }
    }
  }

  void insert_roi(vertex_descriptor vd)   
  {
    roi.push_back(vd);
  }

  // Re-assign handles from begin to end
  void assign_handles(vertex_iterator begin, vertex_iterator end)
  {
    hdl.clear();
    hdl.insert(hdl.end(), begin, end);
  }

  void insert_handle(vertex_descriptor vd)  
  {
    hdl.push_back(vd);
  }

  // compute cotangent weights of all edges 
 

  void compute_edge_weight()
  {
    // iterate over ros vertices and calculate weights for edges which are incident to ros
    size_t next_edge_id = 1;
    for (std::size_t i = 0; i < ros.size(); i++)
    {
      vertex_descriptor vi = ros[i];
      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
      {
        std::size_t e_idx = boost::get(edge_index_map, *e);
        if( e_idx != 0) { continue; } // we have assigned an id already, which means we also calculted the weight
        
        boost::put(edge_index_map, *e, next_edge_id++);
        double weight = cot_weight(*e);
        // replace cotangent weight by mean-value coordinate
        if ( weight < 0 )
        {
          weight = mean_value(*e);
          edge_weight.push_back(weight);
          // assign the weights to opposite edges
          edge_descriptor e_oppo = CGAL::opposite_edge(*e, polyhedron);   
          std::size_t e_oppo_idx = boost::get(edge_index_map, e_oppo);
          if(e_oppo_idx == 0)
          {
            boost::put(edge_index_map, e_oppo, next_edge_id++);
            edge_weight.push_back(mean_value(e_oppo));            
          }
        }
        else
        {
          edge_weight.push_back(weight);
          // assign the weights to opposite edges
          edge_descriptor e_oppo = CGAL::opposite_edge(*e, polyhedron);   
          std::size_t e_oppo_idx = boost::get(edge_index_map, e_oppo);
          if(e_oppo_idx == 0)
          {
            boost::put(edge_index_map, e_oppo, next_edge_id++);
            edge_weight.push_back(weight);            
          }
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
      std::size_t vt_index = boost::get(vertex_index_map, vt);
      if( vt_index == 0 )  // neighboring vertices outside roi && not visited
      {
        push_vector.push_back(vt);
        boost::put(vertex_index_map, vt, next_id++);
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
    // assign ids to ROI - offset is 1
    for(std::size_t i = 0; i < roi.size(); i++)
    {
      boost::put(vertex_index_map, roi[i], i+1);
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
      size_t v_index = boost::get(vertex_index_map, roi[i]);
      is_roi[v_index-1] = true;
    }
    for(std::size_t i = 0; i < hdl.size(); i++)
    {
      size_t v_index = boost::get(vertex_index_map, hdl[i]);
      is_hdl[v_index-1] = true;
    }
  }

  // Before we can model we have to do some precomputation
  ///
  /// @commentheading Template parameters:
  /// @param SparseLinearAlgebraTraits_d Definite positive sparse linear solver.
  void preprocess()
  {
    CGAL_TRACE_STREAM << "Calls preprocess()\n";

    Timer task_timer; task_timer.start();


    CGAL_TRACE_STREAM << "  Creates matrix...\n";

    
    region_of_solution();
    compute_edge_weight(); // compute_edge_weight() has to come later then region_of_solution()

    // Assemble linear system A*X=B
    typename SparseLinearAlgebraTraits_d::Matrix A(ros.size()); // matrix is definite positive, and not necessarily symmetric
    assemble_laplacian(A);		

    CGAL_TRACE_STREAM << "  Creates " << ros.size() << "*" << ros.size() << " matrix: done (" << task_timer.time() << " s)\n";

    CGAL_TRACE_STREAM << "  Pre-factorizing linear system...\n";

    // Pre-factorizing the linear system A*X=B
    task_timer.reset();
    double D;
    if(!m_solver.pre_factor(A, D))
      return;

    CGAL_TRACE_STREAM << "  Pre-factorizing linear system: done (" << task_timer.time() << " s)\n";

  }


  // Assemble Laplacian matrix A of linear system A*X=B
  ///
  /// @commentheading Template parameters:
  /// @param SparseLinearAlgebraTraits_d definite positive sparse linear solver.
  void assemble_laplacian(typename SparseLinearAlgebraTraits_d::Matrix& A)
  {
    /// assign cotangent Laplacian to ros vertices
    for(std::size_t i = 0; i < ros.size(); i++)
      {
        vertex_descriptor vi = ros[i];
        std::size_t vertex_idx_i = boost::get(vertex_index_map, vi) -1;
        if ( is_roi[vertex_idx_i] && !is_hdl[vertex_idx_i] )          // vertices of ( roi - hdl )
          {
            double diagonal = 0;
            in_edge_iterator e, e_end;
            for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
              {
                vertex_descriptor vj = boost::source(*e, polyhedron);
                double wij = edge_weight[ boost::get(edge_index_map, *e) -1];  // cotangent Laplacian weights
                double wji = edge_weight[boost::get(edge_index_map, CGAL::opposite_edge(*e, polyhedron))-1];
                double total_weight = wij + wji;
                std::size_t vj_index = boost::get(vertex_index_map, vj) - 1;
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

  
  // Returns the cotangent weight of specified edge_descriptor
  double cot_weight(edge_descriptor e)
  {
     vertex_descriptor v0 = boost::target(e, polyhedron);
     vertex_descriptor v1 = boost::source(e, polyhedron);
     // Only one triangle for border edges
     if (boost::get(CGAL::edge_is_border, polyhedron, e) ||
         boost::get(CGAL::edge_is_border, polyhedron, CGAL::opposite_edge(e, polyhedron)))
     {
       
       edge_descriptor e_cw = CGAL::next_edge_cw(e, polyhedron);
       vertex_descriptor v2 = boost::source(e_cw, polyhedron);
       if (boost::get(CGAL::edge_is_border, polyhedron, e_cw) ||
           boost::get(CGAL::edge_is_border, polyhedron, CGAL::opposite_edge(e_cw, polyhedron)) )
       {
          edge_descriptor e_ccw = CGAL::next_edge_ccw(e, polyhedron);
          v2 = boost::source(e_ccw, polyhedron);
       }
      
       return ( cot_value(v0, v2, v1)/2.0 );
     }
     else
     {
        edge_descriptor e_cw = CGAL::next_edge_cw(e, polyhedron);
        vertex_descriptor v2 = boost::source(e_cw, polyhedron);     
        edge_descriptor e_ccw = CGAL::next_edge_ccw(e, polyhedron);
        vertex_descriptor v3 = boost::source(e_ccw, polyhedron);

        return ( cot_value(v0, v2, v1)/2.0 + cot_value(v0, v3, v1)/2.0 );
     }
  }

  // Returns the cotangent value of angle v0_v1_v2
  double cot_value(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    
    Vector vec0 = v1->point() - v2->point();
    Vector vec1 = v2->point() - v0->point();
    Vector vec2 = v0->point() - v1->point();
    double e0_square = vec0.squared_length();
    double e1_square = vec1.squared_length();
    double e2_square = vec2.squared_length();
    double e0 = std::sqrt(e0_square); 
    double e2 = std::sqrt(e2_square);
    double cos_angle = ( e0_square + e2_square - e1_square ) / 2.0 / e0 / e2;
    double sin_angle = std::sqrt(1-cos_angle*cos_angle);

    return (cos_angle/sin_angle) / std::sqrt(squared_area(v0->point(), v1->point(), v2->point()));

  }

  // Returns the tangent value of half angle v0_v1_v2/2
  double half_tan_value(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {

    Vector vec0 = v1->point() - v2->point();
    Vector vec1 = v2->point() - v0->point();
    Vector vec2 = v0->point() - v1->point();
    double e0_square = vec0.squared_length();
    double e1_square = vec1.squared_length();
    double e2_square = vec2.squared_length();
    double e0 = std::sqrt(e0_square); 
    double e2 = std::sqrt(e2_square);
    double cos_angle = ( e0_square + e2_square - e1_square ) / 2.0 / e0 / e2;
    double angle = acos(cos_angle);

    return ( tan(angle/2.0) );

  }

  // Returns the mean-value coordinate of specified edge_descriptor
  double mean_value(edge_descriptor e)
  {
    vertex_descriptor v0 = boost::target(e, polyhedron);
    vertex_descriptor v1 = boost::source(e, polyhedron);
    Vector vec = v0->point() - v1->point();
    double norm = std::sqrt( vec.squared_length() );

    // Only one triangle for border edges
    if (boost::get(CGAL::edge_is_border, polyhedron, e) ||
        boost::get(CGAL::edge_is_border, polyhedron, CGAL::opposite_edge(e, polyhedron)))
    {

      edge_descriptor e_cw = CGAL::next_edge_cw(e, polyhedron);
      vertex_descriptor v2 = boost::source(e_cw, polyhedron);
      if (boost::get(CGAL::edge_is_border, polyhedron, e_cw) || 
          boost::get(CGAL::edge_is_border, polyhedron, CGAL::opposite_edge(e_cw, polyhedron)) )
      {
        edge_descriptor e_ccw = CGAL::next_edge_ccw(e, polyhedron);
        v2 = boost::source(e_ccw, polyhedron);
      }

      return ( half_tan_value(v1, v0, v2)/norm );
    }
    else
    {
      edge_descriptor e_cw = CGAL::next_edge_cw(e, polyhedron);
      vertex_descriptor v2 = boost::source(e_cw, polyhedron);     
      edge_descriptor e_ccw = CGAL::next_edge_ccw(e, polyhedron);
      vertex_descriptor v3 = boost::source(e_ccw, polyhedron);

      return ( half_tan_value(v1, v0, v2)/norm + half_tan_value(v1, v0, v3)/norm );
    }
  }

  // Set the number of iterations made in operator()
  void set_iterations(unsigned int ite)
  {
    iterations = ite;
  }

  // Set the tolerance of convergence made in operator()
  void set_tolerance(double tole)
  {
    tolerance = tole;
  }
  
  // The operator will be called in a real time loop from the GUI.
  // assign translation vector to all handles
  void operator()(const Vector& translation)
  {
    for (std::size_t idx = 0; idx < hdl.size(); idx++)
      {
        vertex_descriptor vd = hdl[idx];
        solution[boost::get(vertex_index_map, vd)-1] = 
            original[boost::get(vertex_index_map,vd)-1] + translation;
      }
  }

  // The operator will be called in a real time loop from the GUI.
  // assign translation vector to specific handle
  void operator()(vertex_descriptor vd, const Vector& translation)
  {
    std::size_t idx = boost::get(vertex_index_map, vd);
    solution[idx-1] = original[idx-1] + translation;
  }

#ifdef CGAL_DEFORM_ROTATION

  template <typename Quaternion, typename Vect>
  void operator()(vertex_descriptor vd, const Point& rotation_center, const Quaternion& quat, const Vect& translation)
  {
    std::size_t idx = boost::get(vertex_index_map, vd);
    Point p = CGAL::ORIGIN + ( original[idx] - rotation_center);
    Vect v = quat * Vect(p.x(),p.y(),p.z());
    p = Point(v[0], v[1], v[2]) + ( rotation_center - CGAL::ORIGIN); 
    p = p + Vector(translation[0],translation[1],translation[2]);
   
    solution[idx] = p;
  }
#endif // CGAL_DEFORM_ROTATION
  // Local step of iterations, computing optimal rotation matrices using SVD decomposition, stable
  void optimal_rotations_svd()
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
        size_t vj_index = boost::get(vertex_index_map, vj);

        Vector pij = original[i] - original[vj_index -1];
        Vector qij = solution[i] - solution[vj_index -1];

        double wij = edge_weight[boost::get(edge_index_map, *e) -1];
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
        Vector pij = original[boost::get(vertex_index_map, vi)] - original[boost::get(vertex_index_map, vj)];
        Vector qij = solution[boost::get(vertex_index_map, vi)] - solution[boost::get(vertex_index_map, vj)];
        double wij = edge_weight[boost::get(edge_index_map, *e)];
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


  // Global step of iterations, updating solution
  void update_solution()
  {
    typename SparseLinearAlgebraTraits_d::Vector X(ros.size()), Bx(ros.size());
    typename SparseLinearAlgebraTraits_d::Vector Y(ros.size()), By(ros.size());
    typename SparseLinearAlgebraTraits_d::Vector Z(ros.size()), Bz(ros.size());

    // assemble right columns of linear system
    for ( std::size_t i = 0; i < ros.size(); i++ )
    {
      vertex_descriptor vi = ros[i];
      std::size_t vertex_idx_i = boost::get(vertex_index_map, vi)-1;
      if ( !is_roi[vertex_idx_i] || is_hdl[vertex_idx_i] )   // hard constraints or handle vertices
      {
        Bx[i] = solution[vertex_idx_i].x(); By[i] = solution[vertex_idx_i].y(); Bz[i] = solution[vertex_idx_i].z();
      }
      else  // ( roi - handle ) vertices
      {
        Bx[i] = 0; By[i] = 0; Bz[i] = 0;
        in_edge_iterator e, e_end;
        Point& pi = original[boost::get(vertex_index_map, vi)-1];
        for (boost::tie(e,e_end) = boost::in_edges(vi, polyhedron); e != e_end; e++)
        {
          vertex_descriptor vj = boost::source(*e, polyhedron);
          std::size_t vj_index = boost::get(vertex_index_map, vj) -1; 
          Vector pij =  pi - original[boost::get(vertex_index_map, vj) -1];
          double wij = edge_weight[boost::get(edge_index_map, *e) -1];
          double wji = edge_weight[boost::get(edge_index_map, CGAL::opposite_edge(*e, polyhedron))-1];
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
      solution[boost::get(vertex_index_map, ros[i])-1] = p;
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
        size_t vj_index = boost::get(vertex_index_map, vj);
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
        double wij = edge_weight[boost::get(edge_index_map, *e) -1];
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

  // Deformation on roi vertices
  void deform()
  {
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
      // for now close energy based termination.
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

  // Assign solution to target mesh
  void assign_solution()
  {
    for(std::size_t i = 0; i < roi.size(); ++i){
      roi[i]->point() = solution[boost::get(vertex_index_map, roi[i])-1];
    }
  }

  // Undo: reset P to be the copied polyhedron
  void undo()
  {
    for(size_t i = 0; i < ros.size(); ++i)
    {
      ros[i]->point() = original[i];
    }
  }
};


} //namespace CGAL

#endif  // CGAL_DEFORM_MESH_H

