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


namespace CGAL {


template <class Polyhedron, class SparseLinearAlgebraTraits_d, 
          class PolyhedronVertexDeformationIndexMap, class PolyhedronEdgeDeformationIndexMap>
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

  // Data members.
public:

  Polyhedron* polyhedron;                     // source mesh, can not be modified
  std::vector<vertex_descriptor> roi;
  std::vector<vertex_descriptor> hdl;         // user specified handles, storing the target positions
  std::vector<vertex_descriptor> ros;         // region of solution, including roi and hard constraints outside roi
  PolyhedronVertexDeformationIndexMap vertex_id_pmap;              // storing indices of all vertices 
  PolyhedronEdgeDeformationIndexMap edge_id_pmap;                  // storing indices of all edges
  std::vector<int> ros_id;                    // index of ros vertices
  std::vector<int> is_roi;                    // flag indicating vertex inside roi or not 
  std::vector<int> is_hdl;               
 
  int iterations;                                     // number of maximal iterations
  double tolerance;                                   // tolerance of convergence 

  std::vector<Eigen::Matrix3d>  rot_mtr;                       // rotation matrices of ros vertices
  std::vector<double> edge_weight;                    // weight of edges
  SparseLinearAlgebraTraits_d m_solver;               // linear sparse solver
  std::vector<Point>  solution;                       // storing position of all vertices during iterations

  // Public methods
public:

  // The constructor gets the polyhedron that we will model
  Deform_mesh(Polyhedron* P)
    :polyhedron(P), vertex_id_pmap(*P), edge_id_pmap(*P)
  {

    // initialize index maps
    vertex_iterator vb, ve;
    int idx = 0;
    for(boost::tie(vb,ve) = boost::vertices(*polyhedron); vb != ve; ++vb )
    {
      boost::put(vertex_id_pmap, *vb, idx++);
    }

    edge_iterator eb, ee;
    idx = 0;
    for(boost::tie(eb,ee) = boost::edges(*polyhedron); eb != ee; ++eb )
    {
      boost::put(edge_id_pmap, *eb, idx++);
    }

    // initialize solution
    for(boost::tie(vb,ve) = boost::vertices(*polyhedron); vb != ve; ++vb )
    {
      solution.push_back( (*vb)->point() );
    }

    // initialize flag vectors of roi, handle, ros 
    ros_id.resize(boost::num_vertices(*polyhedron), 0);
    is_roi.resize(boost::num_vertices(*polyhedron), 0);
    is_hdl.resize(boost::num_vertices(*polyhedron), 0);


    // AF: What is a good number of iterations
    // YX: I will add another new threshold according to energy value.
    iterations = 5;
    tolerance = 1e-4;       

  }

  // Release resources
  ~Deform_mesh(void)
  {
  }

  // determine the roi vertices inside k-ring for all the vertices from begin to end
  void region_of_interest(vertex_iterator begin, vertex_iterator end, size_t k)
  {
    roi.clear();
    roi.insert(roi.end(), begin, end);

    is_roi.clear();  
    is_roi.resize( boost::num_vertices(*polyhedron), 0 );    // mark all the vertices as ROI or not
    for (int i = 0; i < roi.size(); i++)
    {
      is_roi[ boost::get(vertex_id_pmap, roi[i]) ] = 1;
    }

    // AF: Why do you insert end?	 
    // Iterator ranges are usually half open

    // YX: These are not iterators of vertex list, they are just pointers to all the vertices.
    //     So we can not forget end. 
    //     Actually This API seems not very reasonable to me. Maybe we can delete it directly,
    //     since we already have the API that allows us select k-ring neighboring vertices.
    //     How do you think?
  
    int idx_lv = 0;    // pointing the neighboring vertices on current level
    int idx_lv_end;
 
    for (size_t lv = 0; lv < k; lv++)
    {
      idx_lv_end = roi.size(); 
      for ( ;idx_lv < idx_lv_end; idx_lv++ )
      {
        vertex_descriptor vd = roi[idx_lv];
        in_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::in_edges(vd, *polyhedron); e != e_end; e++)
        {
          vertex_descriptor vt = boost::source(*e, *polyhedron);
          int idx = boost::get(vertex_id_pmap, vt);
          if ( !is_roi[idx] )       // not visited yet
          {
            roi.push_back(vt);
            is_roi[idx] = 1;
          }
        }
      }
    }
    

  }


  void roi_clear()
  {
    roi.clear();
    is_roi.clear();
    is_roi.resize( boost::num_vertices(*polyhedron), 0 );
  }

  void roi_push(vertex_descriptor vd)   
  {
    int idx = boost::get(vertex_id_pmap, vd);
    if (!is_roi[idx])
    {
      roi.push_back(vd);
      is_roi[idx] = 1;
    }
  }

  // Re-assign handles from begin to end
  void handles(vertex_iterator begin, vertex_iterator end)
  {
    hdl.clear();
    hdl.insert(hdl.end(), begin, end);

    is_hdl.clear();  
    is_hdl.resize( boost::num_vertices(*polyhedron), 0 );    // mark all the vertices as handle or not
    for (int i = 0; i < hdl.size(); i++)
    {
      is_hdl[ boost::get(vertex_id_pmap, hdl[i]) ] = 1;
    }
  }

  void handle_clear()
  {
    hdl.clear();
    is_hdl.clear(); 
    is_hdl.resize( boost::num_vertices(*polyhedron), 0 );
  }

  void handle_push(vertex_descriptor vd)  
  {
    int idx = boost::get(vertex_id_pmap, vd);
    if (!is_hdl[idx])
    {
      hdl.push_back(vd);
      is_hdl[idx] = 1;
    }
  }

  // compute cotangent weights of all edges 
  void compute_edge_weight()
  {
    // refer the information that whether the weight of an edge has already computed
    std::vector<int> edge_weight_computed(boost::num_edges(*polyhedron));  

    edge_weight.clear();
    edge_weight.resize( boost::num_edges(*polyhedron), 0 ); 

    edge_iterator eb, ee;
    for(boost::tie(eb,ee) = boost::edges(*polyhedron); eb != ee; ++eb )
    {
      int e_idx = boost::get(edge_id_pmap, *eb);
      if ( !edge_weight_computed[e_idx] )
      {
        double weight = cot_weight(*eb);
        // replace cotangent weight by mean-value coordinate
        if ( weight < 0 )
        {
          weight = mean_value(*eb);
          edge_weight[e_idx] = weight;
          edge_weight_computed[e_idx] = 1;
          // assign the weights to opposite edges
          edge_descriptor e_oppo = CGAL::opposite_edge(*eb, *polyhedron);   
          int e_oppo_idx = boost::get(edge_id_pmap, e_oppo);
          edge_weight[e_oppo_idx] = mean_value(e_oppo);
          edge_weight_computed[e_oppo_idx] = 1;
        }
        else
        {
          edge_weight[e_idx] = weight;
          edge_weight_computed[e_idx] = 1;
          // assign the weights to opposite edges
          edge_descriptor e_oppo = CGAL::opposite_edge(*eb, *polyhedron);   
          int e_oppo_idx = boost::get(edge_id_pmap, e_oppo);
          edge_weight[e_oppo_idx] = weight;
          edge_weight_computed[e_oppo_idx] = 1;
        }
        
      }


    }
  }

  // find region of solution, including roi and hard constraints, which is the 1-ring vertices out roi
  void region_of_solution()
  {
    ros.clear();
    ros.insert(ros.end(),roi.begin(), roi.end());

    // initialize the indices of ros vertices 
    ros_id.clear();
    ros_id.resize(boost::num_vertices(*polyhedron), -1);
    int ros_idx;
    for ( ros_idx = 0; ros_idx < ros.size(); ros_idx++)
    {
      ros_id[ boost::get(vertex_id_pmap, ros[ros_idx]) ] = ros_idx; 
    }

    for (int i = 0;i < roi.size(); i++)
    {
      vertex_descriptor vd = roi[i];
      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(vd, *polyhedron); e != e_end; e++)
      {
        vertex_descriptor vt = boost::source(*e, *polyhedron);
        int idx = boost::get(vertex_id_pmap, vt);
        if ( !is_roi[idx] && ros_id[idx] == -1 )    // neighboring vertices outside roi && not visited 
        {
          ros.push_back(vt);
          ros_id[idx] = ros_idx++;
        }
      }
    }

    // initialize the rotation matrices with the same size of ROS
    rot_mtr.clear();
    rot_mtr.resize(ros.size());
    for (int i = 0; i < ros.size(); i++)
    {
      rot_mtr[i].setIdentity();
    }
    
  }

  // Before we can model we have to do some precomputation
  ///
  /// @commentheading Template parameters:
  /// @param SparseLinearAlgebraTraits_d Symmetric definite positive sparse linear solver.
  void preprocess()
  {
    CGAL_TRACE_STREAM << "Calls preprocess()\n";

    Timer task_timer; task_timer.start();


    CGAL_TRACE_STREAM << "  Creates matrix...\n";

    compute_edge_weight();
    region_of_solution();

    // Assemble linear system A*X=B
    typename SparseLinearAlgebraTraits_d::Matrix A(ros.size()); // matrix is symmetric definite positive
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
    // initialize the Laplacian matrix
    for (int i = 0; i < ros.size(); i++)
      {
        A.set_coef(i, i, 1.0, true);     
      }

    /// assign cotangent Laplacian to ros vertices
    for(int i = 0; i < ros.size(); i++)
      {
        vertex_descriptor vi = ros[i];
        int vertex_idx_i = boost::get(vertex_id_pmap, vi);
        if ( is_roi[vertex_idx_i] && !is_hdl[vertex_idx_i] )          // vertices of ( roi - hdl )
          {
            double diagonal = 0;
            in_edge_iterator e, e_end;
            for (boost::tie(e,e_end) = boost::in_edges(vi, *polyhedron); e != e_end; e++)
              {
                vertex_descriptor vj = boost::source(*e, *polyhedron);
                double wij = edge_weight[ boost::get(edge_id_pmap, *e) ];  // cotangent Laplacian weights
                int ros_idx_j = ros_id[ boost::get(vertex_id_pmap, vj) ];
                A.set_coef(i, ros_idx_j, -wij, true);	// off-diagonal coefficient
                diagonal += wij;  
              }
            // diagonal coefficient
            A.set_coef(i, i, diagonal);
          } 
      }

  }

  
  // Returns the cotangent weight of specified edge_descriptor
  double cot_weight(edge_descriptor e)
  {
     vertex_descriptor v0 = boost::target(e, *polyhedron);
     vertex_descriptor v1 = boost::source(e, *polyhedron);
     // Only one triangle for border edges
     if (boost::get(CGAL::edge_is_border, *polyhedron, e) ||
         boost::get(CGAL::edge_is_border, *polyhedron, CGAL::opposite_edge(e, *polyhedron)))
     {
       
       edge_descriptor e_cw = CGAL::next_edge_cw(e, *polyhedron);
       vertex_descriptor v2 = boost::source(e_cw, *polyhedron);
       if (boost::get(CGAL::edge_is_border, *polyhedron, e_cw) ||
           boost::get(CGAL::edge_is_border, *polyhedron, CGAL::opposite_edge(e_cw, *polyhedron)) )
       {
          edge_descriptor e_ccw = CGAL::next_edge_ccw(e, *polyhedron);
          v2 = boost::source(e_ccw, *polyhedron);
       }
      
       return ( cot_value(v0, v2, v1)/2.0 );
     }
     else
     {
        edge_descriptor e_cw = CGAL::next_edge_cw(e, *polyhedron);
        vertex_descriptor v2 = boost::source(e_cw, *polyhedron);     
        edge_descriptor e_ccw = CGAL::next_edge_ccw(e, *polyhedron);
        vertex_descriptor v3 = boost::source(e_ccw, *polyhedron);

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

    return (cos_angle/sin_angle);

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
    vertex_descriptor v0 = boost::target(e, *polyhedron);
    vertex_descriptor v1 = boost::source(e, *polyhedron);
    Vector vec = v0->point() - v1->point();
    double norm = std::sqrt( vec.squared_length() );

    // Only one triangle for border edges
    if (boost::get(CGAL::edge_is_border, *polyhedron, e) ||
        boost::get(CGAL::edge_is_border, *polyhedron, CGAL::opposite_edge(e, *polyhedron)))
    {

      edge_descriptor e_cw = CGAL::next_edge_cw(e, *polyhedron);
      vertex_descriptor v2 = boost::source(e_cw, *polyhedron);
      if (boost::get(CGAL::edge_is_border, *polyhedron, e_cw) || 
          boost::get(CGAL::edge_is_border, *polyhedron, CGAL::opposite_edge(e_cw, *polyhedron)) )
      {
        edge_descriptor e_ccw = CGAL::next_edge_ccw(e, *polyhedron);
        v2 = boost::source(e_ccw, *polyhedron);
      }

      return ( half_tan_value(v1, v0, v2)/norm );
    }
    else
    {
      edge_descriptor e_cw = CGAL::next_edge_cw(e, *polyhedron);
      vertex_descriptor v2 = boost::source(e_cw, *polyhedron);     
      edge_descriptor e_ccw = CGAL::next_edge_ccw(e, *polyhedron);
      vertex_descriptor v3 = boost::source(e_ccw, *polyhedron);

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
    for (int idx = 0; idx < hdl.size(); idx++)
      {
        vertex_descriptor vd = hdl[idx];
        solution[boost::get(vertex_id_pmap, vd)] = vd->point() + translation;
      }
  }

  // The operator will be called in a real time loop from the GUI.
  // assign translation vector to specific handle
  void operator()(vertex_descriptor vd, const Vector& translation)
  {
    int idx = boost::get(vertex_id_pmap, vd);
    solution[idx] = vd->point() + translation;
  }


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
    for ( int i = 0; i < ros.size(); i++ )
    {
      vertex_descriptor vi = ros[i];
      // compute covariance matrix
      cov.setZero();

      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(vi, *polyhedron); e != e_end; e++)
      {
        vertex_descriptor vj = boost::source(*e, *polyhedron);
        Vector pij = vi->point() - vj->point();
        Vector qij = solution[boost::get(vertex_id_pmap, vi)] - solution[boost::get(vertex_id_pmap, vj)];
        double wij = edge_weight[boost::get(edge_id_pmap, *e)];
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
            u(0, j0) = -1.0*u(0, j0);
            u(1, j0) = -1.0*u(1, j0);
            u(2, j0) = -1.0*u(2, j0);
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
    int k = -1;
    Eigen::Matrix3d Y;
    double alpha, beta, gamma;
    do 
    {
      k++;
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
    for ( int i = 0; i < ros.size(); i++ )
    {
      vertex_descriptor vi = ros[i];
      // compute covariance matrix
      cov.setZero();

      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(vi, *polyhedron); e != e_end; e++)
      {
        vertex_descriptor vj = boost::source(*e, *polyhedron);
        Vector pij = vi->point() - vj->point();
        Vector qij = solution[boost::get(vertex_id_pmap, vi)] - solution[boost::get(vertex_id_pmap, vj)];
        double wij = edge_weight[boost::get(edge_id_pmap, *e)];
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
            u(0, j0) = -1.0*u(0, j0);
            u(1, j0) = -1.0*u(1, j0);
            u(2, j0) = -1.0*u(2, j0);
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
    for ( int i = 0; i < ros.size(); i++ )
    {
      vertex_descriptor vi = ros[i];
      int vertex_idx_i = boost::get(vertex_id_pmap, vi);
      if ( !is_roi[vertex_idx_i] || is_hdl[vertex_idx_i] )   // hard constraints or handle vertices
      {
        Bx[i] = solution[vertex_idx_i].x(); By[i] = solution[vertex_idx_i].y(); Bz[i] = solution[vertex_idx_i].z();
      }
      else  // ( roi - handle ) vertices
      {
        Bx[i] = 0; By[i] = 0; Bz[i] = 0;
        in_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::in_edges(vi, *polyhedron); e != e_end; e++)
        {
          vertex_descriptor vj = boost::source(*e, *polyhedron);
          int ros_idx_j = ros_id[ boost::get(vertex_id_pmap, vj) ]; 
          Vector pij = vi->point() - vj->point();
          double wij = edge_weight[boost::get(edge_id_pmap, *e)];
          Vector rot_p(0, 0, 0);                  // vector ( r_i + r_j )*p_ij
          for (int j = 0; j < 3; j++)
          {
            double x = ( rot_mtr[i](0, j) + rot_mtr[ros_idx_j](0, j) ) * pij[j];
            double y = ( rot_mtr[i](1, j) + rot_mtr[ros_idx_j](1, j) ) * pij[j];
            double z = ( rot_mtr[i](2, j) + rot_mtr[ros_idx_j](2, j) ) * pij[j];
            Vector v(x, y, z);
            rot_p = rot_p + v;
          }
          Vector vec = wij*rot_p/2.0;
          Bx[i] += vec.x(); By[i] += vec.y(); Bz[i] += vec.z(); 
        }
      }
    }

    // solve "A*X = B".
    m_solver.solve(Bx, X); m_solver.solve(By, Y); m_solver.solve(Bz, Z);

    // copy to solution
    for (int i = 0; i < ros.size(); i++)
    {
      Point p(X[i], Y[i], Z[i]);
      solution[boost::get(vertex_id_pmap, ros[i])] = p;
    }

  }

  // Compute modeling energy
  double energy()
  {
    double sum_of_energy = 0;
    // only accumulate ros vertices
    for( int i = 0; i < ros.size(); i++ )
    {
      vertex_descriptor vi = ros[i];
      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(vi, *polyhedron); e != e_end; e++)
      {
        vertex_descriptor vj = boost::source(*e, *polyhedron);
        Vector pij = vi->point() - vj->point();
        double wij = edge_weight[boost::get(edge_id_pmap, *e)];
        Vector rot_p(0, 0, 0);                 // vector rot_i*p_ij
        for (int j = 0; j < 3; j++)
        {
          double x = rot_mtr[i](0, j) * pij[j];
          double y = rot_mtr[i](1, j) * pij[j];
          double z = rot_mtr[i](2, j) * pij[j];
          Vector v(x, y, z);
          rot_p = rot_p + v;
        }
        Vector qij = solution[boost::get(vertex_id_pmap, vi)] - solution[boost::get(vertex_id_pmap, vj)];
        sum_of_energy += wij*(qij - rot_p).squared_length();
      }
    }
    return sum_of_energy;
  }

  // Deformation on roi vertices
  void deform(Polyhedron* P)
  {
    double energy_this = 0;
    double energy_last;
    // iterations
    CGAL_TRACE_STREAM << "iteration started...\n";
    for ( int ite = 0; ite < iterations; ite ++)
    {
      update_solution();
#ifdef CGAL_DEFORM_EXPERIMENTAL
      optimal_rotations_polar();    // polar decomposition for optimal rotations, faster than SVD but unstable 
#else
      optimal_rotations_svd();
#endif
      energy_last = energy_this;
      energy_this = energy();
      CGAL_TRACE_STREAM << ite << " iterations: energy = " << energy_this << "\n";
      if ( abs((energy_last-energy_this)/energy_this) < tolerance )
      {
        break;
      }
    }
    CGAL_TRACE_STREAM << "iteration end!\n";

    // copy solution to target mesh P
    assign_solution(P);
  }

  // Assign solution to target mesh P
  void assign_solution(Polyhedron* P)
  {
    vertex_iterator vb, ve;
    boost::tie(vb,ve) = boost::vertices(*P);
    for ( int i = 0; i < boost::num_vertices(*P); i++ )
    {
      if (is_roi[i])     // only copy ROI vertices
      {
        (*vb)->point() = solution[i];
      }
      vb++;
    }
  }

  // Undo: reset P to be the copied polyhedron
  void undo(Polyhedron* P)
  {
    vertex_iterator vb, ve;
    boost::tie(vb,ve) = boost::vertices(*P);
    vertex_iterator vb_copy, ve_copy;
    boost::tie(vb_copy,ve_copy) = boost::vertices(*polyhedron);
    for ( int i = 0; i < boost::num_vertices(*P); i++ )
    {
      if (is_roi[i])     // only copy ROI vertices
      {
        (*vb)->point() = (*vb_copy)->point();
        solution[i] = (*vb_copy)->point();
      }
      vb++;
      vb_copy++;
    }
  }
};


} //namespace CGAL

#endif  // CGAL_DEFORM_MESH_H

