#ifndef CGAL_DEFORM_MESH_BGL_H
#define CGAL_DEFORM_MESH_BGL_H


#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <CGAL/svd.h>



#include <iostream>
#include <list>
#include <fstream>


namespace CGAL {

/// @heading Parameters:
/// @param Gt Geometric traits class.

template <class Polyhedron, class SparseLinearAlgebraTraits_d>
class Deform_mesh_BGL
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
  

  int iterations;                                     // number of iterations
  std::map<vertex_descriptor, std::vector<double> > rot_mtr;               // rotation matrices of all vertices
  std::map<edge_descriptor, double> edge_weight;           // weight of edges
  SparseLinearAlgebraTraits_d m_solver;               // linear sparse solver
  std::map<vertex_descriptor, Point> solution;        // storing vertex solution during iterations
  std::map<vertex_descriptor, int> ros_idx;           // index of ros vertices


  // Public methods
public:

  // The constructor gets the polyhedron that we will model
  Deform_mesh_BGL(Polyhedron* P)
    :polyhedron(P)
  {

    vertex_iterator vb, ve;
    for(boost::tie(vb,ve) = boost::vertices(*polyhedron); vb != ve; ++vb )
    {
      rot_mtr[*vb].resize(9);
    }

    for(boost::tie(vb,ve) = boost::vertices(*polyhedron); vb != ve; ++vb )
    {
      solution[*vb] = (*vb)->point();
    }

    iterations = 10;

  }

  // Release resources
  ~Deform_mesh_BGL(void)
  {
  }

  // determine the roi vertices inside k-ring for all the vertices from begin to end
  void region_of_interest(vertex_iterator begin, vertex_iterator end, size_t k)
  {
    roi.clear();
    for (vertex_iterator vit = begin; vit != end; vit ++)
    {
      roi.push_back(*vit);
    }
    roi.push_back(*end);
  
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
          std::vector<vertex_descriptor> ::iterator result = find(roi.begin(), roi.end(), vt);
          if (result == roi.end())
          {
            roi.push_back(vt);
          }
        }
      }
    }
    

  }


  void roi_clear()
  {
    roi.clear();
  }

  void roi_push(vertex_descriptor vd)   
  {
    roi.push_back(vd);
  }


  void handles(vertex_iterator begin, vertex_iterator end)
  {
    hdl.clear();
    for (vertex_iterator vit = begin; vit != end; vit ++)
    {
      hdl.push_back(*vit);
    }
    hdl.push_back(*end);
  }

  void handle_clear()
  {
    hdl.clear();
  }

  void handle_push(vertex_descriptor vd)  
  {
    hdl.push_back(vd);
  }


  // find region of solution, including roi and hard constraints, which is the 1-ring vertices out roi
  void region_of_solution()
  {
    ros.clear();
    for (int i = 0; i < roi.size(); i++)
    {
      ros.push_back(roi[i]);
    }

    std::vector<vertex_descriptor> hard_constraints;
    for (int i = 0;i < roi.size(); i++)
    {
      vertex_descriptor vd = roi[i];
      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(vd, *polyhedron); e != e_end; e++)
      {
        vertex_descriptor vt = boost::source(*e, *polyhedron);
        std::vector<vertex_descriptor> ::iterator result_roi = find(roi.begin(), roi.end(), vt);
        std::vector<vertex_descriptor> ::iterator result = find(hard_constraints.begin(), hard_constraints.end(), vt);
        if ( result == hard_constraints.end() && result_roi == roi.end() )    // neighboring vertices outside roi 
        {
          hard_constraints.push_back(vt);
        }
      }
    }

    for (int i = 0; i < hard_constraints.size(); i++)
    {
      ros.push_back(hard_constraints[i]);
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
    assemble_laplacian(A, "cot");

    CGAL_TRACE_STREAM << "  Creates matrix: done (" << task_timer.time() << " s)\n";

    CGAL_TRACE_STREAM << "  Pre-factorizing linear system...\n";

    // Pre-factorizing the linear system A*X=B
    task_timer.reset();
    double D;
    if(!m_solver.pre_factor(A, D))
      return;

    CGAL_TRACE_STREAM << "  Pre-factorizing linear system: done (" << task_timer.time() << " s)\n";

  }

  // pre-compute cotangent weight of edges
  void compute_edge_weight()
  {
    edge_iterator eb, ee;
    for(boost::tie(eb,ee) = boost::edges(*polyhedron); eb != ee; ++eb )
    {
      edge_weight[*eb] = cot_value(*eb) / 2.0;
    }
  }

	// Assemble Laplacian matrix A of linear system A*X=B
  ///
  /// @commentheading Template parameters:
  /// @param SparseLinearAlgebraTraits_d definite positive sparse linear solver.
  void assemble_laplacian(typename SparseLinearAlgebraTraits_d::Matrix& A, std::string type)
	{
    // initialize the indices of roi vertices and Laplacian matrix
    ros_idx.clear();
    for (int i = 0; i < ros.size(); i++)
    {
      ros_idx[ros[i]] = i;
      A.set_coef(i, i, 1.0, true);     
    }

    /// assign cotangent Laplacian to ros vertices
		for(int i = 0; i < ros.size(); i++)
		{
			vertex_descriptor vi = ros[i];
      std::vector<vertex_descriptor> ::iterator result_roi = find(roi.begin(), roi.end(), vi);
      std::vector<vertex_descriptor> ::iterator result_hdl = find(hdl.begin(), hdl.end(), vi);
      if ( result_roi != roi.end() && result_hdl == hdl.end() )  // vertices of ( roi - hdl )
      {
        double diagonal = 0;
        int idx_i = ros_idx[vi];
        in_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::in_edges(vi, *polyhedron); e != e_end; e++)
        {
          vertex_descriptor vj = boost::source(*e, *polyhedron);
          // only account neighboring vertices inside ros
          std::vector<vertex_descriptor> ::iterator result_neigh = find(ros.begin(), ros.end(), vj);
          if ( result_neigh != ros.end() )
          {
            double wij = 1;
            if (type == "cot")                    // cotangent Laplacian weights
            {
              wij = edge_weight[*e];
            }
            int idx_j = ros_idx[vj];
            A.set_coef(idx_i, idx_j, -wij, true);	// off-diagonal coefficient
            diagonal += wij;
          }    
        }
        // diagonal coefficient
        A.set_coef(idx_i, idx_i, diagonal);
      } 
		}

	}

  
  // Returns the cotanget value of specified edge_descriptor
  double cot_value(edge_descriptor e)
  {
     vertex_descriptor v0 = boost::target(e, *polyhedron);
     vertex_descriptor v1 = boost::source(e, *polyhedron);
     // Only one triangle for border edges
     if (boost::get(CGAL::edge_is_border, *polyhedron, e)||boost::get(CGAL::edge_is_border, *polyhedron, CGAL::opposite_edge(e, *polyhedron)))
     {
       
       edge_descriptor e_cw = CGAL::next_edge_cw(e, *polyhedron);
       vertex_descriptor v2 = boost::source(e_cw, *polyhedron);
       if (boost::get(CGAL::edge_is_border, *polyhedron, e_cw) || boost::get(CGAL::edge_is_border, *polyhedron, CGAL::opposite_edge(e_cw, *polyhedron)) )
       {
          edge_descriptor e_ccw = CGAL::next_edge_ccw(e, *polyhedron);
          v2 = boost::source(e_ccw, *polyhedron);
       }
      
       return ( cot_value(v0, v2, v1) );
     }
     else
     {
        edge_descriptor e_cw = CGAL::next_edge_cw(e, *polyhedron);
        vertex_descriptor v2 = boost::source(e_cw, *polyhedron);
        edge_descriptor e_ccw = CGAL::next_edge_ccw(e, *polyhedron);
        vertex_descriptor v3 = boost::source(e_ccw, *polyhedron);

        return ( cot_value(v0, v2, v1) + cot_value(v0, v3, v1) );
     }
  }

  // Returns the cotanget value of angle v0_v1_v2
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


  // Set the number of iterations made in operator()
  void set_iterations(unsigned int ite)
  {
    iterations = ite;
  }
  
	// The operator will be called in a real time loop from the GUI.
  // assign translation vector to handles
	void operator()(Vector translation)
	{
    for (int idx = 0; idx < hdl.size(); idx++)
    {
      vertex_descriptor vd = hdl[idx];
      solution[vd] = vd->point() + translation;
    }
	}

  // Local step of iterations, computing optimal rotation matrices
  void optimal_rotations()
  {
    // malloc memory for svd decomposition
    float** u = new float*[4];    // covariance matrix
    for (int i = 0; i < 4; i++)
    {
      u[i] = new float[4];
    }
    float** v = new float*[4];
    for (int i = 0; i < 4; i++)
    {
      v[i] = new float[4];
    }
    float* w = new float[4];

    // only accumulate ros vertices
    int idx = 0;
    for ( int i = 0; i < ros.size(); i++ )
    {
      vertex_descriptor vi = ros[i];
      // compute covariance matrix
      for (int i = 0; i < 4; i++)
      {
        w[i] = 0;
        for (int j = 0; j < 4; j++)
        {
          u[i][j] = 0;
          v[i][j] = 0;
        }
      }
      in_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::in_edges(vi, *polyhedron); e != e_end; e++)
      {
        vertex_descriptor vj = boost::source(*e, *polyhedron);
        Vector pij = vi->point() - vj->point();
        Vector qij = solution[vi] - solution[vj];
        double wij = edge_weight[*e];
        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            u[i+1][j+1] += wij*pij[i]*qij[j]; 
          }
        }
      }

      // svd decomposition
      svdcmp(u, 3, 3, w, v);

      // extract rotation matrix
      for (int i = 0; i < 9; i++)
      {
        rot_mtr[vi][i] = 0;
      }
      for (int i = 0; i < 3; i++)      // row index
      {
        for (int j = 0; j < 3; j++)   // column index
        {
          for (int k = 0; k < 3; k++)
          {
            rot_mtr[vi][i*3+j] += v[i+1][k+1] * u[j+1][k+1];
          }  
        }
      }
    }
  
    for (int i = 0; i < 4; i++)
    {
      delete [] u[i];
      delete [] v[i];
    }
    delete [] u;
    delete [] v;
    delete [] w;

  }

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
      std::vector<vertex_descriptor> ::iterator result_roi = find(roi.begin(), roi.end(), vi);
      std::vector<vertex_descriptor> ::iterator result_hdl = find(hdl.begin(), hdl.end(), vi);
      if ( result_roi == roi.end() || result_hdl != hdl.end() )   // handle vertices or hard constraints
      {
        Bx[ros_idx[vi]] = solution[vi].x(); By[ros_idx[vi]] = solution[vi].y(); Bz[ros_idx[vi]] = solution[vi].z();
      }
      else  // roi vertices
      {
        Bx[ros_idx[vi]] = 0; By[ros_idx[vi]] = 0;Bz[ros_idx[vi]] = 0;
        in_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::in_edges(vi, *polyhedron); e != e_end; e++)
        {
          vertex_descriptor vj = boost::source(*e, *polyhedron);
          // only account neighboring vertices inside ros
          std::vector<vertex_descriptor> ::iterator result_neigh = find(ros.begin(), ros.end(), vj);
          if ( result_neigh != ros.end() )
          {
            Vector pij = vi->point() - vj->point();
            double wij = edge_weight[*e];
            Vector rot_p(0, 0, 0);                 // vector ( r_i + r_j )*p_ij
            for (int j = 0; j < 3; j++)
            {
              double x = ( rot_mtr[vi][j] + rot_mtr[vj][j] ) * pij[j];
              double y = ( rot_mtr[vi][3+j] + rot_mtr[vj][3+j] ) * pij[j];
              double z = ( rot_mtr[vi][6+j] + rot_mtr[vj][6+j] ) * pij[j];
              Vector v(x, y, z);
              rot_p = rot_p + v;
            }
            Vector vec = wij*rot_p/2.0;
            Bx[ros_idx[vi]] += vec.x(); By[ros_idx[vi]] += vec.y(); Bz[ros_idx[vi]] += vec.z();
          }  
        }
      }
    }

    // solve "A*X = B".
    m_solver.solve(Bx, X); m_solver.solve(By, Y); m_solver.solve(Bz, Z);

    // copy to solution
    for (int i = 0; i < ros.size(); i++)
    {
      vertex_descriptor vd = ros[i];
      Point p(X[ros_idx[vd]], Y[ros_idx[vd]], Z[ros_idx[vd]]);
      solution[vd] = p;
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
        double wij = edge_weight[*e];
        Vector rot_p(0, 0, 0);                 // vector rot_i*p_ij
        for (int j = 0; j < 3; j++)
        {
          double x = rot_mtr[vi][j] * pij[j];
          double y = rot_mtr[vi][3+j] * pij[j];
          double z = rot_mtr[vi][6+j] * pij[j];
          Vector v(x, y, z);
          rot_p = rot_p + v;
        }
        Vector qij = solution[vi] - solution[vj];
        sum_of_energy += wij*(qij - rot_p).squared_length();
      }
    }
    return sum_of_energy;
  }

  // Deformation on roi vertices
  void deform(Polyhedron* P)
  {
    // iterations
    CGAL_TRACE_STREAM << "iteration started...\n";
    optimal_rotations();
    for ( int ite = 0; ite < iterations; ite ++)
    {
      update_solution();
      optimal_rotations();
      CGAL_TRACE_STREAM << ite << " iterations: energy = " << energy() << "\n";
    }
    CGAL_TRACE_STREAM << "iteration end!\n";

    // copy solution to target mesh P
    vertex_iterator vb_t, ve_t;
    boost::tie(vb_t,ve_t) = boost::vertices(*P);

    vertex_iterator vb, ve;
    for(boost::tie(vb,ve) = boost::vertices(*polyhedron); vb != ve; ++vb )
    {
      (*vb_t)->point() = solution[*vb];
      vb_t++;
    }
  }
};


} //namespace CGAL

#endif  // CGAL_DEFORM_MESH_BGL_H

