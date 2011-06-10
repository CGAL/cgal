#ifndef CGAL_DEFORM_MESH_BGL_H
#define CGAL_DEFORM_MESH_BGL_H


#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <CGAL/Taucs_solver_traits.h>


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
  typedef CGAL::Cartesian<double>::Point_3                                              Point;
  typedef CGAL::Cartesian<double>::Vector_3                                             Vector;

  // Repeat Polyhedron types
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor		vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_iterator		  vertex_iterator;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor		  edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::edge_iterator		    edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::in_edge_iterator		in_edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::out_edge_iterator		out_edge_iterator;


  // Data members.
public:

  Polyhedron polyhedron;                // target mesh
  std::map<vertex_descriptor, vertex_descriptor> s2t;          // access from source mesh to target mesh
  std::vector<vertex_descriptor> roi;
  std::vector<vertex_descriptor> hdl;           // user specified handles, storing the target positions

  SparseLinearAlgebraTraits_d m_solver;


  // Public methods
public:

  // The constructor gets the Polyhedron that we will model
  Deform_mesh_BGL(Polyhedron &P)
    :polyhedron(P)
  {
    int index = 0;
    vertex_iterator vb_s, ve_s, vb_t, ve_t;
    // boost::tie assigns the first and second element of the std::pair
    // returned by boost::vertices to the variables vb and ve
    boost::tie(vb_s,ve_s) = boost::vertices(P);
    for(boost::tie(vb_t,ve_t) = boost::vertices(polyhedron); vb_t != ve_t; ++vb_t ){
      vertex_descriptor  vd_t = *vb_t;
      vd_t->id() = index;
      vertex_descriptor  vd_s = *vb_s;
      vd_s->id() = index;
      s2t[vd_s] = vd_t;
      vb_s++;
      index++;
    }

  }

  // Release resources
  ~Deform_mesh_BGL(void)
  {
  }

  // The region of interest and the handles are a set of vertices on target mesh, iterators come from source mesh
  void region_of_interest(vertex_iterator begin, vertex_iterator end, size_t k)
  {
    roi.clear();
    for (vertex_iterator vit = begin; vit != end; vit ++)
    {
      roi.push_back(s2t[*vit]);
    }
    roi.push_back(s2t[*end]);
  
    int idx_lv = 0;    // pointing the neighboring vertices on current level
    int idx_lv_end;
    
    for (size_t lv = 0; lv < k; lv++)
    {
      idx_lv_end = roi.size(); 
      for ( ;idx_lv < idx_lv_end; idx_lv++ )
      {
        vertex_descriptor vd = roi[idx_lv];
        out_edge_iterator e, e_end;
        for (boost::tie(e,e_end) = boost::out_edges(vd, polyhedron); e != e_end; e++)
        {
          vertex_descriptor vt = boost::target(*e, polyhedron);
          std::vector<vertex_descriptor> ::iterator result = find(roi.begin(), roi.end(), vt);
          if (result == roi.end())
          {
            roi.push_back(vt);
          }
        }
      }
    }
    

  }


  void handles(vertex_iterator begin, vertex_iterator end)
  {
    hdl.clear();
    for (vertex_iterator vit = begin; vit != end; vit ++)
    {
      hdl.push_back(s2t[*vit]);
    }
    hdl.push_back(s2t[*end]);
  }


  // Before we can model we have to do some precomputation
  ///
  /// @commentheading Template parameters:
  /// @param SparseLinearAlgebraTraits_d Symmetric definite positive sparse linear solver.
  void preprocess()
  {
    CGAL_TRACE_STREAM << "Calls preprocess()\n";

    Timer task_timer; task_timer.start();

    // get #variables
    unsigned int nb_variables = boost::num_vertices(polyhedron);

    CGAL_TRACE_STREAM << "  Creates matrix...\n";
    // Assemble linear system A*X=B
    typename SparseLinearAlgebraTraits_d::Matrix A(nb_variables); // matrix is symmetric definite positive
    typename SparseLinearAlgebraTraits_d::Vector X(nb_variables), B(nb_variables);

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



	// Assemble Laplacian matrix A of linear system A*X=B 
  ///
  /// @commentheading Template parameters:
  /// @param SparseLinearAlgebraTraits_d Symmetric definite positive sparse linear solver.
	void assemble_laplacian(typename SparseLinearAlgebraTraits_d::Matrix& A, std::string type)
	{
    vertex_iterator vb, ve;
		for(boost::tie(vb,ve) = boost::vertices(polyhedron); vb != ve; ++vb )
		{
			vertex_descriptor vi = *vb;
      double diagonal = 0;
			int idx_i = vi->id();
      out_edge_iterator e, e_end;
      for (boost::tie(e,e_end) = boost::out_edges(vi, polyhedron); e != e_end; e++)
      {
        vertex_descriptor vj = boost::target(*e, polyhedron);
        //std::cout << vj->id();
        double wij = 1;
        if (type == "cot")   // cotangent Laplacian weights
        {
          wij = cot_value(*e);
        }
        int idx_j = vj->id();
        A.set_coef(idx_i, idx_j, -wij, true);	// off-diagonal coefficient
        diagonal += wij;
      }
			// diagonal coefficient
			A.set_coef(idx_i, idx_i, diagonal, true);
		}

		// handle constraints
		for (int i = 0; i < hdl.size(); i++)
		{
			int idx_i = hdl[i]->id();
			A.set_coef(idx_i, idx_i, 1.0);
		}
	}

  
  // Returns the cotanget value of specified edge_descriptor
  double cot_value(edge_descriptor e)
  {
     vertex_descriptor v0 = boost::source(e, polyhedron);
     vertex_descriptor v1 = boost::target(e, polyhedron);
     // Only one triangle for border edges
     if (boost::get(CGAL::edge_is_border, polyhedron, e)||boost::get(CGAL::edge_is_border, polyhedron, CGAL::opposite_edge(e, polyhedron)))
     {
       
       edge_descriptor e_cw = CGAL::next_edge_cw(e, polyhedron);
       vertex_descriptor v2 = boost::target(e_cw, polyhedron);
       if (boost::get(CGAL::edge_is_border, polyhedron, e_cw) || boost::get(CGAL::edge_is_border, polyhedron, CGAL::opposite_edge(e_cw, polyhedron)) )
       {
          edge_descriptor e_ccw = CGAL::next_edge_ccw(e, polyhedron);
          v2 = boost::target(e_ccw, polyhedron);
       }
      
       //std::cout << v0->id() << " " << v1->id() << " " << v2->id();
       return cot_value(v0, v2, v1);
     }
     else
     {
        edge_descriptor e_cw = CGAL::next_edge_cw(e, polyhedron);
        vertex_descriptor v2 = boost::target(e_cw, polyhedron);
        edge_descriptor e_ccw = CGAL::next_edge_ccw(e, polyhedron);
        vertex_descriptor v3 = boost::target(e_ccw, polyhedron);

        //std::cout << v0->id() << " " << v1->id() << " " << v2->id() << " " << v3->id();

        return ( cot_value(v0, v2, v1)/2.0 + cot_value(v0, v3, v1)/2.0 );
     }
  }

  // Returns the cotanget value of angle v0_v1_v2
  double cot_value(vertex_descriptor v0, vertex_descriptor v1, vertex_descriptor v2)
  {
    Point p0 = v0->point();
    Point p1 = v1->point();
    Point p2 = v2->point();
    Vector vec0 = v1->point() - v2->point();
    Vector vec1 = v2->point() - v0->point();
    Vector vec2 = v0->point() - v1->point();
    double e0 = std::sqrt(vec0.squared_length());
    double e1 = std::sqrt(vec1.squared_length());
    double e2 = std::sqrt(vec2.squared_length());
    double cos_angle = (e0*e0+e2*e2-e1*e1)/2.0/e0/e2;
    double sin_angle = std::sqrt(1-cos_angle*cos_angle);

    return (cos_angle/sin_angle);

  }
  
	// The operator will be called in a real time loop from the GUI.
	void operator()(vertex_iterator vit, Vector v)
	{
    Point p = s2t[*vit]->point();
		p = p-v;
    s2t[*vit]->point() = p;
    std::cout << hdl[0]->point();
	}
};


} //namespace CGAL

#endif  // CGAL_DEFORM_MESH_BGL_H

