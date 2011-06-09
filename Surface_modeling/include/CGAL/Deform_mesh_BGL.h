#ifndef CGAL_DEFORM_MESH_BGL_H
#define CGAL_DEFORM_MESH_BGL_H


#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/Taucs_solver_traits.h>


#include <iostream>
#include <list>
#include <fstream>






namespace CGAL {

/// @heading Parameters:
/// @param Gt Geometric traits class.

template <class Polyhedron>
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
  typedef typename boost::graph_traits<Polyhedron>::out_edge_iterator		out_edge_iterator;


  // Data members.
public:

  Polyhedron polyhedron;                // target mesh
  std::map<vertex_descriptor, vertex_descriptor> s2t;          // access from source mesh to target mesh
  std::vector<vertex_descriptor> roi;
  std::vector<vertex_descriptor> hdl;           // user specified handles, storing the target positions


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
  template <class SparseLinearAlgebraTraits_d>
  void preprocess(
    SparseLinearAlgebraTraits_d& solver = SparseLinearAlgebraTraits_d())
  {
    CGAL_TRACE_STREAM << "Calls preprocess()\n";

    Timer task_timer; task_timer.start();

    // get #variables
    unsigned int nb_variables = boost::num_vertices(polyhedron);

    CGAL_TRACE_STREAM << "  Creates matrix...\n";
    // Assemble linear system A*X=B
    typename SparseLinearAlgebraTraits_d::Matrix A(nb_variables); // matrix is symmetric definite positive
    typename SparseLinearAlgebraTraits_d::Vector X(nb_variables), B(nb_variables);

    assemble_laplacian<SparseLinearAlgebraTraits_d>(A, "uni");

    CGAL_TRACE_STREAM << "  Creates matrix: done (" << task_timer.time() << " s)\n";

    CGAL_TRACE_STREAM << "  Pre-factorizing linear system...\n";

    // Pre-factorizing the linear system A*X=B
    task_timer.reset();
    double D;
    if(!solver.pre_factor(A, D))
      return;

    CGAL_TRACE_STREAM << "  Pre-factorizing linear system: done (" << task_timer.time() << " s)\n";

  }


  void preprocess()
  {
    return preprocess< Taucs_solver_traits<double> >();
  }

  void preprocess( Taucs_solver_traits<double>& solver )
  {
    return preprocess< Taucs_solver_traits<double> >(solver);
  }

	// Assemble Laplacian matrix A of linear system A*X=B 
  ///
  /// @commentheading Template parameters:
  /// @param SparseLinearAlgebraTraits_d Symmetric definite positive sparse linear solver.
  template <class SparseLinearAlgebraTraits_d>
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
        double wij = 1;
        if (type == "cot")   // cotangent Laplacian weights
        {
          ;
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

