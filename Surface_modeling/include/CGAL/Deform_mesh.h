#ifndef CGAL_DEFORM_MESH_H
#define CGAL_DEFORM_MESH_H


#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Taucs_solver_traits.h>


typedef CGAL::Simple_cartesian<double>      Kernel;
typedef Kernel::Vector_3                    Vector;
typedef Kernel::Point_3                     Point;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;

typedef Polyhedron::Vertex_handle                            Vertex_handle;
typedef Polyhedron::Vertex_const_handle                      Vertex_const_handle;
typedef Polyhedron::Vertex_iterator							             Vertex_iterator;
typedef Polyhedron::Vertex_const_iterator                    Vertex_const_iterator;
typedef Polyhedron::Halfedge_around_vertex_const_circulator  HV_circulator;


namespace CGAL {

class Deform_mesh
{
public:
	Polyhedron polyhedron;                // target mesh
  std::map<Vertex_const_handle, Vertex_handle> s2t;          // access from source mesh to target mesh
  std::vector<Vertex_handle> roi;
  std::vector<Vertex_handle> hdl;           // user specified handles, storing the target positions

  Taucs_solver_traits<double> solver;

public:
	// The constructor gets the Polyhedron that we will model
	Deform_mesh(Polyhedron &P)
    :polyhedron(P)
	{
    Vertex_iterator vit_t = polyhedron.vertices_begin();
		for (Vertex_const_iterator vit_s = P.vertices_begin(); vit_s != P.vertices_end(); vit_s++)
		{
      s2t[vit_s] = vit_t;
		}
	}

	// Release resources
	~Deform_mesh(void)
	{
	}

	// The region of interest and the handles are a set of vertices on target mesh, iterators come from source mesh
	void region_of_interest(Vertex_const_iterator begin, Vertex_const_iterator end, size_t k)
	{
    std::vector<Vertex_const_handle> roi_s;
		for (Vertex_const_iterator vit = begin; vit != end; vit ++)
		{
      roi_s.push_back(vit);
		}
    roi_s.push_back(end);

		int idx_lv = 0;    // pointing the neighboring vertices on current level
		int idx_lv_end;

		for (size_t lv = 0; lv < k; lv++)
		{
			idx_lv_end = roi_s.size(); 
			for ( ;idx_lv < idx_lv_end; idx_lv++ )
			{
				Vertex_const_handle vh = roi_s[idx_lv];
				HV_circulator wc = vh->vertex_begin(), done(wc);
				do {
					Vertex_const_handle wh = wc->opposite()->vertex();
          std::vector<Vertex_const_handle> ::iterator result = find(roi_s.begin(), roi_s.end(), wh);
					if (result == roi_s.end())
					{
						roi_s.push_back(wh);
					}
					++wc;
				}while(wc != done);
			}
		}

    // mapping handles from source to target
    roi.clear();
    for (int i = 0; i < roi_s.size(); i++)
    {
      roi.push_back(s2t[roi_s[i]]);
    }

	}


	void handles(Vertex_const_iterator begin, Vertex_const_iterator end)
	{
		hdl.clear();
		for (Vertex_const_iterator vit = begin; vit != end; vit ++)
		{
			hdl.push_back(s2t[vit]);
		}
		hdl.push_back(s2t[end]);
	}


	// Before we can model we have to do some precomputation
	void preprocess()
	{
		CGAL_TRACE_STREAM << "Calls preprocess()\n";

		Timer task_timer; task_timer.start();

		// get #variables
		unsigned int nb_variables = polyhedron.size_of_vertices();

		CGAL_TRACE_STREAM << "  Creates matrix...\n";
		// Assemble linear system A*X=B
    Taucs_solver_traits<double>::Matrix A(nb_variables); // matrix is symmetric definite positive
    Taucs_solver_traits<double>::Vector X(nb_variables), B(nb_variables);

		assemble_laplacian(A, "uni");

		CGAL_TRACE_STREAM << "  Creates matrix: done (" << task_timer.time() << " s)\n";

		CGAL_TRACE_STREAM << "  Pre-factorizing linear system...\n";

		// Pre-factorizing the linear system A*X=B
		task_timer.reset();
		double D;
		if(!solver.pre_factor(A, D))
			return;

		CGAL_TRACE_STREAM << "  Pre-factorizing linear system: done (" << task_timer.time() << " s)\n";

	}

	// Assemble Laplacian matrix A of linear system A*X=B 
	void assemble_laplacian(Taucs_solver_traits<double>::Matrix& A, std::string type)
	{
    std::map<Vertex_const_handle, int> idx;
		int index = 0;
		for (Vertex_const_iterator vit = polyhedron.vertices_begin(); vit != polyhedron.vertices_end(); vit++)
		{
			idx[vit] = index;
			index++;
		}

		for (Vertex_const_iterator vit = polyhedron.vertices_begin(); vit != polyhedron.vertices_end(); vit++)
		{
			Vertex_const_handle vi = vit;
			double diagonal = 0;
			int idx_i = idx[vi];
			HV_circulator wc = vi->vertex_begin(), done(wc);
			do {
				Vertex_const_handle vj = wc->opposite()->vertex();
				double wij = 1;
				if (type == "cot")   // cotangent Laplacian weights
				{
					;
				}
				int idx_j = idx[vj];
				A.set_coef(idx_i, idx_j, -wij, true);	// off-diagonal coefficient
				diagonal += wij;
				++wc;
			}while(wc != done);
			// diagonal coefficient
			A.set_coef(idx_i, idx_i, diagonal, true);
		}

		// handle constraints
		for (int i = 0; i < hdl.size(); i++)
		{
			int idx_i = idx[hdl[i]];
			A.set_coef(idx_i, idx_i, 1.0);
		}
	}

	// The operator will be called in a real time loop from the GUI.
	void operator()(Vertex_const_iterator vit, Vector v)
	{
    Point p = s2t[vit]->point();
		p = p-v;
    s2t[vit]->point() = p;
	}
};


} //namespace CGAL

#endif  // CGAL_DEFORM_MESH_H

