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

typedef Polyhedron::Vertex_const_handle                      Vertex_handle;
typedef Polyhedron::Vertex_iterator							             Vertex_iterator;
typedef Polyhedron::Halfedge_around_vertex_const_circulator  HV_circulator;


namespace CGAL {

class Deform_mesh
{
private:
	Polyhedron polyhedron;                // target mesh
  std::vector<Vertex_handle> roi;
  std::vector<Vertex_handle> hdl;
  std::vector<Vertex_handle> dsplc;         // displacement of handles

  Taucs_solver_traits<double> solver;

public:
	// The constructor gets the Polyhedron that we will model
	Deform_mesh(Polyhedron &P)
	{
		polyhedron = P;
	}

	// Release ressources
	~Deform_mesh(void)
	{
	}

	// The region of interest and the handles are a set of vertices
	void region_of_interest(Vertex_iterator begin, Vertex_iterator end, size_t k)
	{
		roi.clear();
		for (Vertex_iterator vit = begin; vit != end; vit ++)
		{
			Vertex_handle handle = vit;
			roi.push_back(handle);
		}
		roi.push_back(end);

		int idx_lv = 0;    // pointing the neighboring vertices on current level
		int idx_lv_end;

		for (size_t lv = 0; lv < k; lv++)
		{
			idx_lv_end = roi.size(); 
			for ( ;idx_lv < idx_lv_end; idx_lv++ )
			{
				Vertex_handle vh = roi[idx_lv];
				HV_circulator wc = vh->vertex_begin(), done(wc);
				do {
					Vertex_handle wh = wc->opposite()->vertex();
          std::vector<Vertex_handle> ::iterator result = find(roi.begin(), roi.end(), wh);
					if (result == roi.end())
					{
						roi.push_back(wh);
					}
					++wc;
				}while(wc != done);
			}
		}
	}


	void handles(Vertex_iterator begin, Vertex_iterator end)
	{
		hdl.clear();
		for (Vertex_iterator vit = begin; vit != end; vit ++)
		{
			Vertex_handle handle = vit;
			hdl.push_back(handle);
		}
		hdl.push_back(end);
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
    std::map<Vertex_handle, int> idx;
		int index = 0;
		for (Vertex_iterator vit = polyhedron.vertices_begin(); vit != polyhedron.vertices_end(); vit++)
		{
			idx[vit] = index;
			index++;
		}

		for (Vertex_iterator vit = polyhedron.vertices_begin(); vit != polyhedron.vertices_end(); vit++)
		{
			Vertex_handle vi = vit;
			double diagonal = 0;
			int idx_i = idx[vi];
			HV_circulator wc = vi->vertex_begin(), done(wc);
			do {
				Vertex_handle vj = wc->opposite()->vertex();
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
	void operator()(Point p, Vector v)
	{
		p = p - v;
	}
};


} //namespace CGAL

#endif  // CGAL_DEFORM_MESH_H

