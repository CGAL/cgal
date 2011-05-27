#include "Deform_mesh.h"


Deform_mesh::Deform_mesh(Polyhedron &P)
{
	polyhedron = P;
}

Deform_mesh::~Deform_mesh(void)
{
}

void Deform_mesh::region_of_interest(Vertex_iterator begin, Vertex_iterator end, size_t k)
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
				vector<Vertex_handle> ::iterator result = find(roi.begin(), roi.end(), wh);
				if (result == roi.end())
				{
					roi.push_back(wh);
				}
				++wc;
			}while(wc != done);
		}
	}
}

void Deform_mesh::handles(Vertex_iterator begin, Vertex_iterator end)
{
	hdl.clear();
	for (Vertex_iterator vit = begin; vit != end; vit ++)
	{
		Vertex_handle handle = vit;
		hdl.push_back(handle);
	}
	hdl.push_back(end);
}

void Deform_mesh::preprocess()
{
	CGAL_TRACE("Calls preprocess()\n");

	double time_init = clock();

	double duration_assembly = 0.0;
	double duration_prefactor = 0.0;

	// get #variables
	unsigned int nb_variables = polyhedron.size_of_vertices();

	CGAL_TRACE("  Creates matrix...\n");
	// Assemble linear system A*X=B
	Taucs_solver_traits<double>::Matrix A(nb_variables); // matrix is symmetric definite positive
	Taucs_solver_traits<double>::Vector X(nb_variables), B(nb_variables);

	assemble_laplacian(A, "uni");

	duration_assembly = (clock() - time_init)/CLOCKS_PER_SEC;
	CGAL_TRACE("  Creates matrix: done (%.2lf s)\n", duration_assembly);
	CGAL_TRACE("  Pre-factorizing linear system...\n");

	// Pre-factorizing the linear system A*X=B
	time_init = clock();
	double D;
	if(!solver.pre_factor(A, D))
		return;
	duration_prefactor = (clock() - time_init)/CLOCKS_PER_SEC;

	CGAL_TRACE("  Pre-factorizing linear system: done (%.2lf s)\n", duration_prefactor);

}

void Deform_mesh::assemble_laplacian(Taucs_solver_traits<double>::Matrix& A, string type)
{
	map<Vertex_handle, int> idx;
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

void Deform_mesh::operator ()(Point p, Vector v)
{
	p = p - v;
}