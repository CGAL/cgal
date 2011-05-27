#pragma once

#include <CGAL/trace.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Taucs_solver_traits.h>
#include <iostream>
#include <string>
#include <fstream>
#include <ctime>

typedef CGAL::Simple_cartesian<double>      Kernel;
typedef Kernel::Vector_3                    Vector;
typedef Kernel::Point_3                     Point;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;

typedef Polyhedron::Vertex_const_handle                      Vertex_handle;
typedef Polyhedron::Vertex_iterator							 Vertex_iterator;
typedef Polyhedron::Halfedge_around_vertex_const_circulator  HV_circulator;

using namespace std;
using namespace CGAL;

class Deform_mesh
{
private:
	Polyhedron polyhedron;                // target mesh
	vector<Vertex_handle> roi;
	vector<Vertex_handle> hdl;
	vector<Vertex_handle> dsplc;         // displacement of handles

	Taucs_solver_traits<double> solver;

public:
	// The constructor gets the Polyhedron that we will model
	Deform_mesh(Polyhedron &P);

	// Release ressources
	~Deform_mesh(void);

	// The region of interest and the handles are a set of vertices
	void region_of_interest(Vertex_iterator begin, Vertex_iterator end, size_t k);

	void handles(Vertex_iterator begin, Vertex_iterator end);

	// Before we can model we have to do some precomputation
	void preprocess();

	// Assemble Laplacian matrix A of linear system A*X=B 
	void assemble_laplacian(Taucs_solver_traits<double>::Matrix& A, string type);

	// The operator will be called in a real time loop from the GUI.
	void operator()(Point p, Vector v);

};
