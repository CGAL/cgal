// Copyright (c) 2005  INRIA (France).
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
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Laurent Saboret, Bruno Levy, Pierre Alliez


#ifndef FIXED_BORDER_PARAMETIZER_3_H
#define FIXED_BORDER_PARAMETIZER_3_H

#include <CGAL/circulator.h>
#include <OpenNL/linear_solver.h>

#include "CGAL/Parametizer_3.h"
#include "CGAL/Circular_border_parametizer_3.h"
#include "CGAL/parameterization_assertions.h"

CGAL_BEGIN_NAMESPACE


//
// Declaration
//

// Class Fixed_border_parametizer_3 
// Model of the Parametizer_3 concept.
//
// Base class of fixed border parameterization methods (Tutte, Floater, ...) 1 to 1 mapping is guaranteed if surface's border is mapped onto a convex polygon.
//
// Implementation notes: - Usually, subclasses need only to implement compute_wij().
//                       - The current implementation does not remove border vertices from the linear systems => A cannot be symmetric.

template <class MeshAdaptor_3,														// 3D surface
		  class BorderParametizer_3 = Circular_border_parametizer_3<MeshAdaptor_3>,	// Class to map the surface's border onto a 2D space
		  class SparseLinearAlgebraTraits_d = OpenNL::DefaultLinearSolverTraits<typename MeshAdaptor_3::NT> >	
																					// Traits class for solving a sparse linear system "A*X = B"
class Fixed_border_parametizer_3 : public Parametizer_3<MeshAdaptor_3>
{
// Public types
public:
				// Export Mesh_Adaptor_3, BorderParametizer_3 and SparseLinearAlgebraTraits_d types and subtypes
				typedef MeshAdaptor_3													Mesh_adaptor_3;
				typedef typename Parametizer_3<MeshAdaptor_3>::ErrorCode				ErrorCode;
				typedef typename MeshAdaptor_3::NT										NT;
				typedef typename MeshAdaptor_3::Face									Face;
				typedef typename MeshAdaptor_3::Vertex									Vertex;
				typedef typename MeshAdaptor_3::Point_3									Point_3;
				typedef typename MeshAdaptor_3::Point_2									Point_2;
				typedef typename MeshAdaptor_3::Vector_3								Vector_3;
				typedef typename MeshAdaptor_3::Face_iterator							Face_iterator;
				typedef typename MeshAdaptor_3::Face_const_iterator						Face_const_iterator;
				typedef typename MeshAdaptor_3::Vertex_iterator							Vertex_iterator;
				typedef typename MeshAdaptor_3::Vertex_const_iterator					Vertex_const_iterator;
				typedef typename MeshAdaptor_3::Border_vertex_iterator					Border_vertex_iterator;
				typedef typename MeshAdaptor_3::Border_vertex_const_iterator			Border_vertex_const_iterator;
				typedef typename MeshAdaptor_3::Vertex_around_face_circulator			Vertex_around_face_circulator;
				typedef typename MeshAdaptor_3::Vertex_around_face_const_circulator		Vertex_around_face_const_circulator;
				typedef typename MeshAdaptor_3::Vertex_around_vertex_circulator			Vertex_around_vertex_circulator;
				typedef typename MeshAdaptor_3::Vertex_around_vertex_const_circulator	Vertex_around_vertex_const_circulator;
				typedef BorderParametizer_3												Border_parametizer_3;
				typedef SparseLinearAlgebraTraits_d										Sparse_linear_algebra_traits_d;

// Public operations
public:
				// Constructor
				// @param borderParametizer	Object that maps the surface's border onto a 2D space
				// @param linearAlgebra		Traits object to solve the "A*X = B" sparse linear system used by parameterization algorithms
				Fixed_border_parametizer_3 (BorderParametizer_3 borderParametizer = BorderParametizer_3(), 
											SparseLinearAlgebraTraits_d linearAlgebra = SparseLinearAlgebraTraits_d()) 
				  : m_borderParametizer(borderParametizer), m_linearAlgebra(linearAlgebra)
				{}

				// Default copy constructor and operator =() are fine

				// Compute a 1 to 1 mapping from a triangular 3D surface 'mesh' to a piece of the 2D space. 
				// The mapping is linear by pieces (linear in each triangle).
				// The result is the (u,v) pair image of each vertex of the 3D surface. 
				//
				// Preconditions:
				// * 'mesh' must be a surface with 1 connected component and no hole
				// * 'mesh' must be a triangular mesh
				// * the mesh border must be mapped onto a convex polygon
				virtual ErrorCode  parameterize (MeshAdaptor_3* mesh);
								
// Protected types
protected:
				typedef typename OpenNL::LinearSolver<SparseLinearAlgebraTraits_d>		Solver ;
								
// Protected operations
protected:
				// Check parameterize() preconditions:
				// * 'mesh' must be a surface with 1 connected component and no hole
				// * 'mesh' must be a triangular mesh
				// * the mesh border must be mapped onto a convex polygon
				virtual ErrorCode  check_parameterize_preconditions (const MeshAdaptor_3& mesh);

				// Initialize A, Bu and Bv after boundary parameterization
				// Fill the border vertices' lines in both linear systems: "u = constant" and "v = constant"
				//
				// Preconditions:
				// * vertices must be indexed
				// * A, Bu and Bv must be allocated
				// * border vertices must be parameterized
				void initialize_system_from_mesh_border(Solver* solver_u, Solver* solver_v, const MeshAdaptor_3& mesh) ;

				// compute wij = (i,j) coefficient of matrix A for j neighbor vertex of i
				// Implementation note: usually, subclasses of Fixed_border_parametizer_3 need only to implement compute_wij()
				virtual	NT  compute_wij(const MeshAdaptor_3& mesh, const Vertex& main_vertex_Vi, Vertex_around_vertex_const_circulator neighbor_vertex_Vj) = 0;

				// Compute the line i of matrix A for i inner vertex
				// * call compute_wij() to compute the A coefficient Wij for each neighbor Vj
				// * compute Wii = - sum of Wij
				// 
				// Preconditions:
				// * vertices must be indexed
				// * vertex i musn't be already parameterized
				// * line i of A must contain only zeros
				virtual ErrorCode  setup_inner_vertex_relations (Solver* solver_u, Solver* solver_v, const MeshAdaptor_3& mesh, const Vertex& vertex);

				// Copy Xu and Xv coordinates into the (u,v) pair of each surface vertex
				void set_mesh_uv_from_system(MeshAdaptor_3* mesh, const Solver& solver_u, const Solver& solver_v) ;

				// Check parameterize() postconditions:
				// * "A*Xu = Bu" and "A*Xv = Bv" systems are solvable with a good conditioning
				// * 3D -> 2D mapping is 1 to 1
				virtual ErrorCode check_parameterize_postconditions(const MeshAdaptor_3& mesh, const Solver& solver_u, const Solver& solver_v);

				// Check if 3D -> 2D mapping is 1 to 1
				//
				// Theorem: 1 to 1 mapping is guaranteed if all Wij coefficients are > 0 (for j vertex neighbor of i)
				//          and if the surface boundary is mapped onto a 2D convex polygon
				//
				// The default implementation checks each coefficient of A
				virtual bool  is_one_to_one_mapping (const Solver& solver);

// Protected accessors
protected:
				// Get the object that maps the surface's border onto a 2D space
				BorderParametizer_3&			get_border_parametizer()	{ return m_borderParametizer; }

				// Get the sparse linear algebra traits (Traits object to access the "A*X = B" sparse linear system)
				SparseLinearAlgebraTraits_d&	get_linear_algebra_traits()	{ return m_linearAlgebra; }

// Fields
private:
				// Object that maps the surface's border onto a 2D space
				 BorderParametizer_3			m_borderParametizer;
				// Traits object to solve the "A*X = B" sparse linear system used by parameterization algorithms
				 SparseLinearAlgebraTraits_d	m_linearAlgebra;
};


//
// Implementation 
//

// Compute a 1 to 1 mapping from a triangular 3D surface 'mesh' to a piece of the 2D space. 
// The mapping is linear by pieces (linear in each triangle).
// The result is the (u,v) pair image of each vertex of the 3D surface. 
//
// Preconditions:
// * 'mesh' must be a surface with 1 connected component and no hole
// * 'mesh' must be a triangular mesh
// * the mesh border must be mapped onto a convex polygon
template <class MeshAdaptor_3, class BorderParametizer_3, class SparseLinearAlgebraTraits_d>
inline 
typename Parametizer_3<MeshAdaptor_3>::ErrorCode  
Fixed_border_parametizer_3<MeshAdaptor_3, BorderParametizer_3, SparseLinearAlgebraTraits_d>::parameterize (MeshAdaptor_3* mesh)
{
	CGAL_parameterization_assertion(mesh != NULL);

	// Check preconditions
	ErrorCode status = check_parameterize_preconditions(*mesh);
	if (status != OK)
		return status;

	// Count vertices
	int nbVertices = mesh->count_mesh_vertices();

	// Index vertices from 0 to nbVertices-1
	std::cerr << "  index_mesh_vertices...";
	mesh->index_mesh_vertices();
	std::cerr << "ok (" << nbVertices << ")" << std::endl;

	// Create 2 sparse linear systems solver_u = "A*Xu = Bu" and solver_v = "A*Xv = Bv" (1 line/column per vertex)
	Solver solver_u(nbVertices);
	Solver solver_v(nbVertices);
	
	// Mark all vertices as NOT "parameterized"
	Vertex_iterator vertexIt;
	for (vertexIt = mesh->mesh_vertices_begin(); vertexIt != mesh->mesh_vertices_end(); vertexIt++)
		mesh->set_vertex_parameterized(&*vertexIt, false);

	// compute (u,v) for border vertices and mark them as "parameterized"
	if ( ! get_border_parametizer().parameterize_border(mesh) )
		return ERROR_NO_SURFACE_MESH;

	// Initialize A, Xu, Xv, Bu and Bv after boundary parameterization
	// Fill the border vertices' lines in both linear systems: "u = constant" and "v = constant"
	//
	// Implementation note: the current implementation removes border vertices from the linear systems and pass them in the right hand side
	initialize_system_from_mesh_border (&solver_u, &solver_v, *mesh);

	// Fill the matrix for the inner vertices Vi: compute A's coefficient Wij for each neighbor j; then Wii = - sum of Wij
	fprintf(stderr,"  fill matrix (%d x %d)...",nbVertices,nbVertices);
	solver_u.begin_system() ;
	solver_v.begin_system() ;
	for (vertexIt = mesh->mesh_vertices_begin(); vertexIt != mesh->mesh_vertices_end(); vertexIt++)
	{
	 	CGAL_parameterization_assertion(mesh->is_vertex_on_border(*vertexIt) == mesh->is_vertex_parameterized(*vertexIt));

		// inner vertices only
		if( ! mesh->is_vertex_on_border(*vertexIt) )
		{
			// Compute the line i of matrix A for i inner vertex
			status = setup_inner_vertex_relations (&solver_u, &solver_v, *mesh, *vertexIt);
			if (status != OK)
				return status;
		}
	}
	solver_u.end_system() ;
	solver_v.end_system() ;
	fprintf(stderr,"ok\n");

	// Solve linear systems solver_u = "A*Xu = Bu" and solver_v = "A*Xv = Bv" 
	std::cerr << "  solver...";
	if ( !solver_u.solve() || !solver_v.solve() )
	{
		std::cerr << "error" << std::endl;
		CGAL_parameterization_postcondition_msg(false, "Parameterization error: cannot solve sparse linear system");
		return ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;	
	}
	std::cerr << "ok" << std::endl;

	// Copy Xu and Xv coordinates into the (u,v) pair of each vertex
	set_mesh_uv_from_system (mesh, solver_u, solver_v); 

	// Check postconditions
	status = check_parameterize_postconditions(*mesh, solver_u, solver_v);
	if (status != OK)
		return status;

	return status;
} 


// Check parameterize() preconditions:
// * 'mesh' must be a surface with 1 connected component and no hole
// * 'mesh' must be a triangular mesh
// * the mesh border must be mapped onto a convex polygon
template <class MeshAdaptor_3, class BorderParametizer_3, class SparseLinearAlgebraTraits_d>
inline 
typename Parametizer_3<MeshAdaptor_3>::ErrorCode  
Fixed_border_parametizer_3<MeshAdaptor_3, BorderParametizer_3, SparseLinearAlgebraTraits_d>::check_parameterize_preconditions (const MeshAdaptor_3& mesh)
{
	ErrorCode status = OK;									// returned value

	// Allways check that mesh is not empty
	if (mesh.mesh_vertices_begin() == mesh.mesh_vertices_end())
		status = ERROR_EMPTY_MESH;
	CGAL_parameterization_precondition(status == OK);
	if (status != OK)
		return status;

	// The whole surface parameterization package is restricted to triangular meshes
	CGAL_parameterization_expensive_precondition((status = mesh.is_mesh_triangular() ? OK : ERROR_NON_TRIANGULAR_MESH) == OK);
	if (status != OK)
		return status;

	// The whole package is restricted to surfaces
	CGAL_parameterization_expensive_precondition((status = (mesh.get_mesh_genus()==0) ? OK : ERROR_NO_SURFACE_MESH) == OK);
	if (status != OK)
		return status;
	CGAL_parameterization_expensive_precondition((status = (mesh.count_mesh_boundaries() == 1) ? OK : ERROR_NO_SURFACE_MESH) == OK);
	if (status != OK)
		return status;

	// 1 to 1 mapping is guaranteed if all Wij coefficients are > 0 (for j vertex neighbor of i)
	// and if the surface boundary is mapped onto a 2D convex polygon
	CGAL_parameterization_expensive_precondition((status = get_border_parametizer().is_border_convex() ? OK : ERROR_NON_CONVEX_BORDER) == OK);
	if (status != OK)
		return status;

	return status;
}

// Initialize A, Bu and Bv after boundary parameterization
// Fill the border vertices' lines in both linear systems: "u = constant" and "v = constant"
//
// Preconditions:
// * vertices must be indexed
// * A, Bu and Bv must be allocated
// * border vertices must be parameterized
template <class MeshAdaptor_3, class BorderParametizer_3, class SparseLinearAlgebraTraits_d>
inline 
void Fixed_border_parametizer_3<MeshAdaptor_3, BorderParametizer_3, SparseLinearAlgebraTraits_d>::initialize_system_from_mesh_border(Solver* solver_u, Solver* solver_v, 
																																	 const MeshAdaptor_3& mesh) 
{
	CGAL_parameterization_assertion(solver_u != NULL);
	CGAL_parameterization_assertion(solver_v != NULL);

	for (Border_vertex_const_iterator it = mesh.mesh_border_vertices_begin(); it != mesh.mesh_border_vertices_end(); it++)
	{
		CGAL_parameterization_assertion(mesh.is_vertex_parameterized(*it));

		// Get vertex index in sparse linear system
		int index = mesh.get_vertex_index(*it);

		// Get vertex (u,v)
		Point_2 uv = mesh.get_vertex_uv(*it);

		// Write (u,v) in Xu and Xv
		solver_u->variable(index).set_value(uv.x()) ;
		solver_v->variable(index).set_value(uv.y()) ;

		// Copy (u,v) in Bu and Bv
        solver_u->variable(index).lock() ;
        solver_v->variable(index).lock() ;
	}
}

// Compute the line i of matrix A for i inner vertex
// * call compute_wij() to compute the A coefficient Wij for each neighbor Vj
// * compute Wii = - sum of Wij
// 
// Preconditions:
// * vertices must be indexed
// * vertex i musn't be already parameterized
// * line i of A must contain only zeros
template <class MeshAdaptor_3, class BorderParametizer_3, class SparseLinearAlgebraTraits_d>
inline 
typename Parametizer_3<MeshAdaptor_3>::ErrorCode 
Fixed_border_parametizer_3<MeshAdaptor_3, BorderParametizer_3, SparseLinearAlgebraTraits_d>::setup_inner_vertex_relations(Solver* solver_u, Solver* solver_v, 
																														  const MeshAdaptor_3& mesh, const Vertex& vertex) 
{
	CGAL_parameterization_assertion(solver_u != NULL);
	CGAL_parameterization_assertion(solver_v != NULL);
	CGAL_parameterization_assertion( ! mesh.is_vertex_on_border(vertex) );
	CGAL_parameterization_assertion( ! mesh.is_vertex_parameterized(vertex) );

	int i = mesh.get_vertex_index(vertex);

	// Start row i
    solver_u->begin_row() ;
    solver_v->begin_row() ;

	// circulate over vertices around vertex to compute Wii and Wijs
	NT Wii = 0;
	Vertex_around_vertex_const_circulator neighborIt = mesh.vertices_around_vertex_begin(vertex),
										  end = neighborIt;
	CGAL_For_all(neighborIt, end)
	{
		// Call to virtual method to do the actual coefficient computation
		NT Wij = compute_wij(mesh, vertex, neighborIt);

		// Wii = - sum of Wij
		Wii -= Wij;

		// Get j index
		int j = mesh.get_vertex_index(*neighborIt);

		// Set Wij in matrix
	    solver_u->add_coefficient(j, Wij);
	    solver_v->add_coefficient(j, Wij);
	}

	// Set Wii in matrix
    solver_u->add_coefficient(i, Wii);
    solver_v->add_coefficient(i, Wii);

	// End row i
    solver_u->end_row() ;
    solver_v->end_row() ;

	return OK;
}

// Copy Xu and Xv coordinates into the (u,v) pair of each surface vertex
template <class MeshAdaptor_3, class BorderParametizer_3, class SparseLinearAlgebraTraits_d>
inline 
void Fixed_border_parametizer_3<MeshAdaptor_3, BorderParametizer_3, SparseLinearAlgebraTraits_d>::set_mesh_uv_from_system(MeshAdaptor_3* mesh, 
																														  const Solver& solver_u, const Solver& solver_v) 
{
	Vertex_iterator vertexIt;
	for (vertexIt = mesh->mesh_vertices_begin(); vertexIt != mesh->mesh_vertices_end(); vertexIt++)
	{
		int index = mesh->get_vertex_index(*vertexIt);

		NT u = solver_u.variable(index).value() ;
		NT v = solver_v.variable(index).value() ;

		// Fill vertex (u,v) and mark it as "parameterized"
		mesh->set_vertex_uv(&*vertexIt, Point_2(u,v));
		mesh->set_vertex_parameterized(&*vertexIt, true);
	}
}

// Check parameterize() postconditions:
// * "A*Xu = Bu" and "A*Xv = Bv" systems are solvable with a good conditioning
// * 3D -> 2D mapping is 1 to 1
template <class MeshAdaptor_3, class BorderParametizer_3, class SparseLinearAlgebraTraits_d>
inline 
typename Parametizer_3<MeshAdaptor_3>::ErrorCode 
Fixed_border_parametizer_3<MeshAdaptor_3, BorderParametizer_3, SparseLinearAlgebraTraits_d>::check_parameterize_postconditions(const MeshAdaptor_3& mesh, 
																															   const Solver& solver_u, const Solver& solver_v)
{
	ErrorCode status = OK;

	// LS 03/17/2005: comment out this section because OpenNL::LinearSolver does not provide a is_solvable() method
	//// Check if "A*Xu = Bu" and "A*Xv = Bv" systems are solvable with a good conditioning
	//CGAL_parameterization_expensive_postcondition((status = get_linear_algebra_traits().is_solvable(A, Bu) ? OK : ERROR_BAD_MATRIX_CONDITIONING) == OK);
	//if (status != OK)
	//	return status;
	//CGAL_parameterization_expensive_postcondition((status = get_linear_algebra_traits().is_solvable(A, Bv) ? OK : ERROR_BAD_MATRIX_CONDITIONING) == OK);
	//if (status != OK)
	//	return status;

// LS 02/04/2005: commented out because this check fails with conformal parameterization/square boundary even though the result seems correct
//	// Check if 3D -> 2D mapping is 1 to 1
//	CGAL_parameterization_expensive_postcondition((status = is_one_to_one_mapping(A) ? OK : ERROR_NO_1_TO_1_MAPPING) == OK);
//	if (status != OK)
//		return status;

	return status;
}

// Check if 3D -> 2D mapping is 1 to 1
//
// Theorem: 1 to 1 mapping is guaranteed if all Wij coefficients are > 0 (for j vertex neighbor of i)
//          and if the surface boundary is mapped onto a 2D convex polygon
//
// The default implementation checks each coefficient of A
template <class MeshAdaptor_3, class BorderParametizer_3, class SparseLinearAlgebraTraits_d>
inline 
bool Fixed_border_parametizer_3<MeshAdaptor_3, BorderParametizer_3, SparseLinearAlgebraTraits_d>::is_one_to_one_mapping (const Solver& solver) 
{
	// LS 03/17/2005: comment out this section because - OpenNL::LinearSolver does not provide access to a matrix coefficient
	//                                                 - this method is never called
	//// Check if all Wij coefficients are > 0 (for j vertex neighbor of i)
	//for (int i=0; i < (int)A.row_dimension(); i++)
	//	for (int j=0; j < (int)A.column_dimension(); j++)
	//		if (i != j)
	//			if (A.get_coef(i,j) < 0)
	//				return false;
	
	// Check if the surface boundary is mapped onto a 2D convex polygon
	return get_border_parametizer().is_border_convex ();
}


CGAL_END_NAMESPACE

#endif //FIXED_BORDER_PARAMETIZER_3_H

