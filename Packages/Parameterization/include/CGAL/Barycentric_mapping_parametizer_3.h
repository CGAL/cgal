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


#ifndef BARYCENTRIC_MAPPING_PARAMETIZER_3_H
#define BARYCENTRIC_MAPPING_PARAMETIZER_3_H

#include "CGAL/Fixed_border_parametizer_3.h"
#include "CGAL/parameterization_assertions.h"

CGAL_BEGIN_NAMESPACE


// Class Barycentric_mapping_parametizer_3
// Model of the Parametizer_3 concept
// Implement Tutte's barycentric mapping. 1 to 1 mapping is guaranteed if surface's border is mapped onto a convex polygon.
template <class MeshAdaptor_3,														// 3D surface
		  class BorderParametizer_3 = Circular_border_parametizer_3<MeshAdaptor_3>,	// Class to map the surface's border onto a 2D space
		  class SparseLinearAlgebraTraits_d = OpenNL::DefaultLinearSolverTraits<typename MeshAdaptor_3::NT> >	
																					// Traits class for solving a sparse linear system "A*X = B"
class Barycentric_mapping_parametizer_3 
	: public Fixed_border_parametizer_3<MeshAdaptor_3, BorderParametizer_3, SparseLinearAlgebraTraits_d>
{
// Public types
public:
				// Export Mesh_Adaptor_3, BorderParametizer_3 and SparseLinearAlgebraTraits_d types and subtypes
				typedef MeshAdaptor_3											Mesh_adaptor_3;
				typedef typename Parametizer_3<MeshAdaptor_3>::ErrorCode		ErrorCode;
				typedef typename MeshAdaptor_3::NT								NT;
				typedef typename MeshAdaptor_3::Face							Face;
				typedef typename MeshAdaptor_3::Vertex							Vertex;
				typedef typename MeshAdaptor_3::Point_3							Point_3;
				typedef typename MeshAdaptor_3::Point_2							Point_2;
				typedef typename MeshAdaptor_3::Vector_3						Vector_3;
				typedef typename MeshAdaptor_3::Face_iterator					Face_iterator;
				typedef typename MeshAdaptor_3::Vertex_iterator					Vertex_iterator;
				typedef typename MeshAdaptor_3::Border_vertex_iterator			Border_vertex_iterator;
				typedef typename MeshAdaptor_3::Vertex_around_face_circulator	Vertex_around_face_circulator;
				typedef typename MeshAdaptor_3::Vertex_around_vertex_circulator	Vertex_around_vertex_circulator;
				typedef BorderParametizer_3										Border_parametizer_3;
				typedef SparseLinearAlgebraTraits_d								Sparse_linear_algebra_traits_d;
				typedef typename SparseLinearAlgebraTraits_d::Vector			Vector;
				typedef typename SparseLinearAlgebraTraits_d::Matrix			Matrix;

// Public operations
public:
				// Constructor
				// @param borderParametizer	Object that maps the surface's border onto a 2D space
				// @param linearAlgebra		Traits object to access the "A*X = B" sparse linear system used by parameterization algorithms
				Barycentric_mapping_parametizer_3 (BorderParametizer_3 borderParametizer = BorderParametizer_3(), 
												   SparseLinearAlgebraTraits_d linearAlgebra = SparseLinearAlgebraTraits_d()) 
				:	Fixed_border_parametizer_3<MeshAdaptor_3, BorderParametizer_3, SparseLinearAlgebraTraits_d>(borderParametizer, linearAlgebra)
				{}

				// Default copy constructor and operator =() are fine

// Protected stuff
protected:
				// compute wij = (i,j) coefficient of matrix A for j neighbor vertex of i 
				virtual NT  compute_wij(MeshAdaptor_3& mesh, Vertex& main_vertex_Vi, Vertex_around_vertex_circulator neighbor_vertex_Vj) 
				{
					// Tutte algorithm is the most simple one: Wij = 1 for j neighbor vertex of i
					return 1;
				}

				// Check if 3D -> 2D mapping is 1 to 1
				//
				// Theorem: 1 to 1 mapping is guaranteed if all Wij coefficients are > 0 (for j vertex neighbor of i)
				//          and if the surface boundary is mapped onto a 2D convex polygon
				virtual bool  is_one_to_one_mapping (const Matrix& A)
				{
					// All Wij coefficients = 1 (for j vertex neighbor of i), thus  
					// mapping is guaranteed if the surface boundary is mapped onto a 2D convex polygon
					return get_border_parametizer().is_border_convex ();
				}
};


CGAL_END_NAMESPACE

#endif //BARYCENTRIC_MAPPING_PARAMETIZER_3_H

