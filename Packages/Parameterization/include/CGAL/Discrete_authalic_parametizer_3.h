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


#ifndef DISCRETE_AUTHALIC_PARAMETIZER_3_H
#define DISCRETE_AUTHALIC_PARAMETIZER_3_H

#include "CGAL/Fixed_border_parametizer_3.h"

CGAL_BEGIN_NAMESPACE


// Class Discrete_authalic_parametizer_3
// Model of the Parametizer_3 concept.
// Implement Discrete Authalic Parameterization algorithm (Alliez et al). 1 to 1 mapping is guaranteed if surface's border is mapped onto a convex polygon.
// This is an authalic parameterization, i.e. it attempts to preserve areas.
template <class MeshAdaptor_3,														// 3D surface
		  class BorderParametizer_3 = Circular_border_parametizer_3<MeshAdaptor_3>,	// Class to map the surface's border onto a 2D space
		  class SparseLinearAlgebraTraits_d = OpenNL::DefaultLinearSolverTraits<MeshAdaptor_3::NT> >	
																					// Traits class for solving a sparse linear system "A*X = B"
class Discrete_authalic_parametizer_3 
	: public Fixed_border_parametizer_3<MeshAdaptor_3, BorderParametizer_3, SparseLinearAlgebraTraits_d>
{
// Public stuff
public:
				// Constructor
				// @param borderParametizer	Object that maps the surface's border onto a 2D space
				// @param linearAlgebra		Traits object to access the "A*X = B" sparse linear system used by parameterization algorithms
				Discrete_authalic_parametizer_3 (BorderParametizer_3 borderParametizer = BorderParametizer_3(), 
											     SparseLinearAlgebraTraits_d linearAlgebra = SparseLinearAlgebraTraits_d()) 
				:	Fixed_border_parametizer_3<MeshAdaptor_3, BorderParametizer_3, SparseLinearAlgebraTraits_d>(borderParametizer, linearAlgebra)
				{}

				// Default copy constructor and operator =() are fine

// Protected stuff
protected:
				// compute wij = (i,j) coefficient of matrix A for j neighbor vertex of i
				virtual NT  compute_wij(MeshAdaptor_3& mesh, Vertex& main_vertex_Vi, Vertex_around_vertex_circulator neighbor_vertex_Vj) 
				{
					Point_3	position_Vi = mesh.get_vertex_position(main_vertex_Vi);
					Point_3	position_Vj = mesh.get_vertex_position(*neighbor_vertex_Vj);

					// Compute the square norm of Vj -> Vi vector
					Vector_3 edge = position_Vi - position_Vj;
					double square_len = edge*edge;

					// Compute cotangent of corner specified by Vk,Vj,Vi points (ie cotan of Vj corner)
					// if Vk is the vertex before Vj when circulating around Vi
					Vertex_around_vertex_circulator previous_vertex_Vk = neighbor_vertex_Vj; previous_vertex_Vk --;
					Point_3	position_Vk = mesh.get_vertex_position(*previous_vertex_Vk);
					double cotg_psi_ij  = cotangent(position_Vk, position_Vj, position_Vi);

					// Compute cotangent of corner specified by Vi,Vj,Vl points (ie cotan of Vj corner)
					// if Vl is the vertex after Vj when circulating around Vi
					Vertex_around_vertex_circulator next_vertex_Vl = neighbor_vertex_Vj; next_vertex_Vl ++;
					Point_3	position_Vl = mesh.get_vertex_position(*next_vertex_Vl);
					double cotg_theta_ij = cotangent(position_Vi, position_Vj, position_Vl);

					double weight = 0.0;
					assert(square_len != 0.0);										// 2 points are identical!
					if(square_len != 0.0)
						weight = (cotg_psi_ij+cotg_theta_ij)/square_len;

					return weight;
				}
};


CGAL_END_NAMESPACE

#endif //DISCRETE_AUTHALIC_PARAMETIZER_3_H

