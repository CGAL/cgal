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


#ifndef CGAL_PARAMETERIZATION_H
#define CGAL_PARAMETERIZATION_H

#include <CGAL/Mean_value_coordinates_parametizer_3.h>

CGAL_BEGIN_NAMESPACE


//
// Main entry points of the parameterization header as C functions
//

// Compute a 1 to 1 mapping from a triangular 3D surface 'mesh' to 2D circle, 
// using Floater's mean value coordinates algorithm .
// 1 to 1 mapping is guaranteed.
//
// The mapping is linear by pieces (linear in each triangle).
// The result is the (u,v) pair image of each vertex of the 3D surface. 
//
// Preconditions:
// * 'mesh' must be a surface with 1 connected component and no hole
// * 'mesh' must be a triangular mesh
// * the mesh border must be mapped onto a convex polygon
//
template <class MeshAdaptor_3>			
typename Parametizer_3<MeshAdaptor_3>::ErrorCode 
parameterize(MeshAdaptor_3* mesh)	// 3D mesh, model of MeshAdaptor_3 concept
{
    Mean_value_coordinates_parametizer_3<MeshAdaptor_3> parametizer;
    return parametizer.parameterize(mesh);
}

// Compute a 1 to 1 mapping from a triangular 3D surface 'mesh' 
// to a piece of the 2D space. 
// The mapping is linear by pieces (linear in each triangle).
// The result is the (u,v) pair image of each vertex of the 3D surface. 
//
// 1 to 1 mapping may be guaranteed or not, depending of 
// Parametizer_3_ algorithm chosen
//
// Preconditions:
// * 'mesh' must be a surface with 1 connected component and no hole
// * 'mesh' must be a triangular mesh
// * the mesh border must be mapped onto a convex polygon 
// (for fixed border parameterizations)
//
template <class MeshAdaptor_3, class Parametizer_3_>			
typename Parametizer_3<MeshAdaptor_3>::ErrorCode 
parameterize(MeshAdaptor_3* mesh,		// 3D mesh, model of MeshAdaptor_3 concept
			 Parametizer_3_ parametizer)// Parameterization class for MeshAdaptor_3 mesh
{
    return parametizer.parameterize(mesh);
}


CGAL_END_NAMESPACE

#endif //CGAL_PARAMETERIZATION_H

