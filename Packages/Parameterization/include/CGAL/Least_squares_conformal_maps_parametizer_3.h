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


#ifndef LEAST_SQUARES_CONFORMAL_MAPS_PARAMETIZER_3_H
#define LEAST_SQUARES_CONFORMAL_MAPS_PARAMETIZER_3_H

#include "Parametizer_3.h"

CGAL_BEGIN_NAMESPACE


// Class Least_squares_conformal_maps_parametizer_3
// Model of the Parametizer_3 concept
// Implement Least Square Conformal Maps parameterization (Levy et al)
// No need to map the surface's border onto a convex polygon but 1 to 1 mapping not guaranteed.
// This is a conformal parameterization, i.e. it attempts to preserve angles. 
// 
class Least_squares_conformal_maps_parametizer_3 : public Parametizer_3 {
// Public stuff
public:
				// Operations
				// Compute a 1 to 1 mapping from a triangular 3D surface 'mesh' to a piece of the 2D space. The mapping is linear by pieces (linear in each triangle).
				// The result is the (u,v) pair image of each vertex of the 3D surface. 
				// 
				// Preconditions:
				// * 'mesh' must be a surface with 1 connected component and no hole
				// * 'mesh' must be a triangular mesh
				ErrorCode  parameterize (MeshAdaptor_3* mesh) {
								
				}
// Protected stuff
protected:
				// Operations
// Private stuff
private:
				// Operations
};


CGAL_END_NAMESPACE

#endif //LEAST_SQUARES_CONFORMAL_MAPS_PARAMETIZER_3_H

