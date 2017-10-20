// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion
//             Pedro Machado
//             Julien Hazebrouck
//             Damien Leroy

#ifndef CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_3_H
#define CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_3_H

#include <CGAL/license/Circular_kernel_3.h>


namespace CGAL {

template < class AK >
inline 
typename AK::Polynomial_for_spheres_2_3
construct_polynomial_sphere_2_3(const typename AK::FT& xc,
			       const typename AK::FT& yc,
			       const typename AK::FT& zc,
			       const typename AK::FT& r_sq)
{ return AK().construct_polynomial_sphere_2_3_object()(xc, yc, zc, r_sq); }

} //namespace CGAL

#endif // CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_3_H
