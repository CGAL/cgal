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
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion
//             Pedro Machado
//             Julien Hazebrouck
//             Damien Leroy

#ifndef CGAL_ALGEBRAIC_KERNEL_GLOBAL_FUNCTIONS_ON_ROOT_FOR_SPHERE_2_3_H
#define CGAL_ALGEBRAIC_KERNEL_GLOBAL_FUNCTIONS_ON_ROOT_FOR_SPHERE_2_3_H

#include <CGAL/license/Circular_kernel_3.h>


#include <CGAL/enum.h>

namespace CGAL {

template < class AK >
inline 
Comparison_result 
compare_x(const typename AK::Root_for_spheres_2_3& r1,
	   const typename AK::Root_for_spheres_2_3& r2)
{ return AK().compare_x_object()(r1, r2); }

template < class AK >
inline 
Comparison_result 
compare_y(const typename AK::Root_for_spheres_2_3& r1,
	   const typename AK::Root_for_spheres_2_3& r2)
{ return AK().compare_y_object()(r1, r2); }

template < class AK >
inline 
Comparison_result 
compare_z(const typename AK::Root_for_spheres_2_3& r1,
	     const typename AK::Root_for_spheres_2_3& r2)
{ return AK().compare_z_object()(r1, r2); }

template < class AK >
inline 
Comparison_result 
compare_xy(const typename AK::Root_for_spheres_2_3& r1,
	     const typename AK::Root_for_spheres_2_3& r2)
{ return AK().compare_xy_object()(r1, r2); }

template < class AK >
inline 
Comparison_result 
compare_xyz(const typename AK::Root_for_spheres_2_3& r1,
	     const typename AK::Root_for_spheres_2_3& r2)
{ return AK().compare_xyz_object()(r1, r2); }

} //namespace CGAL

#endif //CGAL_ALGEBRAIC_KERNEL_GLOBAL_FUNCTIONS_ON_ROOT_FOR_SPHERE_2_3_H
