// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
// Author(s)     : Monique Teillaud, Sylvain Pion

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_GLOBAL_FUNCTIONS_ON_ROOT_FOR_CIRCLE_2_2_H
#define CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_GLOBAL_FUNCTIONS_ON_ROOT_FOR_CIRCLE_2_2_H

namespace CGAL {

template < class AK >
inline 
Comparison_result 
compare_x(const typename AK::Root_for_circles_2_2& r1,
	   const typename AK::Root_for_circles_2_2& r2)
{ return AK().compare_x_object()(r1, r2); }

template < class AK >
inline 
Comparison_result 
compare_y(const typename AK::Root_for_circles_2_2& r1,
	   const typename AK::Root_for_circles_2_2& r2)
{ return AK().compare_y_object()(r1, r2); }

template < class AK >
inline 
Comparison_result 
compare_xy(const typename AK::Root_for_circles_2_2& r1,
	     const typename AK::Root_for_circles_2_2& r2)
{ return AK().compare_xy_object()(r1, r2); }

} //namespace CGAL

#endif //CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_GLOBAL_FUNCTIONS_ON_ROOT_FOR_CIRCLE_2_2_H
