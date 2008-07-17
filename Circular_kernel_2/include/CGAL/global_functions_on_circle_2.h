// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Julien Hazebrouck

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_KERNEL_GLOBAL_FUNCTIONS_ON_CIRCLE_2_H
#define CGAL_CIRCULAR_KERNEL_GLOBAL_FUNCTIONS_ON_CIRCLE_2_H

namespace CGAL {

template< class CK >
inline
typename CK::Polynomial_for_circles_2_2
get_equation(const typename CK::Circle_2 & c)
{
  return CK().get_equation_object()(c);
}

template< class CK >
inline
typename CK::Circle_2
construct_circle_2(const typename CK::Polynomial_for_circles_2_2 & eq)
{
  return CK().construct_circle_2_object()(eq);
}

template< class CK, class OutputIterator>
inline
OutputIterator
intersect_2( const typename CK::Circle_2 & c1,
			   const typename CK::Circle_2 & c2,
			   OutputIterator res )
{
  return CK().intersect_2_object()(c1,c2,res);
}

template< class CK >
inline
bool
has_on_2(const typename CK::Circle_2 &c, const typename CK::Circular_arc_point_2 &p)
{
  return CK().has_on_2_object()(c, p);
}

} // namespace CGAL

#endif // CGAL_CIRCULAR_KERNEL_GLOBAL_FUNCTIONS_ON_CIRCLE_2_H
