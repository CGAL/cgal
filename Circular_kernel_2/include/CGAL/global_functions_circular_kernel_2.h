// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Julien Hazebrouck, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)


#ifndef CGAL_CIRCULAR_KERNEL_GLOBAL_FUNCTIONS_CIRCULAR_KERNEL_2_H
#define CGAL_CIRCULAR_KERNEL_GLOBAL_FUNCTIONS_CIRCULAR_KERNEL_2_H

#include <CGAL/license/Circular_kernel_2.h>


// global functions

#include <CGAL/Circular_arc_2.h>
#include <CGAL/Circular_arc_point_2.h>
#include <CGAL/Line_arc_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Circular_kernel_2/internal_functions_on_circle_2.h>

namespace CGAL {

template <class CK>
Circular_arc_point_2<CK>
x_extremal_point(const Circle_2<CK> & c, bool i)
{
	return CircularFunctors::x_extremal_point<CK>(c,i);
}
  
template <class CK, class OutputIterator>
OutputIterator
x_extremal_points(const Circle_2<CK> & c, OutputIterator res)
{
 	return CircularFunctors::x_extremal_points<CK>(c,res);
}

template <class CK>
Circular_arc_point_2<CK>
y_extremal_point(const Circle_2<CK> & c, bool i)
{
  return CircularFunctors::y_extremal_point<CK>(c,i);
}

template <class CK, class OutputIterator>
OutputIterator
y_extremal_points(const Circle_2<CK> & c, OutputIterator res)
{
 	return CircularFunctors::y_extremal_points<CK>(c,res);
}

// Not Documented
template< class CK >
inline
CGAL::Comparison_result 
compare_x(const Circular_arc_2<CK> &A1, const bool b1, 
	  const Circular_arc_2<CK> &A2, const bool b2)
{
  return CK().compare_x_2_object()(A1, b1, A2, b2);
}

template< class CK >
inline
CGAL::Comparison_result 
compare_x(const Circular_arc_point_2<CK> &p, const Circular_arc_point_2<CK> &q)
{
  return CK().compare_x_2_object()(p, q);
}

template< class CK >
inline
CGAL::Comparison_result 
compare_y(const Circular_arc_point_2<CK> &p, const Circular_arc_point_2<CK> &q)
{
  return CK().compare_y_2_object()(p, q);
}

template< class CK >
inline
CGAL::Comparison_result 
compare_xy(const Circular_arc_point_2<CK> &p, const Circular_arc_point_2<CK> &q)
{
  return CK().compare_xy_2_object()(p, q);
}

template< class CK >
inline
CGAL::Comparison_result 
compare_y_to_right(const Circular_arc_2<CK> &A1,
		   const Circular_arc_2<CK> &A2,
		   const Circular_arc_point_2<CK> &p)
{
  return CK().compare_y_to_right_2_object()(A1, A2, p);
}

template < class CK >
inline
bool
has_in_x_range(const Circular_arc_2<CK> &A, const Circular_arc_point_2<CK> &p) 
{
  return CK().in_x_range_2_object()(A, p);
}

template < class CK >
CGAL::Comparison_result       
compare_y_at_x(const Circular_arc_point_2<CK> &p, const Circular_arc_2<CK> &a)
{
  return CK().compare_y_at_x_2_object()(p, a);
}

template < class CK, class OutputIterator >
OutputIterator
make_x_monotone(const Circular_arc_2<CK> &A, OutputIterator it)
{
  return CK().make_x_monotone_2_object()(A, it);
}

template < class CK, class OutputIterator >
OutputIterator
make_xy_monotone(const Circular_arc_2<CK> &A, OutputIterator it)
{
  return CK().make_xy_monotone_2_object()(A, it);
}

template< class CK >
inline
bool
has_on(const Circle_2<CK> &c, const Circular_arc_point_2<CK> &p)
{
  return CK().has_on_2_object()(c, p);
}

// Not Documented
template< class CK >
inline
CGAL::Comparison_result 
compare_y_to_right(const Line_arc_2<CK> &A1,
		   const Line_arc_2<CK> &A2,
		   const Circular_arc_point_2<CK> &p)
{
  return CK().compare_y_to_right_2_object()(A1, A2, p);
}

template < class CK >
inline
bool
has_in_x_range(const Line_arc_2<CK> &A, const Circular_arc_point_2<CK> &p) 
{
  return CK().in_x_range_2_object()(A, p);
}

template < class CK >
CGAL::Comparison_result       
compare_y_at_x(const Circular_arc_point_2<CK> &p, const Line_arc_2<CK> &a)
{
  return CK().compare_y_at_x_2_object()(p, a);
}

// Not Documented
template < class CK, class OutputIterator >
OutputIterator
make_x_monotone(const Line_arc_2<CK> &A, OutputIterator it)
{
  return CK().make_x_monotone_2_object()(A, it);
}

template < class CK, class OutputIterator >
OutputIterator
make_xy_monotone(const Line_arc_2<CK> &A, OutputIterator it)
{
  return CK().make_xy_monotone_2_object()(A, it);
}

} //namespace CGAL

#endif // CGAL_CIRCULAR_KERNEL_GLOBAL_FUNCTIONS_CIRCULAR_KERNEL_2_H
