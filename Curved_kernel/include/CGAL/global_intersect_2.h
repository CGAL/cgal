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

#ifndef CGAL_CURVED_KERNEL_GLOBAL_INTERSECT_2
#define CGAL_CURVED_KERNEL_GLOBAL_INTERSECT_2

#include <CGAL/Curved_kernel/function_objects_on_line_2.h>
#include <CGAL/Curved_kernel/function_objects_on_circle_2.h>
#include <CGAL/Curved_kernel/function_objects_polynomial_circular.h>

CGAL_BEGIN_NAMESPACE

// missing: intersection of Line_2 with {Circle_2,Circular_arc_2}

template< class CK, class OutputIterator>
inline
OutputIterator
intersect( const typename CK::Circle_2 & c1,
	   const typename CK::Circle_2 & c2,
	   OutputIterator res )
{
  return CK::Intersect_2()(c1,c2,res);
}

template< class CK, class OutputIterator>
inline
OutputIterator
intersect( const typename CK::Line_arc_2 & c1,
	   const typename CK::Line_arc_2 & c2,
	   OutputIterator res )
{
  return CK::Intersect_2()(c1,c2,res);
}

template< class CK, class OutputIterator>
inline
OutputIterator
intersect( const typename CK::Circular_arc_2 & c1,
	   const typename CK::Circular_arc_2 & c2,
	   OutputIterator res )
{
  return CK::Intersect_2()(c1,c2,res);
}


template< class CK, class OutputIterator>
inline
OutputIterator
intersect( const typename CK::Line_arc_2 & c1,
	   const typename CK::Circle_2 & c2,
	   OutputIterator res )
{
  return CK::Intersect_2()(c1,c2,res);
}
template< class CK, class OutputIterator>
inline
OutputIterator
intersect( const typename CK::Circle_2 & c1,
	   const typename CK::Line_arc_2 & c2,
	   OutputIterator res )
{
  return CK::Intersect_2()(c1,c2,res);
}

template< class CK, class OutputIterator>
inline
OutputIterator
intersect( const typename CK::Line_arc_2 & c1,
	   const typename CK::Circular_arc_2 & c2,
	   OutputIterator res )
{
  return CK::Intersect_2()(c1,c2,res);
}
template< class CK, class OutputIterator>
inline
OutputIterator
intersect( const typename CK::Circular_arc_2 & c1,
	   const typename CK::Line_arc_2 & c2,
	   OutputIterator res )
{
  return CK::Intersect_2()(c1,c2,res);
}

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_GLOBAL_INTERSECT_2
