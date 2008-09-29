// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_KERNEL_GLOBAL_FUNCTIONS_ON_LINE_2_H
#define CGAL_CIRCULAR_KERNEL_GLOBAL_FUNCTIONS_ON_LINE_2_H

CGAL_BEGIN_NAMESPACE

template< class CK >
inline
typename CK::Polynomial_1_2
get_equation(const Line_2<CK> & l)
{
  return CK().get_equation_object()(l);
}

template< class CK >
inline
CGAL::Line_2<CK>
construct_line_2(const typename CK::Polynomial_1_2 & eq)
{
  return CK().construct_line_2_object()(eq);
}

template< class CK, class OutputIterator>
inline
OutputIterator
intersect_2( const Line_2<CK> & l,
			   const Circle_2<CK> & c,
			   OutputIterator res )
{
  return CK().intersect_2_object()(l,c,res);
}


CGAL_END_NAMESPACE

#endif // CGAL_CIRCULAR_KERNEL_GLOBAL_FUNCTIONS_ON_LINE_2_H
