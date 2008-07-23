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

#ifndef CGAL_CIRCULAR_KERNEL_FUNCTIONS_ON_CIRCLE_2_H
#define CGAL_CIRCULAR_KERNEL_FUNCTIONS_ON_CIRCLE_2_H

namespace CGAL {

  // Should we have an iterator based interface, or both ?
  template <class CK>
  typename CK::Circular_arc_point_2
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
  typename CK::Circular_arc_point_2
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

} // namespace CGAL

#endif // CGAL_CIRCULAR_KERNEL_FUNCTIONS_ON_CIRCLE_2_H
