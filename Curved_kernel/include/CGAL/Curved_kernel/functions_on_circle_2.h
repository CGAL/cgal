// Copyright (c) 2005  INRIA Sophia-Antipolis (France) 
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CURVED_KERNEL_FUNCTIONS_ON_CIRCLE_2_H
#define CGAL_CURVED_KERNEL_FUNCTIONS_ON_CIRCLE_2_H

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

#endif // CGAL_CURVED_KERNEL_FUNCTIONS_ON_CIRCLE_2_H
