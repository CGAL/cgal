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
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_ARC_2_H
#define CGAL_CIRCULAR_ARC_2_H

#include <CGAL/config.h>

namespace CGAL {
  
template <class CircularKernel> 
class Circular_arc_2 
  : public CircularKernel::Kernel_base::Circular_arc_2
{
  typedef typename CircularKernel::RT             RT;
  typedef typename CircularKernel::FT             FT;
  typedef typename CircularKernel::Point_2        Point_2;
  typedef typename CircularKernel::Line_2         Line_2;
  typedef typename CircularKernel::Circle_2       Circle_2;
  typedef typename CircularKernel::Circular_arc_point_2
                                                Circular_arc_point_2;
  
  typedef typename CircularKernel::Kernel_base::Circular_arc_2 RCircular_arc_2; 
  // RCircular_arc_2 to avoid clash with self 
public:
  typedef  RCircular_arc_2 Rep;
  typedef  CircularKernel   R; 
  

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }


  Circular_arc_2()
    : RCircular_arc_2(typename R::Construct_circular_arc_2()())
  {}

  Circular_arc_2(const Circle_2 &c)
    : RCircular_arc_2(typename R::Construct_circular_arc_2()(c))
  {}

  // Not Documented
  Circular_arc_2(const Circle_2 &support, 
                 const Line_2 &l1, const bool b_l1,
                 const Line_2 &l2, const bool b_l2)
    : RCircular_arc_2(typename 
		      R::Construct_circular_arc_2()(support,l1,b_l1,l2,b_l2))
  {}

  // Not Documented
  Circular_arc_2(const Circle_2 &c, 
		 const Circle_2 &c1, const bool b_1,
		 const Circle_2 &c2, const bool b_2)
    : RCircular_arc_2(typename 
		      R::Construct_circular_arc_2()(c,c1,b_1,c2,b_2))
  {}

  Circular_arc_2(const Point_2 &start,
                 const Point_2 &middle,
                 const Point_2 &end)
    : RCircular_arc_2(typename 
		      R::Construct_circular_arc_2()(start, middle, end)) 
  {}
  
  Circular_arc_2(const Circle_2 &support,
                 const Circular_arc_point_2 &begin,
                 const Circular_arc_point_2 &end)
    : RCircular_arc_2(typename 
		      R::Construct_circular_arc_2()(support, begin, end)) 
  {}

  Circular_arc_2(const Point_2 &start,
                 const Point_2 &end,
		 const FT &bulge)
    : RCircular_arc_2(typename 
		      R::Construct_circular_arc_2()(start, end, bulge)) 
  {}
  
 Circular_arc_2(const RCircular_arc_2 & a)
    : RCircular_arc_2(a)
  {}


  typename cpp11::result_of<typename R::Construct_circular_source_vertex_2(Circular_arc_2)>::type
  source() const
  {
    return typename R::Construct_circular_source_vertex_2()(*this);
  }

  typename cpp11::result_of<typename R::Construct_circular_target_vertex_2(Circular_arc_2)>::type
  target() const
  {
    return typename R::Construct_circular_target_vertex_2()(*this);
  }

  typename cpp11::result_of<typename R::Construct_circular_min_vertex_2(Circular_arc_2)>::type
  left() const
  {
    return typename R::Construct_circular_min_vertex_2()(*this);
  }

  typename cpp11::result_of<typename R::Construct_circular_max_vertex_2(Circular_arc_2)>::type
  right() const
  {
    return typename R::Construct_circular_max_vertex_2()(*this);
  }

  bool is_x_monotone() const
  {
    return typename R::Is_x_monotone_2()(*this);
  }

  bool is_y_monotone() const
  {
    return typename R::Is_y_monotone_2()(*this);
  }

  typename cpp11::result_of<typename R::Construct_circle_2(Circular_arc_2)>::type
  supporting_circle() const
  {
    return typename R::Construct_circle_2()(*this);
  }

  typename cpp11::result_of<typename R::Construct_center_2(Circular_arc_2)>::type
  center() const
  {
    return typename R::Construct_center_2()(*this);
  }

  typename cpp11::result_of<typename R::Compute_squared_radius_2( Circular_arc_2)>::type
  squared_radius() const
  {
    return typename R::Compute_squared_radius_2()(*this);
  }

  Bbox_2 bbox(void) const
  {
    return typename R::Construct_bbox_2()(*this);
  }

};

  template < typename CircularKernel >
  inline
  bool
  operator==(const Circular_arc_2<CircularKernel> &p,
	     const Circular_arc_2<CircularKernel> &q)
  {
    return CircularKernel().equal_2_object()(p, q);
  }
  
  template < typename CircularKernel >
  inline
  bool
  operator!=(const Circular_arc_2<CircularKernel> &p,
	     const Circular_arc_2<CircularKernel> &q)
  {
    return ! (p == q);
  }

  template < typename CK >
  std::ostream &
  operator<<(std::ostream & os, const Circular_arc_2<CK> &a)
  {
    // The output format is :
    // - supporting circle
    // - circle c1
    // - bool b1
    // - circle c2
    // - bool b2
    return os << a.supporting_circle() << " "
	      << a.source() << " "
	      << a.target() << " ";
  }
  
  template < typename CK >
  std::istream &
  operator>>(std::istream & is, Circular_arc_2<CK> &a)
  {
    typename CK::Circle_2 s;
    typename CK::Circular_arc_point_2 p1;
    typename CK::Circular_arc_point_2 p2;
    is >> s >> p1 >> p2 ;
    if (is)
      a = Circular_arc_2<CK>(s, p1, p2);
    return is;
  }


} //namespace CGAL

#endif // CGAL_CIRCULAR_ARC_2_H
