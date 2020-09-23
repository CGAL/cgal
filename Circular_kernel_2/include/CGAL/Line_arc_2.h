// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Julien Hazebrouck, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_LINE_ARC_2_H
#define CGAL_LINE_ARC_2_H

#include <CGAL/license/Circular_kernel_2.h>


#include <CGAL/result_of.h>
#include <CGAL/Bbox_2.h>

namespace CGAL {

template <class CircularKernel>
class Line_arc_2
  : public CircularKernel::Kernel_base::Line_arc_2
{
  typedef typename CircularKernel::FT                        FT;
  typedef typename CircularKernel::RT                        RT;
  //typedef typename CircularKernel::Linear_kernel::Point_2    Point_2;
  typedef typename CircularKernel::Point_2                   Point_2;
  typedef typename CircularKernel::Line_2                    Line_2;
  typedef typename CircularKernel::Circle_2                  Circle_2;
  typedef typename CircularKernel::Circular_arc_point_2   Circular_arc_point_2;
  typedef typename CircularKernel::Segment_2                 Segment_2;

  typedef typename CircularKernel::Kernel_base::Line_arc_2 RLine_arc_2;
public:
  typedef  RLine_arc_2 Rep;
  typedef  CircularKernel   R;

 const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

   Line_arc_2()
     : RLine_arc_2(typename R::Construct_line_arc_2()())
   {}

   // Not Documented
   Line_arc_2(const Line_2 &support,
              const Circle_2 &c1,const bool b1,
              const Circle_2 &c2,const bool b2)
     : RLine_arc_2(typename R::Construct_line_arc_2()(support, c1, b1, c2, b2))
   {}

   // Not Documented
   Line_arc_2(const Line_2 &support,
               const Line_2 &l1,
               const Line_2 &l2)
     : RLine_arc_2(typename R::Construct_line_arc_2()(support, l1, l2))
   {}

   Line_arc_2(const Line_2 &support,
               const Circular_arc_point_2 &p1,
               const Circular_arc_point_2 &p2)
     : RLine_arc_2(typename R::Construct_line_arc_2()(support, p1, p2))
   {}

   Line_arc_2(const Segment_2 &s)
     : RLine_arc_2(typename R::Construct_line_arc_2()(s))
   {}

    Line_arc_2(const Point_2 &p1,
               const Point_2 &p2)
      : RLine_arc_2(typename R::Construct_line_arc_2()(p1, p2))
   {}

    Line_arc_2(const RLine_arc_2 &a )
     : RLine_arc_2(a)
   {}

  typename cpp11::result_of< typename R::Construct_circular_source_vertex_2(Line_arc_2)>::type
    source() const
  {
        return typename R::Construct_circular_source_vertex_2()(*this);
  }

  typename cpp11::result_of< typename R::Construct_circular_target_vertex_2(Line_arc_2)>::type
    target() const
  {
        return typename R::Construct_circular_target_vertex_2()(*this);
  }

  typename cpp11::result_of< typename R::Construct_circular_min_vertex_2(Line_arc_2)>::type
  left() const
  {
        return typename R::Construct_circular_min_vertex_2()(*this);
  }

  typename cpp11::result_of< typename R::Construct_circular_max_vertex_2(Line_arc_2)>::type
  right() const
  {
        return typename R::Construct_circular_max_vertex_2()(*this);
  }

  Line_2
  supporting_line() const
  {
        return typename R::Construct_line_2()(*this);
  }


  bool is_vertical() const
  {
      return typename R::Is_vertical_2()(*this);
  }

  Bbox_2  bbox() const
  {
        return typename R::Construct_bbox_2()(*this);
  }

 };

template < typename CircularKernel >
inline
bool
operator==(const Line_arc_2<CircularKernel> &p,
           const Line_arc_2<CircularKernel> &q)
{
  return CircularKernel().equal_2_object()(p, q);
}

template < typename CircularKernel >
inline
bool
operator!=(const Line_arc_2<CircularKernel> &p,
           const Line_arc_2<CircularKernel> &q)
{
  return ! (p == q);
}


 template < typename CK >
    std::ostream &
    operator<<(std::ostream & os, const Line_arc_2<CK> &a)
    {

      return os << a.supporting_line() << " "
                << a.source() << " "
                << a.target() << " ";
    }

  template < typename CK >
  std::istream &
  operator>>(std::istream & is, Line_arc_2<CK> &a)
  {
    typename CK::Line_2 l;
    typename CK::Circular_arc_point_2 p1;
    typename CK::Circular_arc_point_2 p2;
    is >> l >> p1 >> p2 ;
    if (is)
      a = Line_arc_2<CK>(l, p1, p2);
    return is;
  }


} //namespace CGAL

#endif // CGAL_LINE_ARC_2_H


