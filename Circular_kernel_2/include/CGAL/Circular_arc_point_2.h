// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_ARC_POINT_2_H
#define CGAL_CIRCULAR_ARC_POINT_2_H

#include <CGAL/license/Circular_kernel_2.h>


#include <CGAL/Bbox_2.h>
#include <CGAL/enum.h>

namespace CGAL {

template < typename CircularKernel >
class Circular_arc_point_2
  : public CircularKernel::Kernel_base::Circular_arc_point_2
{
  typedef typename CircularKernel::Kernel_base::Circular_arc_point_2
                                               RCircular_arc_point_2;
  typedef typename CircularKernel::Point_2     Point_2;
  typedef typename CircularKernel::Circle_2    Circle_2;

  typedef typename CircularKernel::Root_of_2   Root_of_2;

public:
  typedef typename CircularKernel::Root_for_circles_2_2
                                               Root_for_circles_2_2;
  typedef CircularKernel                       R;
  typedef RCircular_arc_point_2                Rep;


  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  Circular_arc_point_2()
    : RCircular_arc_point_2(typename R::Construct_circular_arc_point_2()())
  {}


  Circular_arc_point_2(const Root_for_circles_2_2 & np)
    : RCircular_arc_point_2(typename R::Construct_circular_arc_point_2()(np))
  {}

  Circular_arc_point_2(const RCircular_arc_point_2 & p)
    : RCircular_arc_point_2(p)
  {}

  Circular_arc_point_2(const Point_2 & p)
    : RCircular_arc_point_2(typename R::Construct_circular_arc_point_2()(p))
  {}

  decltype(auto)
  x() const
  {
    return typename R::Compute_circular_x_2()(*this);
  }

  decltype(auto)
  y() const
  {
    return typename R::Compute_circular_y_2()(*this);
  }

  Bbox_2  bbox() const
  {
    return typename R::Construct_bbox_2()(*this);
  }

};

  template < typename CircularKernel >
  inline
  bool
  operator==(const Circular_arc_point_2<CircularKernel> &p,
             const Circular_arc_point_2<CircularKernel> &q)
  {
    return CircularKernel().equal_2_object()(p, q);
  }

  template < typename CircularKernel >
  inline
  bool
  operator!=(const Circular_arc_point_2<CircularKernel> &p,
             const Circular_arc_point_2<CircularKernel> &q)
  {
    return ! (p == q);
  }

  template < typename CircularKernel >
  inline
  bool
  operator<(const Circular_arc_point_2<CircularKernel> &p,
             const Circular_arc_point_2<CircularKernel> &q)
  {
    return CircularKernel().compare_xy_2_object()(p, q) == CGAL::SMALLER;
  }

  template < typename CircularKernel >
  inline
  bool
  operator>(const Circular_arc_point_2<CircularKernel> &p,
     const Circular_arc_point_2<CircularKernel> &q)
  {
    return CircularKernel().compare_xy_2_object()(p, q) == CGAL::LARGER;
  }

  template < typename CircularKernel >
  inline
  bool
  operator<=(const Circular_arc_point_2<CircularKernel> &p,
             const Circular_arc_point_2<CircularKernel> &q)
        {
                CGAL::Comparison_result c = CircularKernel().compare_xy_2_object()(p, q);
    return (c == CGAL::SMALLER) || (c == CGAL::EQUAL);
        }

  template < typename CircularKernel >
  inline
  bool
  operator>=(const Circular_arc_point_2<CircularKernel> &p,
     const Circular_arc_point_2<CircularKernel> &q)
  {
    CGAL::Comparison_result c = CircularKernel().compare_xy_2_object()(p, q);
    return (c == CGAL::LARGER) || (c == CGAL::EQUAL);
  }

  template < typename CK >
  std::istream &
  operator>>(std::istream & is, Circular_arc_point_2<CK> &p)
  {
    typedef typename CK::Root_for_circles_2_2 Root_for_circles_2_2;

    Root_for_circles_2_2 r;
    is >> r;
    if(is)
      p = Circular_arc_point_2<CK>(r);
    return is;
  }

  template < class CK >
  std::ostream&
  operator<<(std::ostream &os, const Circular_arc_point_2<CK> &p)
  {
    return os << p.x() << " " << p.y() << " ";
  }


} //namespace CGAL

#endif // CGAL_CIRCULAR_ARC_POINT_2_H
