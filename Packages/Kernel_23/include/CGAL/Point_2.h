// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri, Stefan Schirra

#ifndef CGAL_POINT_2_H
#define CGAL_POINT_2_H

#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>

CGAL_BEGIN_NAMESPACE


template <class R_>
class Point_2 : public R_::Kernel_base::Point_2
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Vector_2              Vector_2;
  typedef typename R_::Kernel_base::Point_2  RPoint_2;

public:
  typedef RPoint_2 Rep;
  typedef typename R_::Cartesian_const_iterator_2 Cartesian_const_iterator;
  typedef typename R_::Cartesian_coordinate_type Cartesian_coordinate_type;
  typedef typename R_::Homogeneous_coordinate_type Homogeneous_coordinate_type;


  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef  R_   R;

  Point_2() 
    : RPoint_2(typename R::Construct_point_2()().rep())
  {}

  Point_2(const Origin& o)
    : RPoint_2(typename R::Construct_point_2()(o).rep())
  {}

#if 1 // still needed by Min_ellipse_2...
  Point_2(const RPoint_2& p)
    : RPoint_2(p)
  {}
#endif

  Point_2(const RT& x, const RT& y)
    : RPoint_2(typename R::Construct_point_2()(x, y).rep())
  {}

  Point_2(const RT& hx, const RT& hy, const RT& hw)
    : RPoint_2(typename R::Construct_point_2()(hx, hy, hw).rep())
  {}

  typename Qualified_result_of<typename R::Compute_x_2,Point_2>::type
  x() const
  {
    return typename R::Compute_x_2()(*this);
  }

  Cartesian_coordinate_type y() const
  {
    return typename R::Compute_y_2()(*this);
  }

   Cartesian_coordinate_type cartesian(int i) const
  {
    CGAL_kernel_precondition( (i == 0) || (i == 1) );
    return (i==0) ?  x() : y();
  }

  Cartesian_coordinate_type 
  operator[](int i) const
  {
      return cartesian(i);
  }

  Cartesian_const_iterator cartesian_begin() const 
  {
    return typename R::Construct_cartesian_const_iterator_2()(*this);
  }

  Cartesian_const_iterator cartesian_end() const 
  {
    return typename R::Construct_cartesian_const_iterator_2()(*this,1);
  }



  Homogeneous_coordinate_type hx() const
  {
    return typename R::Compute_hx_2()(*this);
  }

  Homogeneous_coordinate_type hy() const
  {
    return typename R::Compute_hy_2()(*this);
  }

  Homogeneous_coordinate_type hw() const
  {
    return typename R::Compute_hw_2()(*this);
  }

  int dimension() const
  {
      return 2;
  }

   Cartesian_coordinate_type homogeneous(int i) const
  {
    CGAL_kernel_precondition( (i >= 0) || (i <= 2) );
    return (i==0) ?  hx() : (i==1)? hy() : hw();
  }

  Bbox_2 bbox() const
  {
    return R().construct_bbox_2_object()(*this);
  }
  
  bool
  operator==(const Point_2& p) const
  {
    return R().equal_2_object()(*this, p);
  }


  bool
  operator!=(const Point_2& p) const
  {
    return !(*this == p);
  }


};



#ifndef CGAL_NO_OSTREAM_INSERT_POINT_2
template < class R >
std::ostream&
operator<<(std::ostream& os, const Point_2<R>& p)
{
  return os << p.rep();
}
#endif // CGAL_NO_OSTREAM_INSERT_POINT_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_POINT_2
template < class R >
std::istream&
operator>>(std::istream& is, Point_2<R>& p)
{
  return is >> p.rep();
}
#endif // CGAL_NO_ISTREAM_EXTRACT_POINT_2

CGAL_END_NAMESPACE

#endif // CGAL_POINT_2_H
