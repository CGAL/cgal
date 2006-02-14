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
// Author(s)     : Stefan Schirra
 
#ifndef CGAL_DIRECTION_2_H
#define CGAL_DIRECTION_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Direction_2 : public R_::Kernel_base::Direction_2
{

  typedef typename R_::RT                    RT;
  typedef typename R_::Vector_2              Vector_2;
  typedef typename R_::Line_2                Line_2;
  typedef typename R_::Ray_2                 Ray_2;
  typedef typename R_::Segment_2             Segment_2;
  typedef typename R_::Kernel_base::Direction_2      RDirection_2;

public:

  typedef RDirection_2 Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef  R_   R;

  Direction_2() {}

  Direction_2(const RDirection_2& d)
    : RDirection_2(d) {}

  Direction_2(const Vector_2& v)
    : RDirection_2(typename R::Construct_direction_2()(v).rep()) {}

  Direction_2(const Line_2& l)
    : RDirection_2(typename R::Construct_direction_2()(l).rep()) {}

  Direction_2(const Ray_2& r)
    : RDirection_2(typename R::Construct_direction_2()(r).rep()) {}

  Direction_2(const Segment_2& s)
    : RDirection_2(typename R::Construct_direction_2()(s).rep()) {}

  Direction_2(const RT &x, const RT &y)
    :  RDirection_2(typename R::Construct_direction_2()(x,y).rep()) {}

  bool
  counterclockwise_in_between(const Direction_2 &d1,
			      const Direction_2 &d2) const 
  {
    return R().counterclockwise_in_between_2_object()(*this, d1, d2);
  }
  
  Direction_2 perpendicular(const Orientation &o) const
  {
    return R().construct_perpendicular_direction_2_object()(*this,o);
  }

  typename Qualified_result_of<typename R::Compute_dx_2, Direction_2>::type
  dx() const
  {
    return R().compute_dx_2_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_dy_2, Direction_2>::type
  dy() const
  {
    return R().compute_dy_2_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_dx_2, Direction_2>::type
  delta(int i) const
  {
    CGAL_kernel_precondition( ( i == 0 ) || ( i == 1 ) );
    return (i==0) ? dx() : dy();
  }
  
  bool
  operator<(const Direction_2 &d) const
  {
    return R().compare_angle_with_x_axis_2_object()(*this, d) == SMALLER;
  }


  bool
  operator>(const Direction_2 &d) const
  {
    return d < *this;
  }


  bool
  operator>=(const Direction_2 &d) const
  {
    return R().compare_angle_with_x_axis_2_object()(*this, d) != SMALLER;
  }
  

  bool
  operator<=(const Direction_2 &d) const
  {
    return R().compare_angle_with_x_axis_2_object()(*this, d) != LARGER;
  }
  
  Direction_2
  operator-() const
  {
    return Direction_2(-dx(), -dy());
  } 
  
  Vector_2 vector() const
  {
  return R().construct_vector_2_object()(*this);
  }
  
  bool
  operator==(const Direction_2& d) const
  {
    return R().equal_2_object()(*this, d);
  }


  bool
  operator!=(const Direction_2& d) const
  {
    return !(*this == d);
  }

};

#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTION_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Direction_2<R> &d)
{

  return os << d.rep();
}
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTION_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTION_2
template < class R >
std::istream &
operator>>(std::istream &is, Direction_2<R> &d)
{
  return is >> d.rep();

}
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTION_2

CGAL_END_NAMESPACE

#endif // CGAL_DIRECTION_2_H
