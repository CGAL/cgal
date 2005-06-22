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
// Author(s)     : Andreas Fabri
//                 Sven Schoenherr

#ifndef CGAL_CIRCLE_2_H
#define CGAL_CIRCLE_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Circle_2 : public R_::Kernel_base::Circle_2
{
  typedef typename R_::FT                    FT;
  typedef typename R_::Point_2               Point_2;
  typedef typename R_::Kernel_base::Circle_2  RCircle_2;

public:
  typedef RCircle_2 Rep;
  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef  R_   R;

  Circle_2()
    : RCircle_2(typename R::Construct_circle_2()().rep())
  {}

  Circle_2(const RCircle_2& t)
    : RCircle_2(t) {}

  Circle_2(const Point_2 &center, const FT &squared_radius,
	   const Orientation &orientation)
    : RCircle_2(typename R::Construct_circle_2()(center, squared_radius, orientation).rep()) {}

  Circle_2(const Point_2 &center, const FT &squared_radius)
    : RCircle_2(typename R::Construct_circle_2()(center, squared_radius, COUNTERCLOCKWISE).rep()) {}

  Circle_2(const Point_2 &p, const Point_2 &q, const Point_2 &r)
    : RCircle_2(typename R::Construct_circle_2()(p,q,r).rep()) {}

  Circle_2(const Point_2 & p, const Point_2 & q,
	   const Orientation &orientation)
    : RCircle_2(typename R::Construct_circle_2()(p,q,orientation).rep()) {}

  Circle_2(const Point_2 & p, const Point_2 & q)
    : RCircle_2(typename R::Construct_circle_2()(p,q,COUNTERCLOCKWISE).rep()) {}

  Circle_2(const Point_2 & center, const Orientation& orientation)
    : RCircle_2(typename R::Construct_circle_2()(center,FT(0),orientation).rep()) {}

  Circle_2(const Point_2 & center)
    : RCircle_2(typename R::Construct_circle_2()(center,FT(0),COUNTERCLOCKWISE).rep()) {}

  
  const Point_2 & center() const
  {
    return R().construct_center_2_object()(*this);
  }

  FT  squared_radius() const
  {
    return R().compute_squared_radius_2_object()(*this);
  }

  Orientation orientation() const
  {
    return R().orientation_2_object()(*this);
  }


  Bounded_side
  bounded_side(const Point_2 &p) const
  {
    return R().bounded_side_2_object()(*this, p);
  }
  
  Oriented_side
  oriented_side(const Point_2 &p) const
  {
    return R().oriented_side_2_object()(*this, p);
  } 

  bool
  has_on_boundary(const Point_2 &p) const
  {
    return bounded_side(p) == ON_BOUNDARY;
  }

  bool
  has_on_bounded_side(const Point_2 &p) const
  {
    return bounded_side(p) == ON_BOUNDED_SIDE;
  }

  bool
  has_on_unbounded_side(const Point_2 &p) const
  {
    return bounded_side(p) == ON_UNBOUNDED_SIDE;
  }

  bool
  has_on_negative_side(const Point_2 &p) const
  {
    if (orientation() == COUNTERCLOCKWISE)
      return has_on_unbounded_side(p);
    return has_on_bounded_side(p);
  }

  bool
  has_on_positive_side(const Point_2 &p) const
  {
    if (orientation() == COUNTERCLOCKWISE)
      return has_on_bounded_side(p);
    return has_on_unbounded_side(p);
  }

  bool
  is_degenerate() const
  {
    return CGAL_NTS is_zero(squared_radius());
  }

  Circle_2
  opposite() const
  {
    return Circle_2(center(),
		    squared_radius(),
		    CGAL::opposite(orientation()) );
  }

  Bbox_2
  bbox() const
  {
    return R().construct_bbox_2_object()(*this);
  }

  bool
  operator==(const Circle_2 &c) const
  { 
    return R().equal_2_object()(*this, c);
  }

  bool
  operator!=(const Circle_2 &c) const
  {
    return !(*this == c);
  }
};

#ifndef CGAL_NO_OSTREAM_INSERT_CIRCLE_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Circle_2<R> &c)
{
  return os << c.rep();
}

#endif // CGAL_NO_OSTREAM_INSERT_CIRCLE_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_CIRCLE_2
template < class R >
std::istream &
operator>>(std::istream &is, Circle_2<R> &c)
{
  return is >> c.rep();
}
#endif // CGAL_NO_ISTREAM_EXTRACT_CIRCLE_2

CGAL_END_NAMESPACE

#endif  // CGAL_CIRCLE_2_H
