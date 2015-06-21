// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_LINE_2_H
#define CGAL_LINE_2_H

#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template <class R_>
class Line_2 : public R_::Kernel_base::Line_2
{
  typedef typename R_::RT                    RT;
  typedef typename R_::FT                    FT;
  typedef typename R_::Point_2               Point_2;
  typedef typename R_::Segment_2             Segment_2;
  typedef typename R_::Ray_2                 Ray_2;
  typedef typename R_::Vector_2              Vector_2;
  typedef typename R_::Direction_2           Direction_2;
  typedef typename R_::Aff_transformation_2  Aff_transformation_2;
  typedef typename R_::Kernel_base::Line_2   RLine_2;

  typedef Line_2                             Self;
  CGAL_static_assertion((boost::is_same<Self, typename R_::Line_2>::value));

public:

  typedef Dimension_tag<2>  Ambient_dimension;
  typedef Dimension_tag<1>  Feature_dimension;

  typedef RLine_2 Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef  R_   R;

  Line_2() {}

  Line_2(const RLine_2& l)  // conversion impl -> interface class
    : RLine_2(l) {}

  Line_2(const Point_2 &p, const Point_2 &q)
    : RLine_2(typename R::Construct_line_2()(Return_base_tag(), p,q)) {}

  Line_2(const RT &a, const RT &b, const RT &c)
    : RLine_2(typename R::Construct_line_2()(Return_base_tag(), a,b,c)) {}

  explicit Line_2(const Segment_2& s)
    : RLine_2(typename R::Construct_line_2()(Return_base_tag(), s)) {}

  explicit Line_2(const Ray_2& r)
    : RLine_2(typename R::Construct_line_2()(Return_base_tag(), r)) {}

  Line_2(const Point_2 &p, const Direction_2 &d)
    : RLine_2(typename R::Construct_line_2()(Return_base_tag(), p,d)) {}

  Line_2(const Point_2 &p, const Vector_2 &v)
    : RLine_2(typename R::Construct_line_2()(Return_base_tag(), p,v)) {}


  // FIXME : Use Qrt<> here.
  RT a() const
  {
    return R().compute_a_2_object()(*this);
  }

  RT b() const
  {
    return R().compute_b_2_object()(*this);
  }

  RT c() const
  {
    return R().compute_c_2_object()(*this);
  }

  Line_2
  transform(const Aff_transformation_2 &t) const
  {
    return Line_2(t.transform(point(0)),
		  t.transform(direction()));
  }

  Line_2
  opposite() const
  {
    return R().construct_opposite_line_2_object()(*this);
  }

  Direction_2
  direction() const
  {
    return R().construct_direction_2_object()(*this);
  }

  Vector_2
  to_vector() const
  {
    return R().construct_vector_2_object()(*this);
  }

  Line_2
  perpendicular(const Point_2 &p) const
  {
    return R().construct_perpendicular_line_2_object()(*this,p);
  }

  Point_2
  projection(const Point_2& p) const
  {
    return R().construct_projected_point_2_object()(*this,p);
  }

  typename R::Boolean
  is_horizontal() const
  {
    return R().is_horizontal_2_object()(*this);
  }

  typename R::Boolean
  is_vertical() const
  {
    return R().is_vertical_2_object()(*this);
  }

  typename R::Boolean
  is_degenerate() const
  { return R().is_degenerate_2_object()(*this); }

  typename R::Oriented_side
  oriented_side(const Point_2 &p) const
  {
    return R().oriented_side_2_object()(*this,p);
  }

  typename R::Boolean
  has_on_boundary(const Point_2 &p) const
  {
    return oriented_side(p) == ON_ORIENTED_BOUNDARY;
  }

  typename R::Boolean
  has_on_positive_side(const Point_2 &p) const
  {
    return oriented_side(p) == ON_POSITIVE_SIDE;
  }

  typename R::Boolean
  has_on_negative_side(const Point_2 &p) const
  {
    return oriented_side(p) == ON_NEGATIVE_SIDE;
  }

  typename R::Boolean
  has_on(const Point_2 &p) const
  {
    return has_on_boundary(p);
  }

  FT
  x_at_y(const FT &y) const
  {
    return R().compute_x_at_y_2_object()(*this, y);
  }

  FT
  y_at_x(const FT &y) const
  {
    return R().compute_y_at_x_2_object()(*this, y);
  }

  Point_2
  point() const
  {
    return R().construct_point_2_object()(*this);
  }

  Point_2
  point(int i) const
  {
    return R().construct_point_2_object()(*this,i);
  }

  typename R::Boolean
  operator==(const Line_2 &l) const
  {
    return R().equal_2_object()(*this, l);
  }

  typename R::Boolean
  operator!=(const Line_2 &l) const
  {
    return !(*this == l);
  }

};


template <class R >
std::ostream&
insert(std::ostream& os, const Line_2<R>& l)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << l.a() << ' ' << l.b() << ' ' << l.c();
    case IO::BINARY :
        write(os, l.a());
        write(os, l.b());
        write(os, l.c());
        return os;
    default:
        return os << "Line_2(" << l.a()
		  << ", " << l.b() << ", " << l.c() <<')';
    }
}

template < class R >
std::ostream &
operator<<(std::ostream &os, const Line_2<R> &l)
{
  return insert(os, l);
}


template <class R >
std::istream&
extract(std::istream& is, Line_2<R>& l)
{
    typename R::RT a, b, c;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> iformat(a) >> iformat(b) >> iformat(c);
        break;
    case IO::BINARY :
        read(is, a);
        read(is, b);
        read(is, c);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
	l = Line_2<R>(a, b, c);
    return is;
}


template < class R >
std::istream &
operator>>(std::istream &is, Line_2<R> &l)
{
  return extract(is, l);
}

} //namespace CGAL

#endif  // CGAL_LINE_2_H
