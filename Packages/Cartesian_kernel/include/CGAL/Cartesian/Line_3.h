// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_LINE_3_H
#define CGAL_CARTESIAN_LINE_3_H

#include <utility>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class LineC3
  : public R_::template Handle<std::pair<typename R_::Point_3,
                                         typename R_::Vector_3> >::type
{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Ray_3                Ray_3;
  typedef typename R_::Line_3               Line_3;
  typedef typename R_::Segment_3            Segment_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  typedef std::pair<Point_3, Vector_3>             rep;
  typedef typename R_::template Handle<rep>::type  base;

public:
  typedef R_                                     R;

  LineC3() {}

  LineC3(const Point_3 &p, const Point_3 &q)
    : base(rep(p, q-p)) {}

  LineC3(const Segment_3 &s)
    : base(R().construct_line_3_object()(s)) {}

  LineC3(const Ray_3 &r)
    : base(R().construct_line_3_object()(r)) {}

  LineC3(const Point_3 &p, const Vector_3 &v)
    : base(rep(p, v)) {}

  LineC3(const Point_3 &p, const Direction_3 &d)
    : base(rep(p, Vector_3(d.dx(), d.dy(), d.dz()))) {}

  bool        operator==(const LineC3 &l) const;
  bool        operator!=(const LineC3 &l) const;

  Plane_3     perpendicular_plane(const Point_3 &p) const;
  Line_3      opposite() const;

  const Point_3 &     point() const
  {
      return Ptr()->first;
  }

  const Vector_3 & to_vector() const
  {
      return Ptr()->second;
  }

  Direction_3 direction() const
  {
      return Direction_3(Ptr()->second);
  }

  Point_3     point(int i) const;

  Point_3     projection(const Point_3 &p) const;

  bool        has_on(const Point_3 &p) const;
  bool        is_degenerate() const;

  Line_3        transform(const Aff_transformation_3 &t) const
  {
    return LineC3<R>(t.transform(point()), t.transform(direction()));
  }
};

template < class R >
inline
bool
LineC3<R>::operator==(const LineC3<R> &l) const
{
  if (identical(l))
      return true;
  return has_on(l.point()) && (direction() == l.direction());
}

template < class R >
inline
bool
LineC3<R>::operator!=(const LineC3<R> &l) const
{
  return !(*this == l);
}

template < class R >
inline
typename LineC3<R>::Point_3
LineC3<R>::point(int i) const
{
  return point_on_line(i, *this);
}

template < class R >
inline
typename LineC3<R>::Plane_3
LineC3<R>::
perpendicular_plane(const typename LineC3<R>::Point_3 &p) const
{
  return Plane_3(p, to_vector());
}

template < class R >
inline
typename LineC3<R>::Line_3
LineC3<R>::opposite() const
{
  return LineC3<R>(point(), -to_vector());
}

template < class R >
inline
typename LineC3<R>::Point_3
LineC3<R>::
projection(const typename LineC3<R>::Point_3 &p) const
{
  return projection_line(p, *this);
}

template < class R >
inline
bool
LineC3<R>::
has_on(const typename LineC3<R>::Point_3 &p) const
{
  return collinear(point(), point()+to_vector(), p);
}

template < class R >
inline
bool
LineC3<R>::is_degenerate() const
{ // FIXME : predicate
  return to_vector() == NULL_VECTOR;
}

#ifndef CGAL_CARTESIAN_NO_OSTREAM_INSERT_LINEC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const LineC3<R> &l)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << l.point(0) << ' ' << l.point(1);
    case IO::BINARY :
        return os << l.point(0) <<  l.point(1);
    default:
        return  os << "LineC3(" << l.point(0) << ", " << l.point(1) << ")";
    }
}
#endif // CGAL_CARTESIAN_NO_OSTREAM_INSERT_LINEC3

#ifndef CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_LINEC3
template < class R >
std::istream &
operator>>(std::istream &is, LineC3<R> &l)
{
    typename R::Point_3 p, q;
    is >> p >> q;
    if (is)
	l = LineC3<R>(p, q);
    return is;
}
#endif // CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_LINEC3

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_LINE_3_H
