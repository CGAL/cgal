// Copyright (c) 2000  Utrecht University (The Netherlands),
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

#ifndef CGAL_CARTESIAN_PLANE_3_H
#define CGAL_CARTESIAN_PLANE_3_H

#include <CGAL/Fourtuple.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class PlaneC3
  : public R_::template Handle<Fourtuple<typename R_::FT> >::type
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Line_3               Line_3;
  typedef typename R_::Ray_3                Ray_3;
  typedef typename R_::Segment_3            Segment_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;
  typedef typename R_::Construct_point_3    Construct_point_3;
  typedef typename R_::Construct_point_2    Construct_point_2;

  typedef Fourtuple<FT>	                           rep;
  typedef typename R_::template Handle<rep>::type  base;

  const base& Base() const { return *this; }
  base& Base() { return *this; }

public:
  typedef R_                                     R;

  PlaneC3() {}

  PlaneC3(const Point_3 &p, const Point_3 &q, const Point_3 &r)
    : base(plane_from_points(p, q, r)) {}

  PlaneC3(const Point_3 &p, const Direction_3 &d)
    : base(plane_from_point_direction(p, d)) {}

  PlaneC3(const Point_3 &p, const Vector_3 &v)
    : base(plane_from_point_direction(p, v.direction())) {}

  PlaneC3(const FT &a, const FT &b, const FT &c, const FT &d)
    : base(a, b, c, d) {}

  PlaneC3(const Line_3 &l, const Point_3 &p)
    : base(plane_from_points(l.point(),
	                     l.point()+l.direction().to_vector(),
			     p)) {}

  PlaneC3(const Segment_3 &s, const Point_3 &p)
    : base(plane_from_points(s.start(), s.end(), p)) {}

  PlaneC3(const Ray_3 &r, const Point_3 &p)
    : base(plane_from_points(r.start(), r.second_point(), p)) {}

  bool         operator==(const PlaneC3 &p) const;
  bool         operator!=(const PlaneC3 &p) const;

  const FT & a() const
  {
      return get(Base()).e0;
  }
  const FT & b() const
  {
      return get(Base()).e1;
  }
  const FT & c() const
  {
      return get(Base()).e2;
  }
  const FT & d() const
  {
      return get(Base()).e3;
  }

  Line_3       perpendicular_line(const Point_3 &p) const;
  Plane_3      opposite() const;

  Point_3      point() const;
  Point_3      projection(const Point_3 &p) const;
  Vector_3     orthogonal_vector() const;
  Direction_3  orthogonal_direction() const;
  Vector_3     base1() const;
  Vector_3     base2() const;

  Point_3      to_plane_basis(const Point_3 &p) const;

  Point_2      to_2d(const Point_3 &p) const;
  Point_3      to_3d(const Point_2 &p) const;

  Plane_3      transform(const Aff_transformation_3 &t) const
  {
    if (t.is_even())
      return PlaneC3<R>(t.transform(point()),
                 t.transpose().inverse().transform(orthogonal_direction()));
    else
      return PlaneC3<R>( t.transform(point()),
               - t.transpose().inverse().transform(orthogonal_direction()));
  }

  Oriented_side oriented_side(const Point_3 &p) const;
  bool         has_on_positive_side(const Point_3 &l) const;
  bool         has_on_negative_side(const Point_3 &l) const;
  bool         has_on(const Point_3 &p) const
  {
    return oriented_side(p) == ON_ORIENTED_BOUNDARY;
  }
  bool         has_on(const Line_3 &l) const
  {
    return has_on(l.point())
       &&  has_on(l.point() + l.direction().to_vector());
  }

  bool         is_degenerate() const;
};

template < class R >
CGAL_KERNEL_INLINE
bool
PlaneC3<R>::operator==(const PlaneC3<R> &p) const
{
  if (CGAL::identical(Base(), p.Base()))
      return true;
  return equal_plane(*this, p);
}

template < class R >
inline
bool
PlaneC3<R>::operator!=(const PlaneC3<R> &p) const
{
  return !(*this == p);
}

template < class R >
inline
typename PlaneC3<R>::Point_3
PlaneC3<R>::point() const
{
  return point_on_plane(*this);
}

template < class R >
inline
typename PlaneC3<R>::Point_3
PlaneC3<R>::
projection(const typename PlaneC3<R>::Point_3 &p) const
{
  return projection_plane(p, *this);
}

template < class R >
inline
typename PlaneC3<R>::Vector_3
PlaneC3<R>::orthogonal_vector() const
{
  return Vector_3(a(), b(), c());
}

template < class R >
inline
typename PlaneC3<R>::Direction_3
PlaneC3<R>::orthogonal_direction() const
{
  return Direction_3(a(), b(), c());
}

template < class R >
typename PlaneC3<R>::Vector_3
PlaneC3<R>::base1() const
{
  return R().construct_base_vector_3_object()(*this, 1);
}

template < class R >
typename PlaneC3<R>::Vector_3
PlaneC3<R>::base2() const
{
  return R().construct_base_vector_3_object()(*this, 2);
}

template < class R >
typename PlaneC3<R>::Point_3
PlaneC3<R>::
to_plane_basis(const typename PlaneC3<R>::Point_3 &p) const
{
  FT alpha, beta, gamma;
  Construct_point_3 construct_point_3;
  solve(base1(), base2(), orthogonal_vector(), p - point(),
	alpha, beta, gamma);

  return construct_point_3(alpha, beta, gamma);
}

template < class R >
typename PlaneC3<R>::Point_2
PlaneC3<R>::
to_2d(const typename PlaneC3<R>::Point_3 &p) const
{
  FT alpha, beta, gamma;
  Construct_point_2 construct_point_2;

  solve(base1(), base2(), orthogonal_vector(), p - point(),
	alpha, beta, gamma);

  return construct_point_2(alpha, beta);
}

template < class R >
inline
typename PlaneC3<R>::Point_3
PlaneC3<R>::
to_3d(const typename PlaneC3<R>::Point_2 &p) const
{
  return R().construct_lifted_point_3_object()(*this, p);
}

template < class R >
inline
typename PlaneC3<R>::Line_3
PlaneC3<R>::
perpendicular_line(const typename PlaneC3<R>::Point_3 &p) const
{
  return Line_3(p, orthogonal_direction());
}

template < class R >
inline
typename PlaneC3<R>::Plane_3
PlaneC3<R>::opposite() const
{
  return PlaneC3<R>(-a(), -b(), -c(), -d());
}

template < class R >
inline
Oriented_side
PlaneC3<R>::
oriented_side(const typename PlaneC3<R>::Point_3 &p) const
{
  return side_of_oriented_plane(*this, p);
}

template < class R >
inline
bool
PlaneC3<R>::
has_on_positive_side(const  typename PlaneC3<R>::Point_3 &p) const
{
  return oriented_side(p) == ON_POSITIVE_SIDE;
}

template < class R >
inline
bool
PlaneC3<R>::
has_on_negative_side(const  typename PlaneC3<R>::Point_3 &p) const
{
  return oriented_side(p) == ON_NEGATIVE_SIDE;
}

template < class R >
inline
bool
PlaneC3<R>::
is_degenerate() const
{ // FIXME : predicate
  return CGAL_NTS is_zero(a()) && CGAL_NTS is_zero(b()) &&
         CGAL_NTS is_zero(c());
}

#ifndef CGAL_NO_OSTREAM_INSERT_PLANEC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const PlaneC3<R> &p)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << p.a() << ' ' << p.b() <<  ' ' << p.c() << ' ' << p.d();
    case IO::BINARY :
        write(os, p.a());
        write(os, p.b());
        write(os, p.c());
        write(os, p.d());
        return os;
        default:
            os << "PlaneC3(" << p.a() <<  ", " << p.b() <<   ", ";
            os << p.c() << ", " << p.d() <<")";
            return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_PLANEC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_PLANEC3
template < class R >
std::istream &
operator>>(std::istream &is, PlaneC3<R> &p)
{
    typename R::FT a, b, c, d;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> a >> b >> c >> d;
        break;
    case IO::BINARY :
        read(is, a);
        read(is, b);
        read(is, c);
        read(is, d);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
	p = PlaneC3<R>(a, b, c, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_PLANEC3

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_PLANE_3_H
