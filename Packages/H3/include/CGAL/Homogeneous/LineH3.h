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
 
#ifndef CGAL_LINEH3_H
#define CGAL_LINEH3_H

#include <utility>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class LineH3
  : public R_::template Handle<std::pair<typename R_::Point_3,
                                         typename R_::Vector_3> >::type
{
  typedef typename R_::RT                   RT;
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Segment_3            Segment_3;
  typedef typename R_::Ray_3                Ray_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  typedef std::pair<Point_3, Vector_3>          rep;
  typedef typename R_::template Handle<rep>::type  base;

  const base& Base() const { return *this; }
  base& Base() { return *this; }

public:
  typedef R_                R;

  LineH3() {}

  LineH3(const Point_3& p, const Point_3& q)
    : base(p, q - p) {}

  LineH3(const Segment_3& s)
    : base(s.start(), s.to_vector()) {}

  LineH3(const Ray_3& r)
    : base(r.start(), r.to_vector()) {}

  LineH3(const Point_3& p, const Direction_3& d)
    : base(p, d.to_vector()) {}

  LineH3(const Point_3& p, const Vector_3& v)
    : base(p, v) {}

  Plane_3  perpendicular_plane(const Point_3& p) const;
  LineH3<R>   opposite() const;
  const Point_3 & point() const;
  Point_3  point(int i) const;
  Point_3  projection(const Point_3& p) const;

  Direction_3 direction() const;
  const Vector_3 & to_vector() const;

  bool            has_on( const Point_3& p ) const;
  bool            is_degenerate() const;

  bool            operator==(const LineH3<R>& l) const ;
  bool            operator!=(const LineH3<R>& l) const ;

  LineH3<R>   transform(const Aff_transformation_3&) const;
};

template < class R >
inline
bool
LineH3<R>::operator!=(const LineH3<R>& l) const
{ return !(*this == l); }

template < class R >
inline
const typename LineH3<R>::Point_3 &
LineH3<R>::point() const
{ return get(Base()).first; }

template < class R >
CGAL_KERNEL_INLINE
typename LineH3<R>::Point_3
LineH3<R>::point(int i) const
{ return point() + to_vector()*RT(i) ; }

template < class R >
inline
const typename LineH3<R>::Vector_3 &
LineH3<R>::to_vector() const
{ return get(Base()).second; }

template < class R >
inline
typename LineH3<R>::Direction_3
LineH3<R>::direction() const
{ return get(Base()).second; }

template < class R >
CGAL_KERNEL_INLINE
typename LineH3<R>::Plane_3
LineH3<R>::perpendicular_plane(const typename LineH3<R>::Point_3& p ) const
{ return Plane_3( p, to_vector() ); }

template < class R >
CGAL_KERNEL_INLINE
LineH3<R>
LineH3<R>::opposite() const
{ return LineH3<R>( point(), -to_vector() ); }

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename LineH3<R>::Point_3
LineH3<R>::projection(const typename LineH3<R>::Point_3& p) const
{
  if ( has_on(p) )
  {
      return p;
  }
  Vector_3  v = p - point();
  const RT  vx = v.hx();
  const RT  vy = v.hy();
  const RT  vz = v.hz();
  const RT  vw = v.hw();
  Vector_3 dir = to_vector();
  const RT  dx = dir.hx();
  const RT  dy = dir.hy();
  const RT  dz = dir.hz();
  const RT  dw = dir.hw();

  RT lambda_num = (vx*dx + vy*dy + vz*dz)*dw; // *dw
  RT lambda_den = (dx*dx + dy*dy + dz*dz)*vw; // *dw

  return point() + ( (lambda_num * dir)/lambda_den );
}

template < class R >
CGAL_KERNEL_INLINE
LineH3<R>
LineH3<R>::transform(const typename LineH3<R>::Aff_transformation_3& t) const
{ return LineH3<R>(t.transform(point() ), t.transform(to_vector() )); }


#ifndef CGAL_NO_OSTREAM_INSERT_LINEH3
template < class R >
std::ostream &operator<<(std::ostream &os, const LineH3<R> &l)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << l.point() << ' ' << l.to_vector();
    case IO::BINARY :
        return os << l.point() <<  l.to_vector();
    default:
        return  os << "LineH3(" << l.point() << ", " << l.to_vector() << ")";
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_LINEH3

#ifndef CGAL_NO_ISTREAM_EXTRACT_LINEH3
template < class R >
std::istream &operator>>(std::istream &is, LineH3<R> &l)
{
  typename R::Point_3 p;
  typename R::Vector_3 v;
  is >> p >> v;
  l = LineH3<R>(p, v);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_LINEH3

template < class R >
CGAL_KERNEL_INLINE
bool
LineH3<R>::has_on( const typename LineH3<R>::Point_3& p ) const
{ return collinear(point(), point()+to_vector(), p); }

template < class R >
CGAL_KERNEL_INLINE
bool
LineH3<R>::is_degenerate() const
{ return to_vector() == NULL_VECTOR; }

template < class R >
CGAL_KERNEL_INLINE
bool
LineH3<R>::operator==(const LineH3<R>& l) const
{
  return l.to_vector() == to_vector() && l.has_on( point() );
}

CGAL_END_NAMESPACE

#endif // CGAL_LINEH3_H
