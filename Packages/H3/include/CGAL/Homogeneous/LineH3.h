// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : include/CGAL/Homogeneous/LineH3.h
// package       : H3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 
#ifndef CGAL_LINEH3_H
#define CGAL_LINEH3_H

#include <utility>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class LineH3
  : public R_::template Handle<std::pair<typename R_::Point_3,
                                         typename R_::Direction_3> >::type
{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::RT                   RT;
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Segment_3            Segment_3;
  typedef typename R_::Ray_3                Ray_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  typedef std::pair<Point_3, Direction_3>          rep;
  typedef typename R_::template Handle<rep>::type  base;

public:
  typedef R_                R;

  LineH3()
    : base(rep()) {}

  LineH3(const Point_3& p, const Point_3& q)
    : base(rep(p, (q - p).direction())) {}

  LineH3(const Segment_3& s)
    : base(rep(s.start(), s.direction())) {}

  LineH3(const Ray_3& r)
    : base(rep(r.start(), r.direction())) {}

  LineH3(const Point_3& p, const Direction_3& d)
    : base(rep(p, d)) {}

  Plane_3  perpendicular_plane(const Point_3& p) const;
  LineH3<R>   opposite() const;
  const Point_3 & point() const;
  Point_3  point(int i) const;
  Point_3  projection(const Point_3& p) const;

  const Direction_3 & direction() const;

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

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
inline
const typename LineH3<R>::Point_3 &
LineH3<R>::point() const
{ return Ptr()->first; }

template < class R >
CGAL_KERNEL_INLINE
typename LineH3<R>::Point_3
LineH3<R>::point(int i) const
{ return point() + RT(i)*direction().to_vector() ; }

template < class R >
inline
const typename LineH3<R>::Direction_3 &
LineH3<R>::direction() const
{ return Ptr()->second; }

template < class R >
CGAL_KERNEL_INLINE
typename LineH3<R>::Plane_3
LineH3<R>::perpendicular_plane(const typename LineH3<R>::Point_3& p ) const
{ return Plane_3( p, direction() ); }

template < class R >
CGAL_KERNEL_INLINE
LineH3<R>
LineH3<R>::opposite() const
{ return LineH3<R>( point(), -direction() ); }

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
  Vector_3 dir = direction().to_vector();
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
{ return LineH3<R>(t.transform(point() ), t.transform(direction() )); }


#ifndef CGAL_NO_OSTREAM_INSERT_LINEH3
template < class R >
std::ostream &operator<<(std::ostream &os, const LineH3<R> &l)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << l.point() << ' ' << l.direction();
    case IO::BINARY :
        return os << l.point() <<  l.direction();
    default:
        return  os << "LineH3(" << l.point() << ", " << l.direction() << ")";
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_LINEH3

#ifndef CGAL_NO_ISTREAM_EXTRACT_LINEH3
template < class R >
std::istream &operator>>(std::istream &is, LineH3<R> &l)
{
  typename R::Point_3 p;
  typename R::Direction_3 d;
  is >> p >> d;
  l = LineH3<R>(p, d);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_LINEH3

template < class R >
CGAL_KERNEL_INLINE
bool
LineH3<R>::has_on( const typename LineH3<R>::Point_3& p ) const
{ return collinear(point(), point()+direction().to_vector(), p); }

template < class R >
CGAL_KERNEL_INLINE
bool
LineH3<R>::is_degenerate() const
{ return direction().is_degenerate(); }

template < class R >
CGAL_KERNEL_INLINE
bool
LineH3<R>::operator==(const LineH3<R>& l) const
{
  return l.direction() == direction() && l.has_on( point() );
}

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_LINEH3_H
