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
// file          : LineH3.h
// package       : H3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_LINEH3_H
#define CGAL_LINEH3_H

#include <CGAL/PVDH3.h>
#include <CGAL/SegmentH3.h>
#include <CGAL/RayH3.h>
#include <CGAL/predicates_on_pointsH3.h>

CGAL_BEGIN_NAMESPACE

template < class R >
class Line_repH3 : public Ref_counted
{
  public:
     Line_repH3() {}
     Line_repH3( const PointH3<R>& p, const DirectionH3<R> d)
      : basepoint(p), direction(d)
     {}

  friend class LineH3<R>;

  private:
    PointH3<R>       basepoint;
    DirectionH3<R>   direction;
};


template < class R_ >
class LineH3
  : public R_::Line_handle_3
{
public:
  typedef R_                R;
  typedef typename R::RT    RT;
  typedef typename R::FT    FT;

  LineH3();
  LineH3(const PointH3<R>& p, const PointH3<R>& q);
  LineH3(const SegmentH3<R>& s);
  LineH3(const RayH3<R>& r);
  LineH3(const PointH3<R>& p, const DirectionH3<R>& d);

  PlaneH3<R>  perpendicular_plane(const PointH3<R>& p) const;
  LineH3<R>   opposite() const;
  PointH3<R>  point() const;
  PointH3<R>  point(int i) const;
  PointH3<R>  projection(const PointH3<R>& p) const;

  DirectionH3<R>
                  direction() const;

  bool            has_on( const PointH3<R>& p ) const;
  bool            is_degenerate() const;

  bool            operator==(const LineH3<R>& l) const ;
  bool            operator!=(const LineH3<R>& l) const ;

  LineH3<R>   transform(const Aff_transformationH3<R>&) const;
};

template < class R >
inline
bool
LineH3<R>::operator!=(const LineH3<R>& l) const
{ return !(*this == l); }


CGAL_END_NAMESPACE

#include <CGAL/PlaneH3.h>

CGAL_BEGIN_NAMESPACE

template < class R >
CGAL_KERNEL_CTOR_INLINE
LineH3<R>::LineH3()
 : Handle_for< Line_repH3<R> >( Line_repH3<R>() )
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
LineH3<R>::LineH3(const PointH3<R>& p,
                      const PointH3<R>& q)
 : Handle_for< Line_repH3<R> >(Line_repH3<R>(p, (q - p).direction()))
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
LineH3<R>::LineH3(const SegmentH3<R>& s)
 : Handle_for< Line_repH3<R> >(Line_repH3<R>(s.start(), s.direction()))
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
LineH3<R>::LineH3(const RayH3<R>& r)
 : Handle_for< Line_repH3<R> >(Line_repH3<R>(r.start(), r.direction()))
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
LineH3<R>::LineH3(const PointH3<R>& p,
                      const DirectionH3<R>& d)
 : Handle_for< Line_repH3<R> >(Line_repH3<R>( p, d ) )
{}

template < class R >
inline
PointH3<R>
LineH3<R>::point() const
{ return ptr->basepoint; }

template < class R >
CGAL_KERNEL_INLINE
PointH3<R>
LineH3<R>::point(int i) const
{ return point() + RT(i)*direction().to_vector() ; }

template < class R >
inline
DirectionH3<R>
LineH3<R>::direction() const
{ return ptr->direction; }

template < class R >
CGAL_KERNEL_INLINE
PlaneH3<R>
LineH3<R>::
perpendicular_plane(const PointH3<R>& p ) const
{ return PlaneH3<R>( p, direction() ); }

template < class R >
CGAL_KERNEL_INLINE
LineH3<R>
LineH3<R>::opposite() const
{ return LineH3<R>( ptr->basepoint, -(ptr->direction ) ); }

template < class R >
CGAL_KERNEL_LARGE_INLINE
PointH3<R>
LineH3<R>::projection(const PointH3<R>& p) const
{
  if ( has_on(p) )
  {
      return p;
  }
  VectorH3<R>  v = p - point();
  const RT  vx = v.hx();
  const RT  vy = v.hy();
  const RT  vz = v.hz();
  const RT  vw = v.hw();
  VectorH3<R> dir = direction().to_vector();
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
LineH3<R>::transform(const Aff_transformationH3<R>& t) const
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
  PointH3<R> p;
  DirectionH3<R> d;
  is >> p >> d;
  l = LineH3<R>(p, d);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_LINEH3

template < class R >
CGAL_KERNEL_INLINE
bool
LineH3<R>::has_on( const PointH3<R>& p ) const
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
  return  (  (l.direction() ==   ptr->direction )
           &&(l.has_on( ptr->basepoint ) ) );
}

CGAL_END_NAMESPACE

#endif // CGAL_LINEH3_H
