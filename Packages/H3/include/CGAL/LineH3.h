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

#ifndef CGAL_PVDH3_H
#include <CGAL/PVDH3.h>
#endif // CGAL_PVDH3_H
#ifndef CGAL_SEGMENTH3_H
#include <CGAL/SegmentH3.h>
#endif // CGAL_SEGMENTH3_H
#ifndef CGAL_RAYH3_H
#include <CGAL/RayH3.h>
#endif // CGAL_RAYH3_H
#ifndef PREDICATES_ON_POINTSH3_H
#include <CGAL/predicates_on_pointsH3.h>
#endif // PREDICATES_ON_POINTSH3_H

CGAL_BEGIN_NAMESPACE

template < class FT, class RT >
class Line_repH3 : public Ref_counted
{
  public:
     Line_repH3() {}
     Line_repH3( const PointH3<FT,RT>& p, const DirectionH3<FT,RT> d)
      : basepoint(p), direction(d)
     {}

  friend class LineH3<FT,RT>;

  private:
    PointH3<FT,RT>       basepoint;
    DirectionH3<FT,RT>   direction;
};


template < class FT, class RT >
class LineH3 : public Handle_for< Line_repH3<FT,RT> >
{
public:
  LineH3();
  LineH3(const PointH3<FT,RT>& p, const PointH3<FT,RT>& q);
  LineH3(const SegmentH3<FT,RT>& s);
  LineH3(const RayH3<FT,RT>& r);
  LineH3(const PointH3<FT,RT>& p, const DirectionH3<FT,RT>& d);

  PlaneH3<FT,RT>  perpendicular_plane(const PointH3<FT,RT>& p) const;
  LineH3<FT,RT>   opposite() const;
  PointH3<FT,RT>  point() const;
  PointH3<FT,RT>  point(int i) const;
  PointH3<FT,RT>  projection(const PointH3<FT,RT>& p) const;

  DirectionH3<FT,RT>
                  direction() const;

  bool            has_on( const PointH3<FT,RT>& p ) const;
  bool            is_degenerate() const;

  bool            operator==(const LineH3<FT,RT>& l) const ;
  bool            operator!=(const LineH3<FT,RT>& l) const ;

  LineH3<FT,RT>   transform(const Aff_transformationH3<FT,RT>&) const;
};

template < class FT, class RT >
inline
bool
LineH3<FT,RT>::operator!=(const LineH3<FT,RT>& l) const
{ return !(*this == l); }


CGAL_END_NAMESPACE


#ifndef CGAL_PLANEH3_H
#include <CGAL/PlaneH3.h>
#endif // CGAL_PLANEH3_H

CGAL_BEGIN_NAMESPACE

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
LineH3<FT,RT>::LineH3()
 : Handle_for< Line_repH3<FT,RT> >( Line_repH3<FT,RT>() )
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
LineH3<FT,RT>::LineH3(const PointH3<FT,RT>& p,
                      const PointH3<FT,RT>& q)
 : Handle_for< Line_repH3<FT,RT> >(Line_repH3<FT,RT>(p, (q - p).direction()))
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
LineH3<FT,RT>::LineH3(const SegmentH3<FT,RT>& s)
 : Handle_for< Line_repH3<FT,RT> >(Line_repH3<FT,RT>(s.start(), s.direction()))
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
LineH3<FT,RT>::LineH3(const RayH3<FT,RT>& r)
 : Handle_for< Line_repH3<FT,RT> >(Line_repH3<FT,RT>(r.start(), r.direction()))
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
LineH3<FT,RT>::LineH3(const PointH3<FT,RT>& p,
                      const DirectionH3<FT,RT>& d)
 : Handle_for< Line_repH3<FT,RT> >(Line_repH3<FT,RT>( p, d ) )
{}

template < class FT, class RT >
inline
PointH3<FT,RT>
LineH3<FT,RT>::point() const
{ return ptr->basepoint; }

template < class FT, class RT >
CGAL_KERNEL_INLINE
PointH3<FT,RT>
LineH3<FT,RT>::point(int i) const
{ return point() + RT(i)*direction().to_vector() ; }

template < class FT, class RT >
inline
DirectionH3<FT,RT>
LineH3<FT,RT>::direction() const
{ return ptr->direction; }

template < class FT, class RT >
CGAL_KERNEL_INLINE
PlaneH3<FT,RT>
LineH3<FT,RT>::
perpendicular_plane(const PointH3<FT,RT>& p ) const
{ return PlaneH3<FT,RT>( p, direction() ); }

template < class FT, class RT >
CGAL_KERNEL_INLINE
LineH3<FT,RT>
LineH3<FT,RT>::opposite() const
{ return LineH3<FT,RT>( ptr->basepoint, -(ptr->direction ) ); }

template < class FT, class RT >
CGAL_KERNEL_LARGE_INLINE
PointH3<FT,RT>
LineH3<FT,RT>::projection(const PointH3<FT,RT>& p) const
{
  if ( has_on(p) )
  {
      return p;
  }
  VectorH3<FT,RT>  v = p - point();
  const RT  vx = v.hx();
  const RT  vy = v.hy();
  const RT  vz = v.hz();
  const RT  vw = v.hw();
  VectorH3<FT,RT> dir = direction().to_vector();
  const RT  dx = dir.hx();
  const RT  dy = dir.hy();
  const RT  dz = dir.hz();
  const RT  dw = dir.hw();

  RT lambda_num = (vx*dx + vy*dy + vz*dz)*dw; // *dw
  RT lambda_den = (dx*dx + dy*dy + dz*dz)*vw; // *dw

  return point() + ( (lambda_num * dir)/lambda_den );
}

template < class FT, class RT >
CGAL_KERNEL_INLINE
LineH3<FT,RT>
LineH3<FT,RT>::transform(const Aff_transformationH3<FT,RT>& t) const
{ return LineH3<FT,RT>(t.transform(point() ), t.transform(direction() )); }


#ifndef NO_OSTREAM_INSERT_LINEH3
template < class FT, class RT >
std::ostream &operator<<(std::ostream &os, const LineH3<FT,RT> &l)
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
#endif // NO_OSTREAM_INSERT_LINEH3

#ifndef NO_ISTREAM_EXTRACT_LINEH3
template < class FT, class RT >
std::istream &operator>>(std::istream &is, LineH3<FT,RT> &l)
{
  PointH3<FT,RT> p;
  DirectionH3<FT,RT> d;
  is >> p >> d;
  l = LineH3<FT,RT>(p, d);
  return is;
}
#endif // NO_ISTREAM_EXTRACT_LINEH3
template < class FT, class RT >
CGAL_KERNEL_INLINE
bool
LineH3<FT,RT>::has_on( const PointH3<FT,RT>& p ) const
{ return collinear(point(), point()+direction().to_vector(), p); }

template < class FT, class RT >
CGAL_KERNEL_INLINE
bool
LineH3<FT,RT>::is_degenerate() const
{ return direction().is_degenerate(); }

template < class FT, class RT >
CGAL_KERNEL_INLINE
bool
LineH3<FT,RT>::operator==(const LineH3<FT,RT>& l) const
{
  return  (  (l.direction() ==   ptr->direction )
           &&(l.has_on( ptr->basepoint ) ) );
}

CGAL_END_NAMESPACE


#endif // CGAL_LINEH3_H
