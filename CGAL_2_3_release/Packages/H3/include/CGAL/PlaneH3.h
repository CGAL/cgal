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
// file          : PlaneH3.h
// package       : H3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_PLANEH3_H
#define CGAL_PLANEH3_H

#include <CGAL/PointH2.h>
#include <CGAL/PointH3.h>
#include <CGAL/LineH3.h>
#include <CGAL/RayH3.h>
#include <CGAL/SegmentH3.h>
#include <CGAL/basic_constructionsH3.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class PlaneH3
  : public R_::Plane_handle_3
{
  public:
    typedef R_                 R;
    typedef typename R::RT     RT;
    typedef typename R::FT     FT;

    typedef typename R::Plane_handle_3            Plane_handle_3_;
    typedef typename Plane_handle_3_::element_type Plane_ref_3;

    PlaneH3()
      : Plane_handle_3_(Plane_ref_3()) {}

    PlaneH3(const PointH3<R>&,
            const PointH3<R>&,
            const PointH3<R>& );
    PlaneH3(const RT& a, const RT& b,
            const RT& c, const RT& d );
    PlaneH3(const PointH3<R>&, const RayH3<R>& );
    PlaneH3(const PointH3<R>&, const LineH3<R>& );
    PlaneH3(const PointH3<R>&, const SegmentH3<R>& );
    PlaneH3(const LineH3<R>&, const PointH3<R>& );
    PlaneH3(const SegmentH3<R>&, const PointH3<R>& );
    PlaneH3(const RayH3<R>&, const PointH3<R>& );
    PlaneH3(const PointH3<R>&, const DirectionH3<R>& );
    PlaneH3(const PointH3<R>&, const VectorH3<R>& );
    PlaneH3(const PointH3<R>&, const DirectionH3<R>&, const DirectionH3<R>& );

    RT             a() const;
    RT             b() const;
    RT             c() const;
    RT             d() const;

    bool           operator==( const PlaneH3<R>& ) const;
    bool           operator!=( const PlaneH3<R>& ) const;

    LineH3<R>  perpendicular_line(const PointH3<R>& ) const;
    PlaneH3<R> opposite() const;  // plane with opposite orientation
    PointH3<R> projection(const PointH3<R>& ) const;

    PointH3<R> point() const;     // same point on the plane
    DirectionH3<R>
                   orthogonal_direction() const;
    VectorH3<R>
                   orthogonal_vector() const;

    Oriented_side  oriented_side(const PointH3<R> &p) const;
    bool           has_on(const PointH3<R> &p) const;
    bool           has_on(const LineH3<R> &p) const;
#ifndef CGAL_NO_DEPRECATED_CODE
    bool           has_on_boundary(const PointH3<R> &p) const;
    bool           has_on_boundary(const LineH3<R> &p) const;
#endif // CGAL_NO_DEPRECATED_CODE
    bool           has_on_positive_side(const PointH3<R>&l) const;
    bool           has_on_negative_side(const PointH3<R>&l) const;

    bool           is_degenerate() const;

    PlaneH3<R> transform(const Aff_transformationH3<R>& ) const;

    Aff_transformationH3<R> transform_to_2d() const;
    PointH2<R>   to_2d(const PointH3<R>& )  const;
    PointH3<R>   to_3d(const PointH2<R>& )  const;
    VectorH3<R>  base1() const;
    VectorH3<R>  base2() const;


protected:
    PointH3<R>   point1() const;   // same point different from point()
    PointH3<R>   point2() const;   // same point different from point()
                                       // and point1()

    void             new_rep(const PointH3<R> &p,
                             const PointH3<R> &q,
                             const PointH3<R> &r);

    void             new_rep(const RT &a, const RT &b,
                             const RT &c, const RT &d);
};

//
//  a() * X + b() * Y + c() * Z() + d() * W() == 0
//
//      |    X        Y       Z       W     |
//      |  p.hx()   p.hy()  p.hz()  p.hw()  |
//      |  q.hx()   q.hy()  q.hz()  q.hw()  |
//      |  r.hx()   r.hy()  r.hz()  r.hw()  |
//
//  Fourtuple<RT> ( a(), b(), c(), d() )

template < class R >
inline
void
PlaneH3<R>::new_rep(const PointH3<R> &p,
                        const PointH3<R> &q,
                        const PointH3<R> &r)
{
  RT phx = p.hx();
  RT phy = p.hy();
  RT phz = p.hz();
  RT phw = p.hw();

  RT qhx = q.hx();
  RT qhy = q.hy();
  RT qhz = q.hz();
  RT qhw = q.hw();

  RT rhx = r.hx();
  RT rhy = r.hy();
  RT rhz = r.hz();
  RT rhw = r.hw();

  initialize_with( Plane_ref_3 (
              phy*( qhz*rhw - qhw*rhz )
            - qhy*( phz*rhw - phw*rhz )     // * X
            + rhy*( phz*qhw - phw*qhz ),

            - phx*( qhz*rhw - qhw*rhz )
            + qhx*( phz*rhw - phw*rhz )     // * Y
            - rhx*( phz*qhw - phw*qhz ),

              phx*( qhy*rhw - qhw*rhy )
            - qhx*( phy*rhw - phw*rhy )     // * Z
            + rhx*( phy*qhw - phw*qhy ),

            - phx*( qhy*rhz - qhz*rhy )
            + qhx*( phy*rhz - phz*rhy )     // * W
            - rhx*( phy*qhz - phz*qhy )          ) );
}

template < class R >
inline
void
PlaneH3<R>::new_rep(const RT &a, const RT &b, const RT &c, const RT &d)
{ initialize_with( Plane_ref_3 (a, b, c, d) ); }

template < class R >
inline
bool
PlaneH3<R>::operator!=(const PlaneH3<R>& l) const
{
 return !(*this == l);
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<R>::PlaneH3(const PointH3<R>& p,
                        const PointH3<R>& q,
                        const PointH3<R>& r)
{ new_rep(p,q,r); }

template < class R >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<R>::PlaneH3(const RT& a, const RT& b,
                        const RT& c, const RT& d)
{ new_rep(a,b,c,d); }

template < class R >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<R>::PlaneH3(const PointH3<R>& p ,
                        const LineH3<R>&  l)
{ new_rep(p, l.point(0), l.point(1) ); }

template < class R >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<R>::PlaneH3(const PointH3<R>& p,
                        const SegmentH3<R>& s)
{ new_rep(p, s.source(), s.target() ); }

template < class R >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<R>::PlaneH3(const PointH3<R>& p ,
                        const RayH3<R>&  r)
{ new_rep(p, r.start(), r.start() + r.direction().to_vector() ); }

template < class R >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<R>::PlaneH3(const LineH3<R>& l ,
                        const PointH3<R>& p)
{ new_rep(l.point(0), p, l.point(1) ); }

template < class R >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<R>::PlaneH3(const SegmentH3<R>& s,
                        const PointH3<R>& p)
{ new_rep(s.source(), p, s.target() ); }

template < class R >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<R>::PlaneH3(const RayH3<R>&  r,
                        const PointH3<R>& p)
{ new_rep(r.start(), p, r.start() + r.direction().to_vector() ); }

template < class R >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<R>::PlaneH3(const PointH3<R>& p,
                        const DirectionH3<R>& d)
{
  VectorH3<R> ov = d.to_vector();
  new_rep( ov.hx()*p.hw(),
           ov.hy()*p.hw(),
           ov.hz()*p.hw(),
          -(ov.hx()*p.hx() + ov.hy()*p.hy() + ov.hz()*p.hz() ) );
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<R>::PlaneH3(const PointH3<R>& p,
                        const VectorH3<R>& ov)
{
  new_rep( ov.hx()*p.hw(),
           ov.hy()*p.hw(),
           ov.hz()*p.hw(),
          -(ov.hx()*p.hx() + ov.hy()*p.hy() + ov.hz()*p.hz() ) );
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<R>::PlaneH3(const PointH3<R>& p,
                        const DirectionH3<R>& d1,
                        const DirectionH3<R>& d2)
{ new_rep( p, p + d1.to_vector(), p + d2.to_vector() ); }

template < class R >
inline
typename PlaneH3<R>::RT
PlaneH3<R>::a() const
{ return Ptr()->e0; }

template < class R >
inline
typename PlaneH3<R>::RT
PlaneH3<R>::b() const
{ return Ptr()->e1; }

template < class R >
inline
typename PlaneH3<R>::RT
PlaneH3<R>::c() const
{ return Ptr()->e2; }

template < class R >
inline
typename PlaneH3<R>::RT
PlaneH3<R>::d() const
{ return Ptr()->e3; }

template < class R >
CGAL_KERNEL_INLINE
LineH3<R>
PlaneH3<R>::perpendicular_line(const PointH3<R>& p) const
{ return LineH3<R>( p, orthogonal_direction() ); }

template < class R >
CGAL_KERNEL_INLINE
PlaneH3<R>
PlaneH3<R>::opposite() const
{ return PlaneH3<R>(-a(), -b(), -c(), -d() ); }

template < class R >
CGAL_KERNEL_INLINE
PointH3<R>
PlaneH3<R>::projection(const PointH3<R>& p) const
{ return _projection( p, *this ); }

template < class R >
CGAL_KERNEL_INLINE
PointH3<R>
PlaneH3<R>::point() const
{
  const RT RT0(0);
  if ( a() != RT0 )
  {
      return PointH3<R>( -d(), RT0, RT0, a() );
  }
  if ( b() != RT0 )
  {
      return PointH3<R>( RT0, -d(), RT0, b() );
  }
  CGAL_kernel_assertion ( c() != RT0);
  return PointH3<R>( RT0, RT0, -d(), c() );
}

template < class R >
CGAL_KERNEL_INLINE
VectorH3<R>
PlaneH3<R>::base1() const
{
 // point():
 // a() != RT0 : PointH3<R>( -d(), RT0, RT0, a() );
 // b() != RT0 : PointH3<R>( RT0, -d(), RT0, b() );
 //            : PointH3<R>( RT0, RT0, -d(), c() );
 // point1():
 // a() != RT0 : PointH3<R>( -b()-d(), a(), RT0, a() );
 // b() != RT0 : PointH3<R>( RT0, -c()-d(), b(), b() );
 //            : PointH3<R>( c(), RT0, -a()-d(), c() );

  const RT RT0(0);
  if ( a() != RT0 )
  {
      return VectorH3<R>( -b(), a(), RT0, a() );
  }
  if ( b() != RT0 )
  {
      return VectorH3<R>( RT0, -c(), b(), b() );
  }
  CGAL_kernel_assertion ( c() != RT(0) );
  return VectorH3<R>( c(), RT0, -a(), c() );
}

template < class R >
inline
VectorH3<R>
PlaneH3<R>::base2() const
{
  VectorH3<R> a = orthogonal_vector();
  VectorH3<R> b = base1();
  return VectorH3<R>(a.hy()*b.hz() - a.hz()*b.hy(),
                         a.hz()*b.hx() - a.hx()*b.hz(),
                         a.hx()*b.hy() - a.hy()*b.hx(),
                         a.hw()*b.hw() );
}
// Actually, the following should work, but bcc doesn't like it:
// { return cross_product( orthogonal_vector(), base1() ); }


template < class R >
inline
PointH3<R>
PlaneH3<R>::point1() const
{ return point() + base1(); }

template < class R >
inline
PointH3<R>
PlaneH3<R>::point2() const
{ return point() + base2(); }

template < class R >
inline
DirectionH3<R>
PlaneH3<R>::orthogonal_direction() const
{ return DirectionH3<R>(a(), b(), c() ); }

template < class R >
inline
VectorH3<R>
PlaneH3<R>::orthogonal_vector() const
{ return VectorH3<R>(a(), b(), c() ); }

template < class R >
PlaneH3<R>
PlaneH3<R>::transform(const Aff_transformationH3<R>& t) const
{
 return t.transform(*this);
}

#ifndef CGAL_NO_OSTREAM_INSERT_PLANE3
template < class R >
std::ostream &operator<<(std::ostream &os, const PlaneH3<R> &p)
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
#endif // CGAL_NO_OSTREAM_INSERT_PLANE3

#ifndef CGAL_NO_ISTREAM_EXTRACT_PLANE3
template < class R  >
std::istream &operator>>(std::istream &is, PlaneH3<R> &p)
{
    typename R::RT a, b, c, d;
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
    p = PlaneH3<R>(a, b, c, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_PLANE3

template < class R >
bool
PlaneH3<R>::is_degenerate() const
{
 const RT RT0(0);
 return ( (a() == RT0 ) && (b() == RT0 ) && (c() == RT0 ) );
}

template < class R >
bool
PlaneH3<R>::has_on_positive_side( const PointH3<R>& p) const
{
 return (a()*p.hx() + b()*p.hy() + c()*p.hz() + d()*p.hw() > RT(0) );
}

template < class R >
bool
PlaneH3<R>::has_on_negative_side( const PointH3<R>& p) const
{
 return (a()*p.hx() + b()*p.hy() + c()*p.hz() + d()*p.hw() < RT(0) );
}


template < class R >
bool
PlaneH3<R>::has_on( const PointH3<R>& p) const
{
 return (a()*p.hx() + b()*p.hy() + c()*p.hz() + d()*p.hw() == RT(0) );
}

template < class R >
bool
PlaneH3<R>::has_on( const LineH3<R>& l) const
{
 PointH3<R>   p   = l.point();
 VectorH3<R>  ld  = l.direction().to_vector();
 VectorH3<R>  ov  = orthogonal_vector();

 return (  ( a()*p.hx() + b()*p.hy() + c()*p.hz() + d()*p.hw()   == RT(0) )
         &&( ld.hx()*ov.hx() + ld.hy()*ov.hy() + ld.hz()*ov.hz() == RT(0) ) );
}

#ifndef CGAL_NO_DEPRECATED_CODE
template < class R >
bool
PlaneH3<R>::has_on_boundary( const PointH3<R>& p) const
{
 return has_on(p);
}

template < class R >
bool
PlaneH3<R>::has_on_boundary( const LineH3<R>& l) const
{
 return has_on(l);
}
#endif

template < class R >
Oriented_side
PlaneH3<R>::oriented_side( const PointH3<R>& p) const
{
 RT value = a()*p.hx() + b()*p.hy() + c()*p.hz() + d()*p.hw() ;
 if (value > RT(0) )
 {
    return ON_POSITIVE_SIDE;
 }
 else
 {
    return
    (value < RT(0) ) ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
 }
}


template < class R >
bool
PlaneH3<R>::operator==(const PlaneH3<R>& l) const
{
 if (  (a() * l.d() != l.a() * d() )
     ||(b() * l.d() != l.b() * d() )
     ||(c() * l.d() != l.c() * d() ) )
 {
    return false;
 }
 int sd  = static_cast<int>(CGAL_NTS sign(d()));
 int sld = static_cast<int>(CGAL_NTS sign(l.d()));
 if ( sd == sld )
 {
    if (sd == 0)
    {
        return (  (a()*l.b() == b()*l.a() )
                &&(a()*l.c() == c()*l.a() )
                &&(b()*l.c() == c()*l.b() )
                &&(CGAL_NTS sign(a() )== CGAL_NTS sign( l.a() ))
                &&(CGAL_NTS sign(b() )== CGAL_NTS sign( l.b() ))
                &&(CGAL_NTS sign(c() )== CGAL_NTS sign( l.c() )) );
    }
    else
    {
        return true;
    }
 }
 else
 {
    return false;
 }
}


CGAL_END_NAMESPACE

#include <CGAL/Aff_transformationH3.h>

CGAL_BEGIN_NAMESPACE

template < class R >
Aff_transformationH3<R>
PlaneH3<R>::transform_to_2d() const
{
  const RT  RT0(0);
  const RT  RT1(1);
  VectorH3<R> nov = orthogonal_vector();
  VectorH3<R> e1v = point1()-point() ;
  VectorH3<R> e2v = point2()-point() ;
  RT orthohx = nov.hx();
  RT orthohy = nov.hy();
  RT orthohz = nov.hz();
  RT e1phx   = e1v.hx();
  RT e1phy   = e1v.hy();
  RT e1phz   = e1v.hz();
  RT e2phx   = e2v.hx();
  RT e2phy   = e2v.hy();
  RT e2phz   = e2v.hz();

  RT t11 =  -( orthohy*e2phz - orthohz*e2phy );
  RT t12 =   ( orthohx*e2phz - orthohz*e2phx );
  RT t13 =  -( orthohx*e2phy - orthohy*e2phx );

  RT t21 =   ( orthohy*e1phz - orthohz*e1phy );
  RT t22 =  -( orthohx*e1phz - orthohz*e1phx );
  RT t23 =   ( orthohx*e1phy - orthohy*e1phx );

  RT t31 =   ( e1phy*e2phz - e1phz*e2phy );
  RT t32 =  -( e1phx*e2phz - e1phz*e2phx );
  RT t33 =   ( e1phx*e2phy - e1phy*e2phx );

  RT scale = det3x3_by_formula( orthohx, orthohy, orthohz,
                                     e1phx,   e1phy,   e1phz,
                                     e2phx,   e2phy,   e2phz );

  Aff_transformationH3<R>
     point_to_origin(TRANSLATION,  - ( point() - ORIGIN ) );
  Aff_transformationH3<R>
     rotate_and_more( t11,    t12,   t13,   RT0,
                      t21,    t22,   t23,   RT0,
                      t31,    t32,   t33,   RT0,
                                            scale);

  PointH3<R> ortho( orthohx, orthohy, orthohz );
  PointH3<R> e1p( e1phx, e1phy, e1phz );
  PointH3<R> e2p( e2phx, e2phy, e2phz );
  CGAL_kernel_assertion((   ortho.transform(rotate_and_more)
        == PointH3<R>( RT(0), RT(0), RT(1)) ));
  CGAL_kernel_assertion((   e1p.transform(rotate_and_more)
        == PointH3<R>( RT(1), RT(0), RT(0)) ));
  CGAL_kernel_assertion((   e2p.transform(rotate_and_more)
        == PointH3<R>( RT(0), RT(1), RT(0)) ));

  return  rotate_and_more * point_to_origin;
}

template < class R >
CGAL_KERNEL_INLINE
PointH2<R>
PlaneH3<R>::to_2d(const PointH3<R>& p) const
{
  PointH3<R> tp = p.transform( transform_to_2d() );
  return PointH2<R>( tp.hx(), tp.hy(), tp.hw());
}


template < class R >
CGAL_KERNEL_INLINE
PointH3<R>
PlaneH3<R>::to_3d(const PointH2<R>& p)  const
{
  PointH3<R> hp( p.hx(), p.hy(), RT(0.0), p.hw());
  return hp.transform( transform_to_2d().inverse() );
}

CGAL_END_NAMESPACE

#endif  // CGAL_PLANEH3_H
