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

#if defined(CGAL_CFG_INCOMPLETE_TYPE_BUG_1) && \
   !defined(CGAL_NO_PLANE_TRANSFORM_IN_AT)
#define CGAL_NO_PLANE_TRANSFORM_IN_AT
#endif // CGAL_CFG_INCOMPLETE_TYPE_BUG_1

#include <CGAL/PointH2.h>
#include <CGAL/PointH3.h>
#include <CGAL/LineH3.h>
#include <CGAL/RayH3.h>
#include <CGAL/SegmentH3.h>
#include <CGAL/basic_constructionsH3.h>

CGAL_BEGIN_NAMESPACE

template < class FT, class RT >
class PlaneH3 : public Handle_for< Fourtuple<RT> >
{
  public:
    PlaneH3();
    PlaneH3(const PointH3<FT,RT>& ,
            const PointH3<FT,RT>& ,
            const PointH3<FT,RT>& );
    PlaneH3(const RT& a, const RT& b,
            const RT& c, const RT& d );
    PlaneH3(const PointH3<FT,RT>& ,
            const RayH3<FT,RT>& );
    PlaneH3(const PointH3<FT,RT>& ,
            const LineH3<FT,RT>& );
    PlaneH3(const PointH3<FT,RT>& ,
            const SegmentH3<FT,RT>& );
    PlaneH3(const LineH3<FT,RT>& ,
            const PointH3<FT,RT>& );
    PlaneH3(const SegmentH3<FT,RT>& ,
            const PointH3<FT,RT>& );
    PlaneH3(const RayH3<FT,RT>& ,
            const PointH3<FT,RT>& );
    PlaneH3(const PointH3<FT,RT>&,
            const DirectionH3<FT,RT>& );
    PlaneH3(const PointH3<FT,RT>&,
            const VectorH3<FT,RT>& );
    PlaneH3(const PointH3<FT,RT>&,
            const DirectionH3<FT,RT>&,
            const DirectionH3<FT,RT>& );

    RT             a() const;
    RT             b() const;
    RT             c() const;
    RT             d() const;

    bool           operator==( const PlaneH3<FT,RT>& ) const;
    bool           operator!=( const PlaneH3<FT,RT>& ) const;

    LineH3<FT,RT>  perpendicular_line(const PointH3<FT,RT>& ) const;
    PlaneH3<FT,RT> opposite() const;  // plane with opposite orientation
    PointH3<FT,RT> projection(const PointH3<FT,RT>& ) const;

    PointH3<FT,RT> point() const;     // same point on the plane
    DirectionH3<FT,RT>
                   orthogonal_direction() const;
    VectorH3<FT,RT>
                   orthogonal_vector() const;

    Oriented_side  oriented_side(const PointH3<FT,RT> &p) const;
    bool           has_on(const PointH3<FT,RT> &p) const;
    bool           has_on(const LineH3<FT,RT> &p) const;
    bool           has_on_boundary(const PointH3<FT,RT> &p) const;
    bool           has_on_boundary(const LineH3<FT,RT> &p) const;
    bool           has_on_positive_side(const PointH3<FT,RT>&l) const;
    bool           has_on_negative_side(const PointH3<FT,RT>&l) const;

    bool           is_degenerate() const;

    PlaneH3<FT,RT> transform(const Aff_transformationH3<FT,RT>& ) const;


    Aff_transformationH3<FT,RT>
                     transform_to_2d() const;
    PointH2<FT,RT>   to_2d(const PointH3<FT,RT>& )  const;
    PointH3<FT,RT>   to_3d(const PointH2<FT,RT>& )  const;
    VectorH3<FT,RT>  base1() const;
    VectorH3<FT,RT>  base2() const;


protected:
    PointH3<FT,RT>   point1() const;   // same point different from point()
    PointH3<FT,RT>   point2() const;   // same point different from point()
                                       // and point1()

    void             new_rep(const PointH3<FT,RT> &p,
                             const PointH3<FT,RT> &q,
                             const PointH3<FT,RT> &r);

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

template < class FT, class RT >
inline
void
PlaneH3<FT,RT>::new_rep(const PointH3<FT,RT> &p,
                        const PointH3<FT,RT> &q,
                        const PointH3<FT,RT> &r)
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

  initialize_with( Fourtuple<RT> (
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

template < class FT, class RT >
inline
void
PlaneH3<FT,RT>::new_rep(const RT &a, const RT &b, const RT &c, const RT &d)
{ initialize_with( Fourtuple<RT>(a, b, c, d) ); }

template < class FT, class RT >
inline
bool
PlaneH3<FT,RT>::operator!=(const PlaneH3<FT,RT>& l) const
{
 return !(*this == l);
}



template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<FT,RT>::PlaneH3()
 : Handle_for< Fourtuple<RT> >(Fourtuple<RT>() )
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<FT,RT>::PlaneH3(const PointH3<FT,RT>& p,
                        const PointH3<FT,RT>& q,
                        const PointH3<FT,RT>& r)
{ new_rep(p,q,r); }

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<FT,RT>::PlaneH3(const RT& a, const RT& b,
                        const RT& c, const RT& d)
{ new_rep(a,b,c,d); }

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<FT,RT>::PlaneH3(const PointH3<FT,RT>& p ,
                        const LineH3<FT,RT>&  l)
{ new_rep(p, l.point(0), l.point(1) ); }

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<FT,RT>::PlaneH3(const PointH3<FT,RT>& p,
                        const SegmentH3<FT,RT>& s)
{ new_rep(p, s.source(), s.target() ); }

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<FT,RT>::PlaneH3(const PointH3<FT,RT>& p ,
                        const RayH3<FT,RT>&  r)
{ new_rep(p, r.start(), r.start() + r.direction().to_vector() ); }

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<FT,RT>::PlaneH3(const LineH3<FT,RT>& l ,
                        const PointH3<FT,RT>& p)
{ new_rep(l.point(0), p, l.point(1) ); }

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<FT,RT>::PlaneH3(const SegmentH3<FT,RT>& s,
                        const PointH3<FT,RT>& p)
{ new_rep(s.source(), p, s.target() ); }

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<FT,RT>::PlaneH3(const RayH3<FT,RT>&  r,
                        const PointH3<FT,RT>& p)
{ new_rep(r.start(), p, r.start() + r.direction().to_vector() ); }

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<FT,RT>::PlaneH3(const PointH3<FT,RT>& p,
                        const DirectionH3<FT,RT>& d)
{
  VectorH3<FT,RT> ov = d.to_vector();
  new_rep( ov.hx()*p.hw(),
           ov.hy()*p.hw(),
           ov.hz()*p.hw(),
          -(ov.hx()*p.hx() + ov.hy()*p.hy() + ov.hz()*p.hz() ) );
}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<FT,RT>::PlaneH3(const PointH3<FT,RT>& p,
                        const VectorH3<FT,RT>& ov)
{
  new_rep( ov.hx()*p.hw(),
           ov.hy()*p.hw(),
           ov.hz()*p.hw(),
          -(ov.hx()*p.hx() + ov.hy()*p.hy() + ov.hz()*p.hz() ) );
}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PlaneH3<FT,RT>::PlaneH3(const PointH3<FT,RT>& p,
                        const DirectionH3<FT,RT>& d1,
                        const DirectionH3<FT,RT>& d2)
{ new_rep( p, p + d1.to_vector(), p + d2.to_vector() ); }

template < class FT, class RT >
inline
RT
PlaneH3<FT,RT>::a() const
{ return ptr->e0; }

template < class FT, class RT >
inline
RT
PlaneH3<FT,RT>::b() const
{ return ptr->e1; }

template < class FT, class RT >
inline
RT
PlaneH3<FT,RT>::c() const
{ return ptr->e2; }

template < class FT, class RT >
inline
RT
PlaneH3<FT,RT>::d() const
{ return ptr->e3; }

template < class FT, class RT >
CGAL_KERNEL_INLINE
LineH3<FT,RT>
PlaneH3<FT,RT>::perpendicular_line(const PointH3<FT,RT>& p) const
{ return LineH3<FT,RT>( p, orthogonal_direction() ); }

template < class FT, class RT >
CGAL_KERNEL_INLINE
PlaneH3<FT,RT>
PlaneH3<FT,RT>::opposite() const
{ return PlaneH3<FT,RT>(-a(), -b(), -c(), -d() ); }

template < class FT, class RT >
CGAL_KERNEL_INLINE
PointH3<FT,RT>
PlaneH3<FT,RT>::projection(const PointH3<FT,RT>& p) const
{ return _projection( p, *this ); }

template < class FT, class RT >
CGAL_KERNEL_INLINE
PointH3<FT,RT>
PlaneH3<FT,RT>::point() const
{
  const RT RT0(0);
  if ( a() != RT0 )
  {
      return PointH3<FT,RT>( -d(), RT0, RT0, a() );
  }
  if ( b() != RT0 )
  {
      return PointH3<FT,RT>( RT0, -d(), RT0, b() );
  }
  CGAL_kernel_assertion ( c() != RT0);
  return PointH3<FT,RT>( RT0, RT0, -d(), c() );
}

template < class FT, class RT >
CGAL_KERNEL_INLINE
VectorH3<FT,RT>
PlaneH3<FT,RT>::base1() const
{
 // point():
 // a() != RT0 : PointH3<FT,RT>( -d(), RT0, RT0, a() );
 // b() != RT0 : PointH3<FT,RT>( RT0, -d(), RT0, b() );
 //            : PointH3<FT,RT>( RT0, RT0, -d(), c() );
 // point1():
 // a() != RT0 : PointH3<FT,RT>( -b()-d(), a(), RT0, a() );
 // b() != RT0 : PointH3<FT,RT>( RT0, -c()-d(), b(), b() );
 //            : PointH3<FT,RT>( c(), RT0, -a()-d(), c() );

  const RT RT0(0);
  if ( a() != RT0 )
  {
      return VectorH3<FT,RT>( -b(), a(), RT0, a() );
  }
  if ( b() != RT0 )
  {
      return VectorH3<FT,RT>( RT0, -c(), b(), b() );
  }
  CGAL_kernel_assertion ( c() != RT(0) );
  return VectorH3<FT,RT>( c(), RT0, -a(), c() );
}

template < class FT, class RT >
inline
VectorH3<FT,RT>
PlaneH3<FT,RT>::base2() const
{
  VectorH3<FT,RT> a = orthogonal_vector();
  VectorH3<FT,RT> b = base1();
  return VectorH3<FT,RT>(a.hy()*b.hz() - a.hz()*b.hy(),
                         a.hz()*b.hx() - a.hx()*b.hz(),
                         a.hx()*b.hy() - a.hy()*b.hx(),
                         a.hw()*b.hw() );
}
// Actually, the following should work, but bcc doesn't like it:
// { return cross_product( orthogonal_vector(), base1() ); }


template < class FT, class RT >
inline
PointH3<FT,RT>
PlaneH3<FT,RT>::point1() const
{ return point() + base1(); }

template < class FT, class RT >
inline
PointH3<FT,RT>
PlaneH3<FT,RT>::point2() const
{ return point() + base2(); }

template < class FT, class RT >
inline
DirectionH3<FT,RT>
PlaneH3<FT,RT>::orthogonal_direction() const
{ return DirectionH3<FT,RT>(a(), b(), c() ); }

template < class FT, class RT >
inline
VectorH3<FT,RT>
PlaneH3<FT,RT>::orthogonal_vector() const
{ return VectorH3<FT,RT>(a(), b(), c() ); }

template < class FT, class RT >
PlaneH3<FT,RT>
PlaneH3<FT,RT>::transform(const Aff_transformationH3<FT,RT>& t) const
{
#ifndef CGAL_NO_PLANE_TRANSFORM_IN_AT
 return t.transform(*this);
#else
 if ( t.is_even() )
 {
     return PlaneH3<FT,RT>(
             t.transform(point() ),
             t.transpose().inverse().transform(orthogonal_direction() ));
 }
 else
 {
     return PlaneH3<FT,RT>(
             t.transform(point() ),
           - t.transpose().inverse().transform(orthogonal_direction() ));
 }
#endif // CGAL_NO_PLANE_TRANSFORM_IN_AT
}



#ifndef CGAL_NO_OSTREAM_INSERT_PLANE3
template < class FT, class RT >
std::ostream &operator<<(std::ostream &os, const PlaneH3<FT,RT> &p)
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
template < class FT, class RT  >
std::istream &operator>>(std::istream &is, PlaneH3<FT,RT> &p)
{
    RT a, b, c, d;
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
    p = PlaneH3<FT,RT>(a, b, c, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_PLANE3

template < class FT, class RT >
bool
PlaneH3<FT,RT>::is_degenerate() const
{
 const RT RT0(0);
 return ( (a() == RT0 ) && (b() == RT0 ) && (c() == RT0 ) );
}

template < class FT, class RT >
bool
PlaneH3<FT,RT>::has_on_positive_side( const PointH3<FT,RT>& p) const
{
 return (a()*p.hx() + b()*p.hy() + c()*p.hz() + d()*p.hw() > RT(0) );
}

template < class FT, class RT >
bool
PlaneH3<FT,RT>::has_on_negative_side( const PointH3<FT,RT>& p) const
{
 return (a()*p.hx() + b()*p.hy() + c()*p.hz() + d()*p.hw() < RT(0) );
}

template < class FT, class RT >
bool
PlaneH3<FT,RT>::has_on_boundary( const PointH3<FT,RT>& p) const
{
 return (a()*p.hx() + b()*p.hy() + c()*p.hz() + d()*p.hw() == RT(0) );
}

template < class FT, class RT >
bool
PlaneH3<FT,RT>::has_on( const PointH3<FT,RT>& p) const
{
 return has_on_boundary(p);
}

template < class FT, class RT >
bool
PlaneH3<FT,RT>::has_on_boundary( const LineH3<FT,RT>& l) const
{
 PointH3<FT,RT>   p   = l.point();
 VectorH3<FT,RT>  ld  = l.direction().to_vector();
 VectorH3<FT,RT>  ov  = orthogonal_vector();

 return (  ( a()*p.hx() + b()*p.hy() + c()*p.hz() + d()*p.hw()   == RT(0) )
         &&( ld.hx()*ov.hx() + ld.hy()*ov.hy() + ld.hz()*ov.hz() == RT(0) ) );
}

template < class FT, class RT >
bool
PlaneH3<FT,RT>::has_on( const LineH3<FT,RT>& l) const
{
 return has_on_boundary(l);
}

template < class FT, class RT >
Oriented_side
PlaneH3<FT,RT>::oriented_side( const PointH3<FT,RT>& p) const
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


template < class FT, class RT >
bool
PlaneH3<FT,RT>::operator==(const PlaneH3<FT,RT>& l) const
{
 if (  (a() * l.d() != l.a() * d() )
     ||(b() * l.d() != l.b() * d() )
     ||(c() * l.d() != l.c() * d() ) )
 {
    return false;
 }
 int sd  = (int)CGAL_NTS sign(d() );
 int sld = (int)CGAL_NTS sign(l.d() );
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

template < class FT, class RT >
Aff_transformationH3<FT,RT>
PlaneH3<FT,RT>::transform_to_2d() const
{
  const RT  RT0(0);
  const RT  RT1(1);
  VectorH3<FT,RT> nov = orthogonal_vector();
  VectorH3<FT,RT> e1v = point1()-point() ;
  VectorH3<FT,RT> e2v = point2()-point() ;
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

  Aff_transformationH3<FT,RT>
     point_to_origin(TRANSLATION,  - ( point() - ORIGIN ) );
  Aff_transformationH3<FT,RT>
     rotate_and_more( t11,    t12,   t13,   RT0,
                      t21,    t22,   t23,   RT0,
                      t31,    t32,   t33,   RT0,
                                            scale);

  PointH3<FT,RT> ortho( orthohx, orthohy, orthohz );
  PointH3<FT,RT> e1p( e1phx, e1phy, e1phz );
  PointH3<FT,RT> e2p( e2phx, e2phy, e2phz );
  CGAL_kernel_assertion((   ortho.transform(rotate_and_more)
        == PointH3<FT,RT>( RT(0), RT(0), RT(1)) ));
  CGAL_kernel_assertion((   e1p.transform(rotate_and_more)
        == PointH3<FT,RT>( RT(1), RT(0), RT(0)) ));
  CGAL_kernel_assertion((   e2p.transform(rotate_and_more)
        == PointH3<FT,RT>( RT(0), RT(1), RT(0)) ));

  return  rotate_and_more * point_to_origin;
}

template < class FT, class RT >
CGAL_KERNEL_INLINE
PointH2<FT,RT>
PlaneH3<FT,RT>::to_2d(const PointH3<FT,RT>& p) const
{
  PointH3<FT,RT> tp = p.transform( transform_to_2d() );
  return PointH2<FT,RT>( tp.hx(), tp.hy(), tp.hw());
}


template < class FT, class RT >
CGAL_KERNEL_INLINE
PointH3<FT,RT>
PlaneH3<FT,RT>::to_3d(const PointH2<FT,RT>& p)  const
{
  PointH3<FT,RT> hp( p.hx(), p.hy(), RT(0.0), p.hw());
  return hp.transform( transform_to_2d().inverse() );
}

CGAL_END_NAMESPACE

#endif  // CGAL_PLANEH3_H
