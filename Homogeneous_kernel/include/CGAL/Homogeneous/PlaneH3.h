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
// Author(s)     : Stefan Schirra

#ifndef CGAL_PLANEH3_H
#define CGAL_PLANEH3_H

#include <CGAL/array.h>

namespace CGAL {

template < class R_ >
class PlaneH3
{
   typedef typename R_::RT                   RT;
   typedef typename R_::FT                   FT;
   typedef typename R_::Point_2              Point_2;
   typedef typename R_::Point_3              Point_3;
   typedef typename R_::Vector_3             Vector_3;
   typedef typename R_::Line_3               Line_3;
   typedef typename R_::Segment_3            Segment_3;
   typedef typename R_::Ray_3                Ray_3;
   typedef typename R_::Direction_3          Direction_3;
   typedef typename R_::Plane_3              Plane_3;
   typedef typename R_::Aff_transformation_3 Aff_transformation_3;

   typedef cpp11::array<RT, 4>               Rep;
   typedef typename R_::template Handle<Rep>::type  Base;

   Base base;

public:

   typedef R_                 R;

    PlaneH3() {}

    PlaneH3(const Point_3&, const Point_3&, const Point_3& );
    PlaneH3(const RT& a, const RT& b,
            const RT& c, const RT& d );
    PlaneH3(const Point_3&, const Ray_3& );
    PlaneH3(const Point_3&, const Line_3& );
    PlaneH3(const Point_3&, const Segment_3& );
    PlaneH3(const Line_3&, const Point_3& );
    PlaneH3(const Segment_3&, const Point_3& );
    PlaneH3(const Ray_3&, const Point_3& );
    PlaneH3(const Point_3&, const Direction_3& );
    PlaneH3(const Point_3&, const Vector_3& );
    PlaneH3(const Point_3&, const Direction_3&, const Direction_3& );

    const RT & a() const;
    const RT & b() const;
    const RT & c() const;
    const RT & d() const;

    bool       operator==( const PlaneH3<R>& ) const;
    bool       operator!=( const PlaneH3<R>& ) const;

    Line_3  perpendicular_line(const Point_3& ) const;
    Plane_3 opposite() const;  // plane with opposite orientation
    Point_3 projection(const Point_3& ) const;

    Point_3 point() const;     // same point on the plane
    Direction_3    orthogonal_direction() const;
    Vector_3       orthogonal_vector() const;

    Oriented_side  oriented_side(const Point_3 &p) const;
    bool           has_on(const Point_3 &p) const;
    bool           has_on(const Line_3 &p) const;
    bool           has_on_positive_side(const Point_3&l) const;
    bool           has_on_negative_side(const Point_3&l) const;

    bool           is_degenerate() const;

    Aff_transformation_3 transform_to_2d() const;
    Point_2   to_2d(const Point_3& )  const;
    Point_3   to_3d(const Point_2& )  const;
    Vector_3  base1() const;
    Vector_3  base2() const;

protected:
    Point_3   point1() const;   // same point different from point()
    Point_3   point2() const;   // same point different from point()
                                       // and point1()

    void             new_rep(const Point_3 &p,
                             const Point_3 &q,
                             const Point_3 &r);

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
//  cpp11::array<RT, 4> ( a(), b(), c(), d() )

template < class R >
inline
void
PlaneH3<R>::new_rep(const typename PlaneH3<R>::Point_3 &p,
                    const typename PlaneH3<R>::Point_3 &q,
                    const typename PlaneH3<R>::Point_3 &r)
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

  base = CGAL::make_array<RT>(
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
            - rhx*( phy*qhz - phz*qhy ));
}

template < class R >
inline
void
PlaneH3<R>::new_rep(const RT &a, const RT &b, const RT &c, const RT &d)
{ base = CGAL::make_array(a, b, c, d); }

template < class R >
inline
bool
PlaneH3<R>::operator!=(const PlaneH3<R>& l) const
{
 return !(*this == l);
}

template < class R >
CGAL_KERNEL_INLINE
PlaneH3<R>::PlaneH3(const typename PlaneH3<R>::Point_3& p,
                    const typename PlaneH3<R>::Point_3& q,
                    const typename PlaneH3<R>::Point_3& r)
{ new_rep(p,q,r); }

template < class R >
CGAL_KERNEL_INLINE
PlaneH3<R>::PlaneH3(const RT& a, const RT& b,
                    const RT& c, const RT& d)
{ new_rep(a,b,c,d); }

template < class R >
CGAL_KERNEL_INLINE
PlaneH3<R>::PlaneH3(const typename PlaneH3<R>::Point_3& p ,
                    const typename PlaneH3<R>::Line_3&  l)
{ new_rep(p, l.point(0), l.point(1) ); }

template < class R >
CGAL_KERNEL_INLINE
PlaneH3<R>::PlaneH3(const typename PlaneH3<R>::Point_3& p,
                        const typename PlaneH3<R>::Segment_3& s)
{ new_rep(p, s.source(), s.target() ); }

template < class R >
CGAL_KERNEL_INLINE
PlaneH3<R>::PlaneH3(const typename PlaneH3<R>::Point_3& p ,
                        const typename PlaneH3<R>::Ray_3&  r)
{ new_rep(p, r.start(), r.start() + r.direction().to_vector() ); }

template < class R >
CGAL_KERNEL_INLINE
PlaneH3<R>::PlaneH3(const typename PlaneH3<R>::Line_3& l ,
                        const typename PlaneH3<R>::Point_3& p)
{ new_rep(l.point(0), p, l.point(1) ); }

template < class R >
CGAL_KERNEL_INLINE
PlaneH3<R>::PlaneH3(const typename PlaneH3<R>::Segment_3& s,
                        const typename PlaneH3<R>::Point_3& p)
{ new_rep(s.source(), p, s.target() ); }

template < class R >
CGAL_KERNEL_INLINE
PlaneH3<R>::PlaneH3(const typename PlaneH3<R>::Ray_3&  r,
                        const typename PlaneH3<R>::Point_3& p)
{ new_rep(r.start(), p, r.start() + r.direction().to_vector() ); }

template < class R >
CGAL_KERNEL_INLINE
PlaneH3<R>::PlaneH3(const typename PlaneH3<R>::Point_3& p,
                        const typename PlaneH3<R>::Direction_3& d)
{
  Vector_3 ov = d.to_vector();
  new_rep( ov.hx()*p.hw(),
           ov.hy()*p.hw(),
           ov.hz()*p.hw(),
          -(ov.hx()*p.hx() + ov.hy()*p.hy() + ov.hz()*p.hz() ) );
}

template < class R >
CGAL_KERNEL_INLINE
PlaneH3<R>::PlaneH3(const typename PlaneH3<R>::Point_3& p,
                        const typename PlaneH3<R>::Vector_3& ov)
{
  new_rep( ov.hx()*p.hw(),
           ov.hy()*p.hw(),
           ov.hz()*p.hw(),
          -(ov.hx()*p.hx() + ov.hy()*p.hy() + ov.hz()*p.hz() ) );
}

template < class R >
CGAL_KERNEL_INLINE
PlaneH3<R>::PlaneH3(const typename PlaneH3<R>::Point_3& p,
                        const typename PlaneH3<R>::Direction_3& d1,
                        const typename PlaneH3<R>::Direction_3& d2)
{ new_rep( p, p + d1.to_vector(), p + d2.to_vector() ); }

template < class R >
inline
const typename PlaneH3<R>::RT &
PlaneH3<R>::a() const
{ return get(base)[0]; }

template < class R >
inline
const typename PlaneH3<R>::RT &
PlaneH3<R>::b() const
{ return get(base)[1]; }

template < class R >
inline
const typename PlaneH3<R>::RT &
PlaneH3<R>::c() const
{ return get(base)[2]; }

template < class R >
inline
const typename PlaneH3<R>::RT &
PlaneH3<R>::d() const
{ return get(base)[3]; }

template < class R >
CGAL_KERNEL_INLINE
typename PlaneH3<R>::Line_3
PlaneH3<R>::perpendicular_line(const typename PlaneH3<R>::Point_3& p) const
{ return Line_3( p, orthogonal_direction() ); }

template < class R >
CGAL_KERNEL_INLINE
typename PlaneH3<R>::Plane_3
PlaneH3<R>::opposite() const
{ return PlaneH3<R>(-a(), -b(), -c(), -d() ); }

template < class R >
CGAL_KERNEL_INLINE
typename PlaneH3<R>::Point_3
PlaneH3<R>::projection(const typename PlaneH3<R>::Point_3& p) const
{ return _projection( p, *this ); }

template < class R >
CGAL_KERNEL_INLINE
typename PlaneH3<R>::Point_3
PlaneH3<R>::point() const
{
  const RT RT0(0);
  if ( a() != RT0 )
  {
      return Point_3( -d(), RT0, RT0, a() );
  }
  if ( b() != RT0 )
  {
      return Point_3( RT0, -d(), RT0, b() );
  }
  CGAL_kernel_assertion ( c() != RT0);
  return Point_3( RT0, RT0, -d(), c() );
}

template < class R >
CGAL_KERNEL_INLINE
typename PlaneH3<R>::Vector_3
PlaneH3<R>::base1() const
{
 // point():
 // a() != RT0 : Point_3( -d(), RT0, RT0, a() );
 // b() != RT0 : Point_3( RT0, -d(), RT0, b() );
 //            : Point_3( RT0, RT0, -d(), c() );
 // point1():
 // a() != RT0 : Point_3( -b()-d(), a(), RT0, a() );
 // b() != RT0 : Point_3( RT0, -c()-d(), b(), b() );
 //            : Point_3( c(), RT0, -a()-d(), c() );

  const RT RT0(0);
  if ( a() != RT0 )
  {
      return Vector_3( -b(), a(), RT0, a() );
  }
  if ( b() != RT0 )
  {
      return Vector_3( RT0, -c(), b(), b() );
  }
  CGAL_kernel_assertion ( c() != RT(0) );
  return Vector_3( c(), RT0, -a(), c() );
}

template < class R >
inline
typename PlaneH3<R>::Vector_3
PlaneH3<R>::base2() const
{
  Vector_3 a = orthogonal_vector();
  Vector_3 b = base1();
  return Vector_3(a.hy()*b.hz() - a.hz()*b.hy(),
                         a.hz()*b.hx() - a.hx()*b.hz(),
                         a.hx()*b.hy() - a.hy()*b.hx(),
                         a.hw()*b.hw() );
}
// Actually, the following should work, but bcc doesn't like it:
// { return cross_product( orthogonal_vector(), base1() ); }


template < class R >
inline
typename PlaneH3<R>::Point_3
PlaneH3<R>::point1() const
{ return point() + base1(); }

template < class R >
inline
typename PlaneH3<R>::Point_3
PlaneH3<R>::point2() const
{ return point() + base2(); }

template < class R >
inline
typename PlaneH3<R>::Direction_3
PlaneH3<R>::orthogonal_direction() const
{ return Direction_3(a(), b(), c() ); }

template < class R >
inline
typename PlaneH3<R>::Vector_3
PlaneH3<R>::orthogonal_vector() const
{ return Vector_3(a(), b(), c() ); }

template < class R >
bool
PlaneH3<R>::is_degenerate() const
{
 const RT RT0(0);
 return ( (a() == RT0 ) && (b() == RT0 ) && (c() == RT0 ) );
}

template < class R >
bool
PlaneH3<R>::has_on_positive_side( const typename PlaneH3<R>::Point_3& p) const
{
 return (a()*p.hx() + b()*p.hy() + c()*p.hz() + d()*p.hw() > RT(0) );
}

template < class R >
bool
PlaneH3<R>::has_on_negative_side( const typename PlaneH3<R>::Point_3& p) const
{
 return (a()*p.hx() + b()*p.hy() + c()*p.hz() + d()*p.hw() < RT(0) );
}


template < class R >
bool
PlaneH3<R>::has_on( const typename PlaneH3<R>::Point_3& p) const
{
 return (a()*p.hx() + b()*p.hy() + c()*p.hz() + d()*p.hw() == RT(0) );
}

template < class R >
bool
PlaneH3<R>::has_on( const typename PlaneH3<R>::Line_3& l) const
{
 Point_3   p   = l.point();
 Vector_3  ld  = l.direction().to_vector();
 Vector_3  ov  = orthogonal_vector();

 return (  ( a()*p.hx() + b()*p.hy() + c()*p.hz() + d()*p.hw()   == RT(0) )
         &&( ld.hx()*ov.hx() + ld.hy()*ov.hy() + ld.hz()*ov.hz() == RT(0) ) );
}

template < class R >
Oriented_side
PlaneH3<R>::oriented_side( const typename PlaneH3<R>::Point_3& p) const
{
 return CGAL_NTS sign( a()*p.hx() + b()*p.hy() + c()*p.hz() + d()*p.hw() );
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

template < class R >
typename PlaneH3<R>::Aff_transformation_3
PlaneH3<R>::transform_to_2d() const
{
  const RT  RT0(0);
  const RT  RT1(1);
  Vector_3 nov = orthogonal_vector();
  Vector_3 e1v = point1()-point() ;
  Vector_3 e2v = point2()-point() ;
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

  RT scale = determinant( orthohx, orthohy, orthohz,
                                     e1phx,   e1phy,   e1phz,
                                     e2phx,   e2phy,   e2phz );

  Aff_transformation_3
     point_to_origin(TRANSLATION,  - ( point() - ORIGIN ) );
  Aff_transformation_3
     rotate_and_more( t11,    t12,   t13,   RT0,
                      t21,    t22,   t23,   RT0,
                      t31,    t32,   t33,   RT0,
                                            scale);

  Point_3 ortho( orthohx, orthohy, orthohz );
  Point_3 e1p( e1phx, e1phy, e1phz );
  Point_3 e2p( e2phx, e2phy, e2phz );
  CGAL_kernel_assertion((   ortho.transform(rotate_and_more)
        == Point_3( RT(0), RT(0), RT(1)) ));
  CGAL_kernel_assertion((   e1p.transform(rotate_and_more)
        == Point_3( RT(1), RT(0), RT(0)) ));
  CGAL_kernel_assertion((   e2p.transform(rotate_and_more)
        == Point_3( RT(0), RT(1), RT(0)) ));

  return  rotate_and_more * point_to_origin;
}

template < class R >
CGAL_KERNEL_INLINE
typename PlaneH3<R>::Point_2
PlaneH3<R>::to_2d(const typename PlaneH3<R>::Point_3& p) const
{
  Point_3 tp = p.transform( transform_to_2d() );
  return Point_2( tp.hx(), tp.hy(), tp.hw());
}


template < class R >
CGAL_KERNEL_INLINE
typename PlaneH3<R>::Point_3
PlaneH3<R>::to_3d(const typename PlaneH3<R>::Point_2& p)  const
{
  Point_3 hp( p.hx(), p.hy(), RT(0.0), p.hw());
  return hp.transform( transform_to_2d().inverse() );
}

} //namespace CGAL

#endif  // CGAL_PLANEH3_H
