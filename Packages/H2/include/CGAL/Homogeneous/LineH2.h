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
 

#ifndef CGAL_LINEH2_H
#define CGAL_LINEH2_H

#include <CGAL/Threetuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class LineH2
  : public R_::template Handle<Threetuple<typename R_::RT> >::type
{
CGAL_VC7_BUG_PROTECTED
    typedef typename R_::FT                   FT;
    typedef typename R_::RT                   RT;
    typedef typename R_::Point_2              Point_2;
    typedef typename R_::Vector_2             Vector_2;
    typedef typename R_::Direction_2          Direction_2;
    typedef typename R_::Segment_2            Segment_2;
    typedef typename R_::Ray_2                Ray_2;
    typedef typename R_::Line_2               Line_2;
    typedef typename R_::Aff_transformation_2 Aff_transformation_2;

    typedef Threetuple<RT>                           rep;
    typedef typename R_::template Handle<rep>::type  base;

    const base& Base() const { return *this; }
    base& Base() { return *this; }

public:
    typedef R_                                    R;

    LineH2() {}
    LineH2(const Point_2& p, const Point_2& q);
    LineH2(const RT& a, const RT& b, const RT& c);
    LineH2(const Segment_2& s);
    LineH2(const Ray_2& r);
    LineH2(const Point_2& p, const Direction_2& d);
    LineH2(const Point_2& p, const Vector_2& v);

    bool           operator==(const LineH2<R>& l) const ;
    bool           operator!=(const LineH2<R>& l) const ;

    const RT &     a() const { return get(Base()).e0; }
    const RT &     b() const { return get(Base()).e1; }
    const RT &     c() const { return get(Base()).e2; }

    FT             x_at_y(FT y) const;
    FT             y_at_x(FT x) const;

    Line_2     perpendicular(const Point_2& p ) const;
    Line_2     opposite() const;
    Point_2    point() const;
    Point_2    point(int i) const;
    Point_2    projection(const Point_2& p) const;
    Direction_2 direction() const;
    Vector_2   to_vector() const;
    Oriented_side  oriented_side( const Point_2& p ) const;
    bool           has_on( const Point_2& p ) const;
    bool           has_on_boundary( const Point_2& p ) const;
    bool           has_on_positive_side( const Point_2& p ) const;
    bool           has_on_negative_side( const Point_2& p ) const;
    bool           is_horizontal() const;
    bool           is_vertical()   const;
    bool           is_degenerate() const;

    Line_2     transform(const Aff_transformation_2&) const;
};

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
LineH2<R>::LineH2(const typename LineH2<R>::Point_2& p,
	          const typename LineH2<R>::Point_2& q)
 : base(
  //  a() * X + b() * Y + c() * W() == 0
  //      |    X        Y       W     |
  //      |  p.hx()   p.hy()  p.hw()  |
  //      |  q.hx()   q.hy()  q.hw()  |

            p.hy()*q.hw() - p.hw()*q.hy(),
            p.hw()*q.hx() - p.hx()*q.hw(),
            p.hx()*q.hy() - p.hy()*q.hx())
{}

template < class R >
CGAL_KERNEL_INLINE
LineH2<R>::LineH2(const RT& a, const RT& b, const RT& c)
 : base(a, b, c)
{}

template < class R >
CGAL_KERNEL_INLINE
LineH2<R>::LineH2(const typename LineH2<R>::Segment_2& s)
{
  Point_2 p = s.start();
  Point_2 q = s.end();
  Base() = rep (
            p.hy()*q.hw() - p.hw()*q.hy(),
            p.hw()*q.hx() - p.hx()*q.hw(),
            p.hx()*q.hy() - p.hy()*q.hx());
}

template < class R >
CGAL_KERNEL_INLINE
LineH2<R>::LineH2(const typename LineH2<R>::Ray_2& r)
{
  Point_2 p = r.start();
  Point_2 q = r.second_point();
  Base() = rep (
            p.hy()*q.hw() - p.hw()*q.hy(),
            p.hw()*q.hx() - p.hx()*q.hw(),
            p.hx()*q.hy() - p.hy()*q.hx() );
}

template < class R >
CGAL_KERNEL_INLINE
LineH2<R>::LineH2(const typename LineH2<R>::Point_2& p,
		  const typename LineH2<R>::Vector_2& v)
{
  Point_2 q = p + v;
  Base() = rep (
            p.hy()*q.hw() - p.hw()*q.hy(),
            p.hw()*q.hx() - p.hx()*q.hw(),
            p.hx()*q.hy() - p.hy()*q.hx() );
}

template < class R >
CGAL_KERNEL_INLINE
LineH2<R>::LineH2(const typename LineH2<R>::Point_2& p,
		  const typename LineH2<R>::Direction_2& d)
{
  Point_2 q = p + d.to_vector();
  Base() = rep (
            p.hy()*q.hw() - p.hw()*q.hy(),
            p.hw()*q.hx() - p.hx()*q.hw(),
            p.hx()*q.hy() - p.hy()*q.hx() );
}

template < class R >
CGAL_KERNEL_INLINE
typename LineH2<R>::FT
LineH2<R>::x_at_y(FT y) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return (FT(-b())*y - FT(c()) )/FT(a());
}

template < class R >
CGAL_KERNEL_INLINE
typename LineH2<R>::FT
LineH2<R>::y_at_x(FT x) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return (FT(-a())*x - FT(c()) )/FT(b());
}

template < class R >
CGAL_KERNEL_INLINE
typename R::Line_2
LineH2<R>::perpendicular(const typename LineH2<R>::Point_2& p ) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return typename R::Line_2( -b()*p.hw(), a()*p.hw(), b()*p.hx() - a()*p.hy());
}

template < class R >
inline
typename R::Line_2
LineH2<R>::opposite() const
{ return typename R::Line_2( -a(), -b(), -c() ); }

template < class R >
CGAL_KERNEL_INLINE
typename LineH2<R>::Point_2
LineH2<R>::point() const
{
  CGAL_kernel_precondition( !is_degenerate() );
  if (is_vertical() )
  {
      return Point_2(-c(), RT(0)  , a() );
  }
  else
  {
      return Point_2(RT(0)  , -c(), b() );
  }
}

template < class R >
CGAL_KERNEL_INLINE
typename LineH2<R>::Point_2
LineH2<R>::point(int i) const
{ return point() + RT(i) * to_vector(); }

template < class R >
CGAL_KERNEL_INLINE
typename LineH2<R>::Point_2
LineH2<R>::projection(const typename LineH2<R>::Point_2& p) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  LineH2<R>  l( p, Direction_2( a(), b() ));
  return Point_2( b()*l.c() - l.b()*c(),
                  l.a()*c() - a()*l.c(),
                  a()*l.b() - l.a()*b() );
}

template < class R >
CGAL_KERNEL_INLINE
typename LineH2<R>::Vector_2
LineH2<R>::to_vector() const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return Vector_2( b(), -a() );
}

template < class R >
CGAL_KERNEL_INLINE
typename LineH2<R>::Direction_2
LineH2<R>::direction() const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return Direction_2( b(), -a() );
}

template < class R >
CGAL_KERNEL_INLINE
typename R::Line_2
LineH2<R>::transform(const typename LineH2<R>::Aff_transformation_2& t) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  Point_2 p = point() + to_vector();
  return typename R::Line_2( t.transform(point() ), t.transform(p) );
}

#ifndef CGAL_NO_OSTREAM_INSERT_LINEH2
template < class R >
std::ostream &
operator<<(std::ostream &os, const LineH2<R> &l)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << l.a() << ' ' << l.b() << ' ' << l.c();
    case IO::BINARY :
        write(os, l.a());
        write(os, l.b());
        write(os, l.c());
        return os;
    default:
       return os << "LineH2(" << l.a() << ", " << l.b() << ", " << l.c() <<')';
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_LINEH2

#ifndef CGAL_NO_ISTREAM_EXTRACT_LINEH2
template < class R >
std::istream &
operator>>(std::istream &is, LineH2<R> &p)
{
  typename R::RT a, b, c;
  switch(is.iword(IO::mode))
  {
    case IO::ASCII :
        is >> a >> b >> c;
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
  p = LineH2<R>(a, b, c);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_LINEH2

template < class R >
CGAL_KERNEL_INLINE
bool
LineH2<R>::has_on( const typename LineH2<R>::Point_2& p ) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return ( ( a()*p.hx() + b()*p.hy() + c()*p.hw() ) == RT(0)   );
}

template < class R >
CGAL_KERNEL_INLINE
bool
LineH2<R>::has_on_boundary( const typename LineH2<R>::Point_2& p ) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return ( ( a()*p.hx() + b()*p.hy() + c()*p.hw() ) == RT(0)   );
}

template < class R >
CGAL_KERNEL_INLINE
bool
LineH2<R>::has_on_positive_side( const typename LineH2<R>::Point_2& p ) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return ( ( a()*p.hx() + b()*p.hy() + c()*p.hw() ) > RT(0)   );
}

template < class R >
CGAL_KERNEL_INLINE
bool
LineH2<R>::has_on_negative_side( const typename LineH2<R>::Point_2& p ) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return ( ( a()*p.hx() + b()*p.hy() + c()*p.hw() ) < RT(0)   );
}

template < class R >
CGAL_KERNEL_INLINE
Oriented_side
LineH2<R>::oriented_side( const typename LineH2<R>::Point_2& p ) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  RT v = a()*p.hx() + b()*p.hy() + c()*p.hw();
  if (v > RT(0)   )
  {
      return ON_POSITIVE_SIDE;
  }
  else
  {
      return (v < RT(0)   ) ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
  }
}

template < class R >
inline
bool
LineH2<R>::is_horizontal() const
{ return ( a() == RT(0)   ); }

template < class R >
inline
bool
LineH2<R>::is_vertical() const
{ return ( b() == RT(0)   ); }

template < class R >
inline
bool
LineH2<R>::is_degenerate() const
{ return (a() == RT(0)  )&&(b() == RT(0)  ) ; }

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool
LineH2<R>::operator==(const LineH2<R>& l) const
{
  if (  (a() * l.c() != l.a() * c() )
      ||(b() * l.c() != l.b() * c() ) )
  {
      return false;
  }
  int sc  = static_cast<int>(CGAL_NTS sign(c()));
  int slc = static_cast<int>(CGAL_NTS sign(l.c()));
  if ( sc == slc )
  {
      if (sc == 0)
      {
          return (  (a()*l.b() == b()*l.a() )
                  &&(CGAL_NTS sign(a() )== CGAL_NTS sign( l.a() ))
                  &&(CGAL_NTS sign(b() )== CGAL_NTS sign( l.b() )) );
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
inline
bool
LineH2<R>::operator!=(const LineH2<R>& l) const
{ return !(*this == l); }

CGAL_END_NAMESPACE

#endif // CGAL_LINEH2_H
