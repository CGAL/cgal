// Copyright (c) 1999  Utrecht University (The Netherlands),
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
// Author(s)     : Stefan Schirra
 

#ifndef CGAL_PREDICATES_ON_POINTSH2_H
#define CGAL_PREDICATES_ON_POINTSH2_H

#include <CGAL/Homogeneous/PointH2.h>

CGAL_BEGIN_NAMESPACE

template < class R>
CGAL_KERNEL_INLINE
bool
less_x(const PointH2<R>& p,
       const PointH2<R>& q)
{ return p.hx()*q.hw() < q.hx()*p.hw(); }

template < class R>
CGAL_KERNEL_INLINE
bool
less_y(const PointH2<R>& p,
       const PointH2<R>& q)
{ return p.hy()*q.hw() < q.hy()*p.hw(); }

template < class R>
CGAL_KERNEL_INLINE
bool
equal_xy(const PointH2<R>& p,
         const PointH2<R>& q)
{
  typedef typename R::RT RT;

  // Using these references allows to spare calls to [pq].hw().
  const RT& phw = p.hw();
  const RT& qhw = q.hw();

  return (p.hx()*qhw == q.hx()*phw) && (p.hy()*qhw == q.hy()*phw);
}

template < class R>
CGAL_KERNEL_INLINE
Comparison_result
compare_xy(const PointH2<R>& p, const PointH2<R>& q)
{
  typedef typename R::RT RT;

  const RT& phx = p.hx();
  const RT& phy = p.hy();
  const RT& phw = p.hw();
  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();

  RT pV = phx*qhw;
  RT qV = qhx*phw;
  if ( pV == qV )
  {
      pV = phy*qhw;
      qV = qhy*phw;
  }
  return CGAL_NTS compare(pV, qV);
}


template < class R>
CGAL_KERNEL_INLINE
bool
lexicographically_xy_smaller_or_equal(const PointH2<R>& p, const PointH2<R>& q)
{
  typedef typename R::RT RT;

  const RT& phx = p.hx();
  const RT& phy = p.hy();
  const RT& phw = p.hw();
  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();

  RT pV = phx * qhw;
  RT qV = qhx * phw;
  if ( qV < pV )
  {
      return false;
  }
  else if ( pV < qV )
  {
      return true;
  }

  pV = phy * qhw;
  qV = qhy * phw;
  return ( pV <= qV );
}

template < class R>
CGAL_KERNEL_INLINE
bool
lexicographically_xy_smaller(const PointH2<R>& p, const PointH2<R>& q)
{
  typedef typename R::RT RT;

  const RT& phx = p.hx();
  const RT& phy = p.hy();
  const RT& phw = p.hw();
  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();

  RT pV = phx * qhw;
  RT qV = qhx * phw;
  if ( qV < pV )
  {
      return false;
  }
  else if ( pV < qV )
  {
      return true;
  }
  pV = phy * qhw;
  qV = qhy * phw;
  return ( pV < qV );
}

template < class R>
CGAL_KERNEL_INLINE
bool
lexicographically_xy_larger_or_equal(const PointH2<R>& p, const PointH2<R>& q)
{
    return !lexicographically_xy_smaller(p,q);
}

template < class R>
CGAL_KERNEL_INLINE
bool
lexicographically_xy_larger(const PointH2<R>& p, const PointH2<R>& q)
{
  typedef typename R::RT RT;

  const RT& phx = p.hx();
  const RT& phy = p.hy();
  const RT& phw = p.hw();
  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();

  RT pV = phx * qhw;
  RT qV = qhx * phw;
  if ( pV < qV )
  {
      return false;
  }
  else if ( qV < pV )
  {
      return true;
  }
  pV = phy * qhw;
  qV = qhy * phw;
  return ( qV < pV  );
}

template < class R>
CGAL_KERNEL_INLINE
Comparison_result
compare_yx(const PointH2<R>& p, const PointH2<R>& q)
{
  typedef typename R::RT RT;

  const RT& phx = p.hx();
  const RT& phy = p.hy();
  const RT& phw = p.hw();
  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();

  RT pV = phy*qhw;
  RT qV = qhy*phw;
  if ( pV == qV )
  {
      pV = phx*qhw;
      qV = qhx*phw;
  }
  if ( pV < qV )
  {
      return SMALLER;
  }
  else
  {
      return ( qV < pV ) ? LARGER : EQUAL;
  }
}

template < class R>
CGAL_KERNEL_INLINE
bool
lexicographically_yx_smaller_or_equal(const PointH2<R>& p, const PointH2<R>& q)
{
  typedef typename R::RT RT;

  const RT& phx = p.hx();
  const RT& phy = p.hy();
  const RT& phw = p.hw();
  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();

  RT pV = phy * qhw;
  RT qV = qhy * phw;
  if ( qV < pV )
  {
      return false;
  }
  else if ( pV < qV )
  {
      return true;
  }
  pV = phx * qhw;
  qV = qhx * phw;
  return ( pV <= qV );
}

template < class R>
CGAL_KERNEL_INLINE
bool
lexicographically_yx_smaller(const PointH2<R>& p, const PointH2<R>& q)
{
  typedef typename R::RT RT;

  const RT& phx = p.hx();
  const RT& phy = p.hy();
  const RT& phw = p.hw();
  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();

  RT pV = phy * qhw;
  RT qV = qhy * phw;
  if ( qV < pV )
  {
      return false;
  }
  else if ( pV < qV )
  {
      return true;
  }
  pV = phx * qhw;
  qV = qhx * phw;
  return ( pV < qV );
}

template < class R>
CGAL_KERNEL_INLINE
bool
lexicographically_yx_larger_or_equal(const PointH2<R>& p, const PointH2<R>& q)
{
    return !lexicographically_yx_smaller(p, q);
}

template < class R>
CGAL_KERNEL_INLINE
bool
lexicographically_yx_larger(const PointH2<R>& p,
                                 const PointH2<R>& q)
{
  typedef typename R::RT RT;

  const RT& phx = p.hx();
  const RT& phy = p.hy();
  const RT& phw = p.hw();
  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();

  RT pV = phy * qhw;
  RT qV = qhy * phw;
  if ( pV < qV )
  {
      return false;
  }
  else if ( qV < pV )
  {
      return true;
  }
  pV = phx * qhw;
  qV = qhx * phw;
  return ( qV < pV );
}


template < class R>
CGAL_KERNEL_INLINE
Comparison_result
compare_x(const PointH2<R>& p, const PointH2<R>& q)
{
  return CGAL_NTS compare(p.hx() * q.hw(), q.hx() * p.hw());
}

template < class R>
CGAL_KERNEL_INLINE
Comparison_result
compare_y(const PointH2<R>& p, const PointH2<R>& q)
{
  return CGAL_NTS compare(p.hy() * q.hw(), q.hy() * p.hw());
}

template < class R>
CGAL_KERNEL_MEDIUM_INLINE
Orientation
orientation( const PointH2<R>& p, const PointH2<R>& q, const PointH2<R>& r)
{
  typedef typename R::RT RT;

  const RT& phx = p.hx();
  const RT& phy = p.hy();
  const RT& phw = p.hw();
  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();
  const RT& rhx = r.hx();
  const RT& rhy = r.hy();
  const RT& rhw = r.hw();
  const RT  RT0 = RT(0);

  // | A B |
  // | C D |

  RT  A = phx*rhw - phw*rhx;
  RT  B = phy*rhw - phw*rhy;
  RT  C = qhx*rhw - qhw*rhx;
  RT  D = qhy*rhw - qhw*rhy;

  RT  det =  A*D - B*C;

/*
  RT det_old =   p.hx() * (q.hy()*r.hw() - q.hw()*r.hy() )
               + p.hy() * (q.hw()*r.hx() - q.hx()*r.hw() )
               + p.hw() * (q.hx()*r.hy() - q.hy()*r.hx() );

  if ( !(CGAL_NTS sign(det) == CGAL_NTS sign(det_old)) )
  {
      std::cerr << "det: " << det << " det_old: " << det_old << flush;
  }
*/

  if (det < RT0  )
  {
      return CLOCKWISE;
  }
  else
  {
      return (RT0 < det) ? COUNTERCLOCKWISE : COLLINEAR;
  }
}

template < class R>
CGAL_KERNEL_MEDIUM_INLINE
Angle
angle( const PointH2<R>& p, const PointH2<R>& q, const PointH2<R>& r)
{
  return (Angle) CGAL_NTS sign((p-q)*(r-q));
}

template < class R>
CGAL_KERNEL_MEDIUM_INLINE
bool
left_turn( const PointH2<R>& p, const PointH2<R>& q, const PointH2<R>& r)
{
  typedef typename R::RT RT;

  const RT& phx = p.hx();
  const RT& phy = p.hy();
  const RT& phw = p.hw();
  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();
  const RT& rhx = r.hx();
  const RT& rhy = r.hy();
  const RT& rhw = r.hw();

  // | A B |
  // | C D |

  RT  A = phx*rhw - phw*rhx;
  RT  B = phy*rhw - phw*rhy;
  RT  C = qhx*rhw - qhw*rhx;
  RT  D = qhy*rhw - qhw*rhy;

  RT  det =  A*D - B*C;

/*
  RT det_old =   p.hx() * (q.hy()*r.hw() - q.hw()*r.hy() )
               + p.hy() * (q.hw()*r.hx() - q.hx()*r.hw() )
               + p.hw() * (q.hx()*r.hy() - q.hy()*r.hx() );

  if ( !(CGAL_NTS sign(det) == CGAL_NTS sign(det_old)) )
  {
      std::cerr << "det: " << det << " det_old: " << det_old << flush;
  }
*/


  return CGAL_NTS is_positive(det);
}

template < class R>
inline
bool
right_turn( const PointH2<R>& p, const PointH2<R>& q, const PointH2<R>& r)
{
  return left_turn(p, r, q);
}

template < class R>
CGAL_KERNEL_MEDIUM_INLINE
bool
collinear( const PointH2<R>& p, const PointH2<R>& q, const PointH2<R>& r)
{
  typedef typename R::RT RT;

  const RT& phx = p.hx();
  const RT& phy = p.hy();
  const RT& phw = p.hw();
  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();
  const RT& rhx = r.hx();
  const RT& rhy = r.hy();
  const RT& rhw = r.hw();

  // | A B |
  // | C D |

  RT  A = phx*rhw - phw*rhx;
  RT  B = phy*rhw - phw*rhy;
  RT  C = qhx*rhw - qhw*rhx;
  RT  D = qhy*rhw - qhw*rhy;

  RT  det =  A*D - B*C;

/*
  RT det_old =   p.hx() * (q.hy()*r.hw() - q.hw()*r.hy() )
               + p.hy() * (q.hw()*r.hx() - q.hx()*r.hw() )
               + p.hw() * (q.hx()*r.hy() - q.hy()*r.hx() );

  if ( !(CGAL_NTS sign(det) == CGAL_NTS sign(det_old)) )
  {
      std::cerr << "det: " << det << " det_old: " << det_old << flush;
  }
*/

  return CGAL_NTS is_zero(det);
}

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
side_of_bounded_circle( const PointH2<R>& p,
                        const PointH2<R>& q,
                        const PointH2<R>& t)
{
  typedef typename R::RT RT;

  const RT& phx = p.hx();
  const RT& phy = p.hy();
  const RT& phw = p.hw();
  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();
  const RT& thx = t.hx();
  const RT& thy = t.hy();
  const RT& thw = t.hw();

  return Bounded_side( CGAL_NTS compare((thx*phw-phx*thw)*(qhx*thw-thx*qhw),
	                                (thy*phw-phy*thw)*(thy*qhw-qhy*thw)) );
}

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
side_of_bounded_circle( const PointH2<R>& q,
                        const PointH2<R>& r,
                        const PointH2<R>& s,
                        const PointH2<R>& t)
{
  typedef typename R::RT RT;

  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();
  const RT& rhx = r.hx();
  const RT& rhy = r.hy();
  const RT& rhw = r.hw();
  const RT& shx = s.hx();
  const RT& shy = s.hy();
  const RT& shw = s.hw();
  const RT& thx = t.hx();
  const RT& thy = t.hy();
  const RT& thw = t.hw();
  const RT  RT0 = RT(0);

  CGAL_kernel_precondition( ! collinear(q,r,s) );

  // compute sign of      |qx  qy  qx^2+qy^2  1 |   | a b c d |
  //                      |      --  r  --      | = | e f g h |
  //     determinant      |      --  s  --      | = | i j k l |
  //                      |      --  t  --      |   | m n o p |
  //           where

  RT a = qhx*qhw;
  RT b = qhy*qhw;
  RT c = qhx*qhx + qhy*qhy;
  RT d = qhw*qhw;

  RT e = rhx*rhw;
  RT f = rhy*rhw;
  RT g = rhx*rhx + rhy*rhy;
  RT h = rhw*rhw;

  RT i = shx*shw;
  RT j = shy*shw;
  RT k = shx*shx + shy*shy;
  RT l = shw*shw;

  RT m = thx*thw;
  RT n = thy*thw;
  RT o = thx*thx + thy*thy;
  RT p = thw*thw;

  RT det =   a * ( f*(k*p - l*o) + j*(h*o - g*p) + n*(g*l - h*k) )
           - e * ( b*(k*p - l*o) + j*(d*o - c*p) + n*(c*l - d*k) )
           + i * ( b*(g*p - h*o) + f*(d*o - c*p) + n*(c*h - d*g) )
           - m * ( b*(g*l - h*k) + f*(d*k - c*l) + j*(c*h - d*g) );


  if ( det == RT0 )
  {
      return ON_BOUNDARY;
  }
  else
  {
      if (orientation(q,r,s) == CLOCKWISE)
      {
          det = -det;
      }
      return (RT0 < det ) ? ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
  }
}

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
Oriented_side
side_of_oriented_circle( const PointH2<R>& q,
                         const PointH2<R>& r,
                         const PointH2<R>& s,
                         const PointH2<R>& t)
{
  typedef typename R::RT RT;

  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();
  const RT& rhx = r.hx();
  const RT& rhy = r.hy();
  const RT& rhw = r.hw();
  const RT& shx = s.hx();
  const RT& shy = s.hy();
  const RT& shw = s.hw();
  const RT& thx = t.hx();
  const RT& thy = t.hy();
  const RT& thw = t.hw();
  const RT  RT0 = RT(0);

  CGAL_kernel_precondition( ! collinear(q,r,s) );

  // compute sign of      |qx  qy  qx^2+qy^2  1 |   | a b c d |
  //                      |      --  r  --      | = | e f g h |
  //     determinant      |      --  s  --      | = | i j k l |
  //                      |      --  t  --      |   | m n o p |
  //           where

  RT a = qhx*qhw;
  RT b = qhy*qhw;
  RT c = qhx*qhx + qhy*qhy;
  RT d = qhw*qhw;

  RT e = rhx*rhw;
  RT f = rhy*rhw;
  RT g = rhx*rhx + rhy*rhy;
  RT h = rhw*rhw;

  RT i = shx*shw;
  RT j = shy*shw;
  RT k = shx*shx + shy*shy;
  RT l = shw*shw;

  RT m = thx*thw;
  RT n = thy*thw;
  RT o = thx*thx + thy*thy;
  RT p = thw*thw;

  RT det =   a * ( f*(k*p - l*o) + j*(h*o - g*p) + n*(g*l - h*k) )
           - e * ( b*(k*p - l*o) + j*(d*o - c*p) + n*(c*l - d*k) )
           + i * ( b*(g*p - h*o) + f*(d*o - c*p) + n*(c*h - d*g) )
           - m * ( b*(g*l - h*k) + f*(d*k - c*l) + j*(c*h - d*g) );


  if ( det < RT0 )
  {
      return ON_NEGATIVE_SIDE;
  }
  else
  {
      return (RT0 < det ) ? ON_POSITIVE_SIDE : ON_ORIENTED_BOUNDARY;
  }
}
template <class R>
CGAL_KERNEL_MEDIUM_INLINE
bool
collinear_are_ordered_along_line( const PointH2<R>& p,
                                  const PointH2<R>& q,
                                  const PointH2<R>& r )
{
  typedef typename R::RT RT;

  const RT& phx = p.hx();
  const RT& phy = p.hy();
  const RT& phw = p.hw();
  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();
  const RT& rhx = r.hx();
  const RT& rhy = r.hy();
  const RT& rhw = r.hw();

  if ( !(phx * rhw == rhx * phw ) )          // non-vertical ?
  {
     return !( (  ( phx * qhw < qhx * phw)
                &&( rhx * qhw < qhx * rhw))
             ||(  ( qhx * phw < phx * qhw)
                &&( qhx * rhw < rhx * qhw)) );
  }
  else if ( !(phy * rhw == rhy * phw ) )
  {
     return !( (  ( phy * qhw < qhy * phw)
                &&( rhy * qhw < qhy * rhw))
             ||(  ( qhy * phw < phy * qhw)
                &&( qhy * rhw < rhy * qhw)) );
  }
  else
     return (( phx*qhw == qhx*phw) && ( phy*qhw == qhy*phw));
}


template <class R>
CGAL_KERNEL_INLINE
bool
are_ordered_along_line( const PointH2<R>& p,
                        const PointH2<R>& q,
                        const PointH2<R>& r )
{
  if ( collinear(p,q,r) )
  {
     return collinear_are_ordered_along_line(p,q,r);
  }
  else
  {
     return false;
  }
}

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
bool
collinear_are_strictly_ordered_along_line( const PointH2<R>& p,
                                           const PointH2<R>& q,
                                           const PointH2<R>& r )
{
  typedef typename R::RT RT;

  const RT& phx = p.hx();
  const RT& phy = p.hy();
  const RT& phw = p.hw();
  const RT& qhx = q.hx();
  const RT& qhy = q.hy();
  const RT& qhw = q.hw();
  const RT& rhx = r.hx();
  const RT& rhy = r.hy();
  const RT& rhw = r.hw();

  if ( !(phx * rhw == rhx * phw ) )
  {
     return (   ( phx * qhw < qhx * phw)
              &&( qhx * rhw < rhx * qhw))
          ||(   ( qhx * phw < phx * qhw)    // ( phx * qhw > qhx * phw)
              &&( rhx * qhw < qhx * rhw));  // ( qhx * rhw > rhx * qhw)
  }
  else
  {
     return (   ( phy * qhw < qhy * phw)
              &&( qhy * rhw < rhy * qhw))
          ||(   ( qhy * phw < phy * qhw)    // ( phy * qhw > qhy * phw)
              &&( rhy * qhw < qhy * rhw));  // ( qhy * rhw > rhy * qhw)
  }
}


template <class R>
CGAL_KERNEL_INLINE
bool
are_strictly_ordered_along_line( const PointH2<R>& p,
                                 const PointH2<R>& q,
                                 const PointH2<R>& r )
{
  return collinear(p, q, r) &&
         collinear_are_strictly_ordered_along_line(p, q, r);
}

template <class R>
inline
bool
x_equal( const PointH2<R>& p, const PointH2<R>& q )
{ return p.hx()*q.hw() == q.hx()*p.hw(); }

template <class R>
inline
bool
y_equal( const PointH2<R>& p, const PointH2<R>& q )
{ return p.hy()*q.hw() == q.hy()*p.hw(); }

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
Oriented_side
_where_wrt_L_wedge( const PointH2<R>& p, const PointH2<R>& q )
{
  Sign xs = CGAL_NTS sign( q.hx()*p.hw() - p.hx()*q.hw() );  // sign( qx - px )
  Sign ys = CGAL_NTS sign( q.hy()*p.hw() - p.hy()*q.hw() );  // sign( qy - py )

  if (( xs == NEGATIVE ) || ( ys == NEGATIVE ))
  {
      return ON_NEGATIVE_SIDE;
  }
  if (( xs == POSITIVE ) && ( ys == POSITIVE ))
  {
      return ON_POSITIVE_SIDE;
  }
  return ON_ORIENTED_BOUNDARY;
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_deltax_deltay(const PointH2<R>& p,
                      const PointH2<R>& q,
                      const PointH2<R>& r,
                      const PointH2<R>& s)
{
  return CGAL_NTS compare(
                  CGAL_NTS abs(p.hx()*q.hw() - q.hx()*p.hw()) * r.hw()*s.hw(),
                  CGAL_NTS abs(r.hy()*s.hw() - s.hy()*r.hw()) * p.hw()*q.hw());
}

#ifndef CGAL_NO_DEPRECATED_CODE
template < class R>
inline
Comparison_result
compare_lexicographically_xy(const PointH2<R>& p, const PointH2<R>& q)
{
  bool THIS_FUNCTION_IS_DEPRECATED;
  return compare_xy(p, q);
}

template < class R>
inline
Comparison_result
compare_lexicographically_yx(const PointH2<R>& p, const PointH2<R>& q)
{
  bool THIS_FUNCTION_IS_DEPRECATED;
  return compare_yx(p, q);
}

template < class R>
inline
bool
leftturn( const PointH2<R>& p, const PointH2<R>& q, const PointH2<R>& r)
{
  bool THIS_FUNCTION_IS_DEPRECATED;
  return left_turn(p, q, r);
}

template < class R>
inline
bool
rightturn( const PointH2<R>& p, const PointH2<R>& q, const PointH2<R>& r)
{
  bool THIS_FUNCTION_IS_DEPRECATED;
  return right_turn(p, q, r);
}
#endif

CGAL_END_NAMESPACE

#endif // CGAL_PREDICATES_ON_POINTSH2_H
