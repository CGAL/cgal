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
 

#ifndef CGAL_PREDICATES_ON_POINTSH3_H
#define CGAL_PREDICATES_ON_POINTSH3_H

#include <CGAL/Homogeneous/PointH3.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
bool
x_equal(const PointH3<R> &p,
        const PointH3<R> &q)
{ return p.x()*q.hw() == q.x()*p.hw(); }

template < class R >
inline
bool
y_equal(const PointH3<R> &p,
        const PointH3<R> &q)
{ return p.y()*q.hw() == q.y()*p.hw(); }

template < class R >
inline
bool
z_equal(const PointH3<R> &p,
        const PointH3<R> &q)
{ return p.z()*q.hw() == q.z()*p.hw(); }

template < class R>
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_xyz(const PointH3<R>& p, const PointH3<R>& q)
{
  typedef typename R::RT RT;
  RT pV = p.hx()*q.hw();
  RT qV = q.hx()*p.hw();
  if ( pV < qV )
  {
      return SMALLER;
  }
  if ( qV < pV )    //   ( pV > qV )
  {
      return LARGER;
  }
  // same x
  pV = p.hy()*q.hw();
  qV = q.hy()*p.hw();
  if ( pV < qV )
  {
      return SMALLER;
  }
  if ( qV < pV )    //   ( pV > qV )
  {
      return LARGER;
  }
  // same x and y
  pV = p.hz()*q.hw();
  qV = q.hz()*p.hw();
  if ( pV < qV )
  {
      return SMALLER;
  }
  else
  {
      return (qV < pV) ? LARGER : EQUAL;
  }
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool
lexicographically_xyz_smaller_or_equal(const PointH3<R> &p,
                                            const PointH3<R> &q)
{
  typedef typename R::RT RT;
  RT pV = p.hx()*q.hw();
  RT qV = q.hx()*p.hw();
  if ( qV < pV )    //   ( pV > qV )
  {
      return false;
  }
  if ( !(pV == qV) )
  {
      return true;
  }
  // same x
  pV = p.hy()*q.hw();
  qV = q.hy()*p.hw();
  if ( qV < pV )    //   ( pV > qV )
  {
      return false;
  }
  if ( !(pV == qV) )
  {
      return true;
  }
  // same x and y
  pV = p.hz()*q.hw();
  qV = q.hz()*p.hw();
  if ( qV < pV )    //   ( pV > qV )
  {
      return false;
  }
  return true;
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool lexicographically_xyz_smaller(const PointH3<R> &p,
                                   const PointH3<R> &q)
{
  typedef typename R::RT RT;
  RT pV = p.hx()*q.hw();
  RT qV = q.hx()*p.hw();
  if ( pV < qV )
  {
      return true;
  }
  if ( qV < pV )
  {
      return false;
  }
  // same x
  pV = p.hy()*q.hw();
  qV = q.hy()*p.hw();
  if ( pV < qV )
  {
      return true;
  }
  if ( qV < pV )
  {
      return false;
  }
  // same x and y
  pV = p.hz()*q.hw();
  qV = q.hz()*p.hw();
  if ( pV < qV )
  {
      return true;
  }
  return false;
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool lexicographically_xy_smaller(const PointH3<R> &p,
                                  const PointH3<R> &q)
{
  typedef typename R::RT RT;
  RT pV = p.hx()*q.hw();
  RT qV = q.hx()*p.hw();
  if ( pV < qV )
  {
      return true;
  }
  if ( qV < pV )
  {
      return false;
  }
  // same x
  pV = p.hy()*q.hw();
  qV = q.hy()*p.hw();
  if ( pV < qV )
  {
      return true;
  }
  return false;
}

template < class R>
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_xy(const PointH3<R>& p, const PointH3<R>& q)
{
  typedef typename R::RT RT;
  RT pV = p.hx()*q.hw();
  RT qV = q.hx()*p.hw();
  if ( pV < qV )
  {
      return SMALLER;
  }
  if ( qV < pV )    //   ( pV > qV )
  {
      return LARGER;
  }
  // same x
  pV = p.hy()*q.hw();
  qV = q.hy()*p.hw();
  if ( pV < qV )
  {
      return SMALLER;
  }
  if ( qV < pV )    //   ( pV > qV )
  {
      return LARGER;
  }
  // same x and y
  return EQUAL;
}

template < class R >
inline
Comparison_result
compare_x(const PointH3<R> &p, const PointH3<R> &q)
{ return CGAL_NTS compare(p.hx() * q.hw(), q.hx() * p.hw() ); }

template < class R >
inline
Comparison_result
compare_y(const PointH3<R> &p, const PointH3<R> &q)
{ return CGAL_NTS compare(p.hy() * q.hw(), q.hy() * p.hw() ); }

template < class R >
inline
Comparison_result
compare_z(const PointH3<R> &p, const PointH3<R> &q)
{ return CGAL_NTS compare(p.hz() * q.hw(), q.hz() * p.hw() ); }

template < class R >
CGAL_KERNEL_INLINE
bool
equal_xy(const PointH3<R> &p, const PointH3<R> &q)
{
  return   (p.hx() * q.hw() == q.hx() * p.hw() )
        && (p.hy() * q.hw() == q.hy() * p.hw() );
}

template < class R >  // ???  ->   ==
CGAL_KERNEL_INLINE
bool
equal_xyz(const PointH3<R> &p, const PointH3<R> &q)
{
  return   (p.hx() * q.hw() == q.hx() * p.hw() )
        && (p.hy() * q.hw() == q.hy() * p.hw() )
        && (p.hz() * q.hw() == q.hz() * p.hw() );
}

template < class R >
CGAL_KERNEL_INLINE
bool
less_x(const PointH3<R> &p, const PointH3<R> &q)
{ return   (p.hx() * q.hw() < q.hx() * p.hw() ); }


template < class R >
CGAL_KERNEL_INLINE
bool
less_y(const PointH3<R> &p, const PointH3<R> &q)
{ return   (p.hy() * q.hw() < q.hy() * p.hw() ); }

template < class R >
CGAL_KERNEL_INLINE
bool
less_z(const PointH3<R> &p, const PointH3<R> &q)
{ return   (p.hz() * q.hw() < q.hz() * p.hw() ); }

template < class R >
CGAL_KERNEL_INLINE
bool
collinear_are_ordered_along_line(const PointH3<R> &p,
                                      const PointH3<R> &q,
                                      const PointH3<R> &r)
{
  typedef typename R::RT RT;
  CGAL_kernel_exactness_precondition( collinear(p, q, r) );
  const RT phx = p.hx();
  const RT phw = p.hw();
  const RT qhx = q.hx();
  const RT qhw = q.hw();
  const RT rhx = r.hx();
  const RT rhw = r.hw();

  const RT pqx = phx*qhw;
  const RT qpx = qhx*phw;
  const RT prx = phx*rhw;
  const RT qrx = qhx*rhw;
  const RT rqx = rhx*qhw;
  const RT rpx = rhx*phw;

  if ( prx != rpx )   // px != rx
  {
                //    (px <= qx)&&(qx <= rx) || (px >= qx)&&(qx >= rx)
                // !(((qx <  px)||(rx <  qx))&&((px <  qx)||(qx <  rx)))
      return ! (   ((qpx < pqx) || (rqx < qrx))
                && ((pqx < qpx) || (qrx < rqx))  );
  }

  const RT phy = p.hy();
  const RT qhy = q.hy();
  const RT rhy = r.hy();

  const RT pqy = phy*qhw;
  const RT qpy = qhy*phw;
  const RT pry = phy*rhw;
  const RT qry = qhy*rhw;
  const RT rqy = rhy*qhw;
  const RT rpy = rhy*phw;

  if ( pry != rpy )
  {
      return ! (   ((qpy < pqy) || (rqy < qry))
                && ((pqy < qpy) || (qry < rqy))  );
  }

  const RT phz = p.hz();
  const RT qhz = q.hz();
  const RT rhz = r.hz();

  const RT pqz = phz*qhw;
  const RT qpz = qhz*phw;
  const RT prz = phz*rhw;
  const RT qrz = qhz*rhw;
  const RT rqz = rhz*qhw;
  const RT rpz = rhz*phw;

  if ( prz != rpz )
  {
      return ! (   ((qpz < pqz) || (rqz < qrz))
                && ((pqz < qpz) || (qrz < rqz))  );
  }
  // p == r
  return  ((rqx == qrx) && (rqy == qry) && (rqz == qrz));
}

CGAL_END_NAMESPACE

#endif // CGAL_PREDICATES_ON_POINTSH3_H
