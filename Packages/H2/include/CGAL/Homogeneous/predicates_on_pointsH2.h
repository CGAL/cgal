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
 

#ifndef CGAL_PREDICATES_ON_POINTSH2_H
#define CGAL_PREDICATES_ON_POINTSH2_H

#include <CGAL/Homogeneous/PointH2.h>

CGAL_BEGIN_NAMESPACE

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
      return SMALLER;
  else
      return ( qV < pV ) ? LARGER : EQUAL;
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
      return false;
  else if ( pV < qV )
      return true;
  pV = phx * qhw;
  qV = qhx * phw;
  return ( pV <= qV );
}

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
Oriented_side
_where_wrt_L_wedge( const PointH2<R>& p, const PointH2<R>& q )
{
  Sign xs = CGAL_NTS sign( q.hx()*p.hw() - p.hx()*q.hw() );  // sign( qx - px )
  Sign ys = CGAL_NTS sign( q.hy()*p.hw() - p.hy()*q.hw() );  // sign( qy - py )

  if (( xs == NEGATIVE ) || ( ys == NEGATIVE ))
      return ON_NEGATIVE_SIDE;
  if (( xs == POSITIVE ) && ( ys == POSITIVE ))
      return ON_POSITIVE_SIDE;
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

CGAL_END_NAMESPACE

#endif // CGAL_PREDICATES_ON_POINTSH2_H
