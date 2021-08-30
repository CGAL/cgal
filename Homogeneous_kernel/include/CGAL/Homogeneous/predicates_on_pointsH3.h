// Copyright (c) 1999,2016
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra, Olivier Devillers, Mariette Yvinec


#ifndef CGAL_PREDICATES_ON_POINTSH3_H
#define CGAL_PREDICATES_ON_POINTSH3_H

#include <CGAL/Homogeneous/PointH3.h>
#include <CGAL/predicates/kernel_ftC3.h>

namespace CGAL {

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



template <class RT>
Oriented_side
power_side_of_oriented_power_sphereH3(
             const RT &phx, const RT &phy, const RT &phz, const RT &phw, const Quotient<RT> &pwt,
             const RT &qhx, const RT &qhy, const RT &qhz, const RT &qhw, const Quotient<RT> &qwt,
             const RT &rhx, const RT &rhy, const RT &rhz, const RT &rhw, const Quotient<RT> &rwt,
             const RT &shx, const RT &shy, const RT &shz, const RT &shw, const Quotient<RT> &swt,
             const RT &thx, const RT &thy, const RT &thz, const RT &thw, const Quotient<RT> &twt)
{
    RT npwt = pwt.numerator();
    RT dpwt = pwt.denominator();
    RT nqwt = qwt.numerator();
    RT dqwt = qwt.denominator();
    RT nrwt = rwt.numerator();
    RT drwt = rwt.denominator();
    RT nswt = swt.numerator();
    RT dswt = swt.denominator();
    RT ntwt = twt.numerator();
    RT dtwt = twt.denominator();

    RT dphx = phx*phw;
    RT dphy = phy*phw;
    RT dphz = phz*phw;
    RT dphw = CGAL_NTS square(phw);
    RT dpz = (CGAL_NTS square(phx) + CGAL_NTS square(phy) +
              CGAL_NTS square(phz))*dpwt - npwt*dphw;

    RT dqhx = qhx*qhw;
    RT dqhy = qhy*qhw;
    RT dqhz = qhz*qhw;
    RT dqhw = CGAL_NTS square(qhw);
    RT dqz = (CGAL_NTS square(qhx) + CGAL_NTS square(qhy) +
              CGAL_NTS square(qhz))*dqwt - nqwt*dqhw;

    RT drhx = rhx*rhw;
    RT drhy = rhy*rhw;
    RT drhz = rhz*rhw;
    RT drhw = CGAL_NTS square(rhw);
    RT drz = (CGAL_NTS square(rhx) + CGAL_NTS square(rhy) +
              CGAL_NTS square(rhz))*drwt - nrwt*drhw;

    RT dshx = shx*shw;
    RT dshy = shy*shw;
    RT dshz = shz*shw;
    RT dshw = CGAL_NTS square(shw);
    RT dsz = (CGAL_NTS square(shx) + CGAL_NTS square(shy) +
              CGAL_NTS square(shz))*dswt - nswt*dshw;

    RT dthx = thx*thw;
    RT dthy = thy*thw;
    RT dthz = thz*thw;
    RT dthw = CGAL_NTS square(thw);
    RT dtz = (CGAL_NTS square(thx) + CGAL_NTS square(thy) +
              CGAL_NTS square(thz))*dtwt - ntwt*dthw;

    dthx *= dtwt;
    dthy *= dtwt;
    dthz *= dtwt;
    dthw *= dtwt;

    return - sign_of_determinant(dphx, dphy, dphz, dpz, dphw,
                                 dqhx, dqhy, dqhz, dqz, dqhw,
                                 drhx, drhy, drhz, drz, drhw,
                                 dshx, dshy, dshz, dsz, dshw,
                                 dthx, dthy, dthz, dtz, dthw);
}

// The 2 degenerate are not speed critical, and they are quite boring and error
// prone to write, so we use the Cartesian version, using FT.

} //namespace CGAL

#endif // CGAL_PREDICATES_ON_POINTSH3_H
