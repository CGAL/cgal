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


#ifndef CGAL_PREDICATES_ON_POINTSH2_H
#define CGAL_PREDICATES_ON_POINTSH2_H

#include <CGAL/Homogeneous/PointH2.h>

namespace CGAL {

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

template <class RT>
Comparison_result
compare_power_distanceH2(const RT& phx, const RT& phy, const RT& phw,
                         const Quotient<RT>& pwt,
                         const RT& qhx, const RT& qhy, const RT& qhw,
                         const Quotient<RT>& qwt,
                         const RT& rhx, const RT& rhy, const RT& rhw)
{
  // returns SMALLER if r is closer to p w.r.t. the power metric
  RT dphx = rhx * phw - phx * rhw;
  RT dphy = rhy * phw - phy * rhw;
  RT dqhx = rhx * qhw - qhx * rhw;
  RT dqhy = rhy * qhw - qhy * rhw;
  RT dphw = CGAL_NTS square(phw);
  RT dqhw = CGAL_NTS square(qhw);
  RT drhw = CGAL_NTS square(rhw);

  RT npwt = pwt.numerator();
  RT dpwt = pwt.denominator();
  RT nqwt = qwt.numerator();
  RT dqwt = qwt.denominator();


  RT dh1 = (CGAL_NTS square(dphx) + CGAL_NTS square(dphy))*dpwt - npwt * dphw * drhw;
  RT dh2 = (CGAL_NTS square(dqhx) + CGAL_NTS square(dqhy))*dqwt - nqwt * dqhw * drhw;
  return CGAL_NTS compare(dh1 * dqhw * dqwt, dh2 * dphw * dpwt );
}


template <class RT>
Oriented_side
power_testH2( const RT &phx, const RT &phy, const RT &phw, const Quotient<RT> &pwt,
              const RT &qhx, const RT &qhy, const RT &qhw, const Quotient<RT> &qwt,
              const RT &rhx, const RT &rhy, const RT &rhw, const Quotient<RT> &rwt,
              const RT &thx, const RT &thy, const RT &thw, const Quotient<RT> &twt)
{
    RT npwt = pwt.numerator();
    RT dpwt = pwt.denominator();
    RT nqwt = qwt.numerator();
    RT dqwt = qwt.denominator();
    RT nrwt = rwt.numerator();
    RT drwt = rwt.denominator();
    RT ntwt = twt.numerator();
    RT dtwt = twt.denominator();

    RT dphx = phx*phw;
    RT dphy = phy*phw;
    RT dphw = CGAL_NTS square(phw);
    RT dpz = (CGAL_NTS square(phx) + CGAL_NTS square(phy))*dpwt  - npwt*dphw;

    RT dqhx = qhx*qhw;
    RT dqhy = qhy*qhw;
    RT dqhw = CGAL_NTS square(qhw);
    RT dqz = (CGAL_NTS square(qhx) + CGAL_NTS square(qhy))*dqwt - nqwt*dqhw;

    RT drhx = rhx*rhw;
    RT drhy = rhy*rhw;
    RT drhw = CGAL_NTS square(rhw);
    RT drz = (CGAL_NTS square(rhx) + CGAL_NTS square(rhy)) *drwt - nrwt*drhw;

    RT dthx = thx*thw;
    RT dthy = thy*thw;
    RT dthw = CGAL_NTS square(thw);
    RT dtz = (CGAL_NTS square(thx) + CGAL_NTS square(thy))*dtwt - ntwt*dthw;

    dthx *= dtwt;
    dthy *= dtwt;
    dthw *= dtwt;

    return sign_of_determinant(dphx, dphy, dpz, dphw,
                               dqhx, dqhy, dqz, dqhw,
                               drhx, drhy, drz, drhw,
                               dthx, dthy, dtz, dthw);
}


template <class RT>
Oriented_side
power_testH2( const RT &phx, const RT &phy, const RT &phw, const Quotient<RT> &pwt,
              const RT &qhx, const RT &qhy, const RT &qhw, const Quotient<RT> &qwt,
              const RT &thx, const RT &thy, const RT &thw, const Quotient<RT> &twt)
{
    RT npwt = pwt.numerator();
    RT dpwt = pwt.denominator();
    RT nqwt = qwt.numerator();
    RT dqwt = qwt.denominator();
    RT ntwt = twt.numerator();
    RT dtwt = twt.denominator();

    // Test if we can project on the (x) axis.  If not, then on the
    // (y) axis
    RT pa, qa, ta;

    if (phx * qhw != qhx * phw )
    {
        pa = phx*phw;
        qa = qhx*qhw;
        ta = thx*thw;
    }
    else
    {
        pa = phy*phw;
        qa = qhy*qhw;
        ta = thy*thw;
    }

    RT dphw = CGAL_NTS square(phw);
    RT dpz = (CGAL_NTS square(phx) + CGAL_NTS square(phy))*dpwt - npwt*dphw;

    RT dqhw = CGAL_NTS square(qhw);
    RT dqz = (CGAL_NTS square(qhx) + CGAL_NTS square(qhy))*dqwt - nqwt*dqhw;

    RT dthw = CGAL_NTS square(thw);
    RT dtz = (CGAL_NTS square(thx) + CGAL_NTS square(thy))*dtwt - ntwt*dthw;

    pa *= dtwt;
    qa *= dtwt;
    ta *= dtwt;

    return CGAL_NTS compare(pa, qa) * sign_of_determinant(pa, dpz, dphw,
                                                          qa, dqz, dqhw,
                                                          ta, dtz, dthw);
}

#if 0
// Unused, undocumented, un-functorized.
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
#endif

} //namespace CGAL

#endif // CGAL_PREDICATES_ON_POINTSH2_H
