// Copyright (c) 1999  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

#ifndef CGAL_REGULAR_TRIANGULATION_FTC2_H
#define CGAL_REGULAR_TRIANGULATION_FTC2_H

// This file contains the low level cartesian predicates
// used by the 2D regular triangulation.

CGAL_BEGIN_NAMESPACE

template <class FT>
Comparison_result
compare_power_distanceC2(const FT& px, const FT& py, const FT& pwt,
			 const FT& qx, const FT& qy, const FT& qwt,
			 const FT& rx, const FT& ry)
{
  // returns SMALLER if r is closer to p w.r.t. the power metric
  FT d1 = CGAL_NTS square(rx - px) + CGAL_NTS square(ry - py) - pwt;
  FT d2 = CGAL_NTS square(rx - qx) + CGAL_NTS square(ry - qy) - qwt;
  return CGAL_NTS compare(d1, d2);
}


template <class FT>
Oriented_side
power_testC2( const FT &px, const FT &py, const FT &pwt,
              const FT &qx, const FT &qy, const FT &qwt,
              const FT &rx, const FT &ry, const FT &rwt,
              const FT &tx, const FT &ty, const FT &twt)
{
    // Note: maybe this can be further optimized like the usual in_circle() ?

    // We translate the 4 points so that T becomes the origin.
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = CGAL_NTS square(dpx) + CGAL_NTS square(dpy) - pwt + twt;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = CGAL_NTS square(dqx) + CGAL_NTS square(dqy) - qwt + twt;
    FT drx = rx - tx;
    FT dry = ry - ty;
    FT drz = CGAL_NTS square(drx) + CGAL_NTS square(dry) - rwt + twt;

    return Oriented_side(sign_of_determinant3x3(dpx, dpy, dpz,
                                                dqx, dqy, dqz,
                                                drx, dry, drz));
}


template <class FT>
Oriented_side
power_testC2( const FT &px, const FT &py, const FT &pwt,
	      const FT &qx, const FT &qy, const FT &qwt,
	      const FT &tx, const FT &ty, const FT &twt)
{
    // Same translation as above.
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = CGAL_NTS square(dpx) + CGAL_NTS square(dpy) - pwt + twt;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = CGAL_NTS square(dqx) + CGAL_NTS square(dqy) - qwt + twt;

    // We do an orthogonal projection on the (x) axis, if possible.
    Comparison_result cmpx = CGAL_NTS compare(px, qx);
    if (cmpx != EQUAL)
	return Oriented_side(cmpx * sign_of_determinant2x2(dpx, dpz, dqx, dqz));

    // If not possible, then on the (y) axis.
    Comparison_result cmpy = CGAL_NTS compare(py, qy);
    return Oriented_side(cmpy * sign_of_determinant2x2(dpy, dpz, dqy, dqz));
}

CGAL_END_NAMESPACE

#ifdef CGAL_ARITHMETIC_FILTER_H
#ifndef CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_FTC2_H
#include <CGAL/Arithmetic_filter/predicates/Regular_triangulation_ftC2.h>
#endif // CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_FTC2_H
#endif

#endif // CGAL_REGULAR_TRIANGULATION_FTC2_H
