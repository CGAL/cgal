// Copyright (c) 1999  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_REGULAR_TRIANGULATION_FTC3_H
#define CGAL_REGULAR_TRIANGULATION_FTC3_H

#include <CGAL/number_utils.h>
#include <CGAL/predicates/sign_of_determinant.h>

// This file contains the low level cartesian predicates
// used by the 3D regular triangulation.

namespace CGAL {

// return minus the sign of the 5x5 determinant [P,Q,R,S,T]
// where column [P] = transpose[px,py,pz,p^2 -wp,1]
template <class FT>
Oriented_side
power_testC3( const FT &px, const FT &py, const FT &pz, const FT &pwt,
              const FT &qx, const FT &qy, const FT &qz, const FT &qwt,
              const FT &rx, const FT &ry, const FT &rz, const FT &rwt,
              const FT &sx, const FT &sy, const FT &sz, const FT &swt,
              const FT &tx, const FT &ty, const FT &tz, const FT &twt)
{
    // We translate the points so that T becomes the origin.
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = pz - tz;
    FT dpt = CGAL_NTS square(dpx) + CGAL_NTS square(dpy) +
             CGAL_NTS square(dpz) + (twt - pwt);
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = qz - tz;
    FT dqt = CGAL_NTS square(dqx) + CGAL_NTS square(dqy) +
             CGAL_NTS square(dqz) + (twt - qwt);
    FT drx = rx - tx;
    FT dry = ry - ty;
    FT drz = rz - tz;
    FT drt = CGAL_NTS square(drx) + CGAL_NTS square(dry) +
             CGAL_NTS square(drz) + (twt - rwt);
    FT dsx = sx - tx;
    FT dsy = sy - ty;
    FT dsz = sz - tz;
    FT dst = CGAL_NTS square(dsx) + CGAL_NTS square(dsy) +
             CGAL_NTS square(dsz) + (twt - swt);

    return - sign_of_determinant(dpx, dpy, dpz, dpt,
				 dqx, dqy, dqz, dqt,
				 drx, dry, drz, drt,
				 dsx, dsy, dsz, dst);
}


template <class FT>
Oriented_side
power_testC3( const FT &px, const FT &py, const FT &pz, const FT &pwt,
	      const FT &qx, const FT &qy, const FT &qz, const FT &qwt,
	      const FT &rx, const FT &ry, const FT &rz, const FT &rwt,
	      const FT &tx, const FT &ty, const FT &tz, const FT &twt)
{
    // Same translation as above.
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = pz - tz;
    FT dpt = CGAL_NTS square(dpx) + CGAL_NTS square(dpy) +
             CGAL_NTS square(dpz) + (twt - pwt);
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = qz - tz;
    FT dqt = CGAL_NTS square(dqx) + CGAL_NTS square(dqy) +
             CGAL_NTS square(dqz) + (twt - qwt);
    FT drx = rx - tx;
    FT dry = ry - ty;
    FT drz = rz - tz;
    FT drt = CGAL_NTS square(drx) + CGAL_NTS square(dry) +
             CGAL_NTS square(drz) + (twt - rwt);
    Sign cmp;

    // Projection on the (xy) plane.
    cmp = sign_of_determinant(dpx, dpy, dpt,
		                 dqx, dqy, dqt,
				 drx, dry, drt);
    if (cmp != ZERO)
	return cmp * sign_of_determinant(px-rx, py-ry,
		                         qx-rx, qy-ry);

    // Projection on the (xz) plane.
    cmp = sign_of_determinant(dpx, dpz, dpt,
		                 dqx, dqz, dqt,
				 drx, drz, drt);
    if (cmp != ZERO)
	return cmp * sign_of_determinant(px-rx, pz-rz,
		                         qx-rx, qz-rz);

    // Projection on the (yz) plane.
    cmp = sign_of_determinant(dpy, dpz, dpt,
		                 dqy, dqz, dqt,
				 dry, drz, drt);
    return cmp * sign_of_determinant(py-ry, pz-rz,
		                     qy-ry, qz-rz);
}


template <class FT>
Oriented_side
power_testC3( const FT &px, const FT &py, const FT &pz, const FT &pwt,
	      const FT &qx, const FT &qy, const FT &qz, const FT &qwt,
	      const FT &tx, const FT &ty, const FT &tz, const FT &twt)
{
    // Same translation as above.
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = pz - tz;
    FT dpt = CGAL_NTS square(dpx) + CGAL_NTS square(dpy) +
             CGAL_NTS square(dpz) + (twt - pwt);
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = qz - tz;
    FT dqt = CGAL_NTS square(dqx) + CGAL_NTS square(dqy) +
             CGAL_NTS square (dqz) + (twt - qwt);
    Comparison_result cmp;

    // We do an orthogonal projection on the (x) axis, if possible.
    cmp = CGAL_NTS compare(px, qx);
    if (cmp != EQUAL)
        return cmp * sign_of_determinant(dpx, dpt, dqx, dqt);

    // We do an orthogonal projection on the (y) axis, if possible.
    cmp = CGAL_NTS compare(py, qy);
    if (cmp != EQUAL)
        return cmp * sign_of_determinant(dpy, dpt, dqy, dqt);

    // We do an orthogonal projection on the (z) axis.
    cmp = CGAL_NTS compare(pz, qz);
    return cmp * sign_of_determinant(dpz, dpt, dqz, dqt);
}

template <class FT>
Oriented_side
power_testC3(const FT &pwt, const FT &qwt)
{
    return CGAL_NTS compare(qwt, pwt);
}

} //namespace CGAL

#endif // CGAL_REGULAR_TRIANGULATION_FTC3_H
