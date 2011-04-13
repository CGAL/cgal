// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/predicates/Regular_triangulation_ftC3.h
// maintainer    : Monique.Teillaud@sophia.inria.fr
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_REGULAR_TRIANGULATION_FTC3_H
#define CGAL_REGULAR_TRIANGULATION_FTC3_H

// This file contains the low level cartesian predicates
// used by the 3D regular triangulation.

CGAL_BEGIN_NAMESPACE

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
             CGAL_NTS square(dpz) - pwt + twt;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = qz - tz;
    FT dqt = CGAL_NTS square(dqx) + CGAL_NTS square(dqy) +
             CGAL_NTS square(dqz) - qwt + twt;
    FT drx = rx - tx;
    FT dry = ry - ty;
    FT drz = rz - tz;
    FT drt = CGAL_NTS square(drx) + CGAL_NTS square(dry) + 
             CGAL_NTS square(drz) - rwt + twt;
    FT dsx = sx - tx;
    FT dsy = sy - ty;
    FT dsz = sz - tz;
    FT dst = CGAL_NTS square(dsx) + CGAL_NTS square(dsy) + 
             CGAL_NTS square(dsz) - swt + twt;

    return Oriented_side( - sign_of_determinant4x4(dpx, dpy, dpz, dpt,
						   dqx, dqy, dqz, dqt,
						   drx, dry, drz, drt,
						   dsx, dsy, dsz, dst));
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
             CGAL_NTS square(dpz) - pwt + twt;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = qz - tz;
    FT dqt = CGAL_NTS square(dqx) + CGAL_NTS square(dqy) + 
             CGAL_NTS square(dqz) - qwt + twt;
    FT drx = rx - tx;
    FT dry = ry - ty;
    FT drz = rz - tz;
    FT drt = CGAL_NTS square(drx) + CGAL_NTS square(dry) + 
             CGAL_NTS square(drz) - rwt + twt;
    Sign cmp;

    // Projection on the (xy) plane.
    cmp = sign_of_determinant3x3(dpx, dpy, dpt,
		                 dqx, dqy, dqt,
				 drx, dry, drt);
    if (cmp != ZERO)
	return Oriented_side(cmp * sign_of_determinant2x2(px-rx, py-ry,
		                                          qx-rx, qy-ry));

    // Projection on the (xz) plane.
    cmp = sign_of_determinant3x3(dpx, dpz, dpt,
		                 dqx, dqz, dqt,
				 drx, drz, drt);
    if (cmp != ZERO)
	return Oriented_side(cmp * sign_of_determinant2x2(px-rx, pz-rz,
		                                          qx-rx, qz-rz));

    // Projection on the (yz) plane.
    cmp = sign_of_determinant3x3(dpy, dpz, dpt,
		                 dqy, dqz, dqt,
				 dry, drz, drt);
    return Oriented_side(cmp * sign_of_determinant2x2(py-ry, pz-rz,
		                                      qy-ry, qz-rz));
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
             CGAL_NTS square(dpz) - pwt + twt;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = qz - tz;
    FT dqt = CGAL_NTS square(dqx) + CGAL_NTS square(dqy) + 
             CGAL_NTS square (dqz) - qwt + twt;
    Comparison_result cmp;

    // We do an orthogonal projection on the (x) axis, if possible.
    cmp = CGAL_NTS compare(px, qx);
    if (cmp != EQUAL)
        return Oriented_side(cmp * sign_of_determinant2x2(dpx, dpt, dqx, dqt));

    // We do an orthogonal projection on the (y) axis, if possible.
    cmp = CGAL_NTS compare(py, qy);
    if (cmp != EQUAL)
        return Oriented_side(cmp * sign_of_determinant2x2(dpy, dpt, dqy, dqt));

    // We do an orthogonal projection on the (z) axis.
    cmp = CGAL_NTS compare(pz, qz);
    return Oriented_side(cmp * sign_of_determinant2x2(dpz, dpt, dqz, dqt));
}

CGAL_END_NAMESPACE

#ifdef CGAL_ARITHMETIC_FILTER_H
#ifndef CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_FTC3_H
#include <CGAL/Arithmetic_filter/predicates/Regular_triangulation_ftC3.h>
#endif // CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_FTC3_H
#endif

#endif // CGAL_REGULAR_TRIANGULATION_FTC3_H
