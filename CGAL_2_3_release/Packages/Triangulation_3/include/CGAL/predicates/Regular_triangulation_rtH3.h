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
// file          : include/CGAL/predicates/Regular_triangulation_rtH3.h
// maintainer    : Monique.Teillaud@sophia.inria.fr
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_REGULAR_TRIANGULATION_RTH3_H
#define CGAL_REGULAR_TRIANGULATION_RTH3_H

// This file contains the low level homogeneous predicates
// used by the 3D regular triangulation.

#include <CGAL/predicates/Regular_triangulation_ftC3.h>

CGAL_BEGIN_NAMESPACE

template <class RT>
Oriented_side
power_testH3(
    const RT &phx, const RT &phy, const RT &phz, const RT &phw, const RT &pwt,
    const RT &qhx, const RT &qhy, const RT &qhz, const RT &qhw, const RT &qwt,
    const RT &rhx, const RT &rhy, const RT &rhz, const RT &rhw, const RT &rwt,
    const RT &shx, const RT &shy, const RT &shz, const RT &shw, const RT &swt,
    const RT &thx, const RT &thy, const RT &thz, const RT &thw, const RT &twt)
{
    RT dphx = phx*phw;
    RT dphy = phy*phw;
    RT dphz = phz*phw;
    RT dphw = CGAL_NTS square(phw);
    RT dpz = CGAL_NTS square(phx) + CGAL_NTS square(phy) + 
             CGAL_NTS square(phz) - pwt*dphw;

    RT dqhx = qhx*qhw;
    RT dqhy = qhy*qhw;
    RT dqhz = qhz*qhw;
    RT dqhw = CGAL_NTS square(qhw);
    RT dqz = CGAL_NTS square(qhx) + CGAL_NTS square(qhy) + 
             CGAL_NTS square(qhz) - qwt*dqhw;

    RT drhx = rhx*rhw;
    RT drhy = rhy*rhw;
    RT drhz = rhz*rhw;
    RT drhw = CGAL_NTS square(rhw);
    RT drz = CGAL_NTS square(rhx) + CGAL_NTS square(rhy) + 
             CGAL_NTS square(rhz) - rwt*drhw;

    RT dshx = shx*shw;
    RT dshy = shy*shw;
    RT dshz = shz*shw;
    RT dshw = CGAL_NTS square(shw);
    RT dsz = CGAL_NTS square(shx) + CGAL_NTS square(shy) + 
             CGAL_NTS square(shz) - swt*dshw;

    RT dthx = thx*thw;
    RT dthy = thy*thw;
    RT dthz = thz*thw;
    RT dthw = CGAL_NTS square(thw);
    RT dtz = CGAL_NTS square(thx) + CGAL_NTS square(thy) + 
             CGAL_NTS square(thz) - twt*dthw;

    return Oriented_side(- sign_of_determinant5x5(dphx, dphy, dphz, dpz, dphw,
	                                        dqhx, dqhy, dqhz, dqz, dqhw,
	                                        drhx, drhy, drhz, drz, drhw,
	                                        dshx, dshy, dshz, dsz, dshw,
	                                        dthx, dthy, dthz, dtz, dthw));
}

// The 2 degenerate are not speed critical, and they are quite boring and error
// prone to write, so we use the Cartesian version, using FT.

CGAL_END_NAMESPACE

#ifdef CGAL_ARITHMETIC_FILTER_H
#include <CGAL/Arithmetic_filter/predicates/Regular_triangulation_rtH3.h>
#endif

#endif // CGAL_REGULAR_TRIANGULATION_RTH3_H
