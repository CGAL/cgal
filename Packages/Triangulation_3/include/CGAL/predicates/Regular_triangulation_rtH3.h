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
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_REGULAR_TRIANGULATION_RTH3_H
#define CGAL_REGULAR_TRIANGULATION_RTH3_H

// This file contains the low level homogeneous predicates
// used by the 3D regular triangulation.

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
    RT dphw = square(phw);
    RT dpz = square(phx) + square(phy) + square(phz) - pwt*dphw;

    RT dqhx = qhx*qhw;
    RT dqhy = qhy*qhw;
    RT dqhz = qhz*qhw;
    RT dqhw = square(qhw);
    RT dqz = square(qhx) + square(qhy) + square(qhz) - qwt*dqhw;

    RT drhx = rhx*rhw;
    RT drhy = rhy*rhw;
    RT drhz = rhz*rhw;
    RT drhw = square(rhw);
    RT drz = square(rhx) + square(rhy) + square(rhz) - rwt*drhw;

    RT dshx = shx*shw;
    RT dshy = shy*shw;
    RT dshz = shz*shw;
    RT dshw = square(shw);
    RT dsz = square(shx) + square(shy) + square(shz) - swt*dshw;

    RT dthx = thx*thw;
    RT dthy = thy*thw;
    RT dthz = thz*thw;
    RT dthw = square(thw);
    RT dtz = square(thx) + square(thy) + square(thz) - twt*dthw;

    return Oriented_side(sign_of_determinant5x5(dphx, dphy, dphz, dpz, dphw,
	                                        dqhx, dqhy, dqhz, dqz, dqhw,
	                                        drhx, drhy, drhz, drz, drhw,
	                                        dshx, dshy, dshz, dsz, dshw,
	                                        dthx, dthy, dthz, dtz, dthw));
}

#error This one is not yet ported to 3D.
template <class RT>
Oriented_side
power_testH3( const RT &phx, const RT &phy, const RT &phw, const RT &pwt,
              const RT &qhx, const RT &qhy, const RT &qhw, const RT &qwt,
              const RT &thx, const RT &thy, const RT &thw, const RT &twt)
{
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

    RT dphw = square(phw);
    RT dpz = square(phx) + square(phy) - pwt*dphw;

    RT dqhw = square(qhw);
    RT dqz = square(qhx) + square(qhy) - qwt*dqhw;

    RT dthw = square(thw);
    RT dtz = square(thx) + square(thy) - twt*dthw;

    return Oriented_side(CGAL::compare(pa, qa) *
	                 sign_of_determinant3x3(pa, dpz, dphw,
				                qa, dqz, dqhw,
				                ta, dtz, dthw));
}

#error one more predicate needed in 3D...

CGAL_END_NAMESPACE

#ifdef CGAL_ARITHMETIC_FILTER_H
#ifndef CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_RTH3_H
#include <CGAL/Arithmetic_filter/predicates/Regular_triangulation_rtH3.h>
#endif // CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_RTH3_H
#endif

#endif // CGAL_REGULAR_TRIANGULATION_RTH3_H
