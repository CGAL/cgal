// Copyright (c) 1999  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

#ifndef CGAL_REGULAR_TRIANGULATION_RTH3_H
#define CGAL_REGULAR_TRIANGULATION_RTH3_H

// This file contains the low level homogeneous predicates
// used by the 3D regular triangulation.

#include <CGAL/predicates/kernel_ftC3.h>

namespace CGAL {

template <class RT>
Oriented_side
power_testH3(
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

#endif // CGAL_REGULAR_TRIANGULATION_RTH3_H
