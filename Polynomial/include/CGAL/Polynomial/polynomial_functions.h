// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : CGAL
// File          : include/CGAL/Polynomial/polynomial_functions.h
// CGAL_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Arno Eigenwillig <arno@mpi-inf.mpg.de>
//
// ============================================================================

/*! \file CGAL/Polynomial/polynomial_functions.h
 *   \brief Various small helper functions for polynomials
 */

#ifndef CGAL_POLYNOMIAL_FUNCTIONS_H
#define CGAL_POLYNOMIAL_FUNCTIONS_H

#include <CGAL/config.h>
#include <CGAL/Polynomial.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class NT>
bool _check_leadcoeff(int x0dm2_degree, int xd_degree,
                      const CGAL::Polynomial< NT > &p) {
    return (x0dm2_degree + p.degree() <= xd_degree);
}


template <class NT>
bool _check_leadcoeff(int xdmj_inner_degree, int xd_degree,
                      const CGAL::Polynomial< CGAL::Polynomial< NT > > &p) {
    
    int d = p.degree();
    while (d >= 0) {
        if (!_check_leadcoeff(xdmj_inner_degree + d, xd_degree, p[d])) {
            return false;
        }
        d--;
    }
    return true;
}

} // namespace CGALi

/*!\brief 
 * check that \c p has a non-zero coefficient
 * of <I>x_d</I><SUP>deg(<I>p</I>)</SUP>
 */
template <class NT>
bool check_leadcoeff(const CGAL::Polynomial< CGAL::Polynomial< NT > > &p) {
    if (CGAL::is_zero(p)) {
        return true;
    }

    return CGALi::_check_leadcoeff(0, p.degree(), p);
}

CGAL_END_NAMESPACE

#endif // CGAL_POLYNOMIAL_FUNCTIONS_H
// EOF
