// Copyright (c) 2002-2008 Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Ralf Schindlbeck <rschindl@mpi-inf.mpg.de>

/*!\file include/CGAL/Polynomial/polynomial_functions.h
 * \brief Various small helper functions for polynomials
 */

#ifndef CGAL_POLYNOMIAL_FUNCTIONS_H
#define CGAL_POLYNOMIAL_FUNCTIONS_H

#include <CGAL/config.h>

#include <vector>

#include <CGAL/Polynomial.h>

namespace CGAL {

// leading coefficient
//////////////////////

namespace internal {

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

} // namespace internal

/*!\brief
 * check that \c p has a non-zero coefficient
 * of <I>x_d</I><SUP>deg(<I>p</I>)</SUP>
 */
template <class NT>
bool check_leadcoeff(const CGAL::Polynomial< CGAL::Polynomial< NT > > &p) {
    if (CGAL::is_zero(p)) {
        return true;
    }

    return internal::_check_leadcoeff(0, p.degree(), p);
}

// partial substitutions
////////////////////////

/*!\brief substitute innermost variable
 *
 *  The substituted number \e x may have a more general number type
 *  NTX than the coefficient type NT, provided there is an explicit
 *  conversion NTX(NT).
 */
template < class NT, class NTX >
typename CGAL::Polynomial< typename CGAL::Coercion_traits< NT, NTX >::Type >
substitute_x(CGAL::Polynomial< CGAL::Polynomial< NT > > p, const NTX& x) {

    typedef CGAL::Polynomial_traits_d<
    CGAL::Polynomial < CGAL::Polynomial< NT > > > PT_d;

    const int d = PT_d::d;

    CGAL_precondition(d >= 1);

    typedef typename CGAL::Coercion_traits< CGAL::Polynomial< NT >,NTX > CT;
    typedef typename CT::Type Coercion;
    typedef typename CGAL::Coercion_traits< NT, NTX > CTi;
    typedef typename CTi::Type Coercion_i;
    typedef CGAL::Polynomial_traits_d< Coercion > PT_d1;

    std::vector< Coercion > replacements;
    replacements.push_back(typename CT::Cast()(x));
    for (int i = 0; i < d-1; i++) {
        Coercion repl(Coercion_i(0), Coercion_i(1));
        typename PT_d1::Move move;
        repl = move(repl, d-2,i);
        replacements.push_back(repl);
    }

    typename PT_d::Substitute substitute;
    Coercion sub =
        substitute(p, replacements.begin(), replacements.end());

    //typename PT_d::Get_coefficient coeff;
    //Coercion ret = coeff(sub,0,d-1);

    return sub;
}

/*!\brief substitute two innermost variables
 *
 *  The substituted numbers \e x and \e y may have a more general number type
 *  NTX than the coefficient type NT, provided there is an explicit
 *  conversion NTX(NT).
 */
template <class NT, class NTX>
typename CGAL::Polynomial< typename CGAL::Coercion_traits< NT, NTX >::Type >
substitute_xy(
    const CGAL::Polynomial< CGAL::Polynomial< CGAL::Polynomial< NT > > >& p,
    const NTX& x, const NTX& y
) {

    typedef CGAL::Polynomial_traits_d<
    CGAL::Polynomial< CGAL::Polynomial < CGAL::Polynomial< NT > > > > PT_d;

    const int d = PT_d::d;

    CGAL_precondition(d >= 2);

    typedef typename CGAL::Coercion_traits< CGAL::Polynomial< NT >, NTX > CT;
    typedef typename CT::Type Coercion;
//    typedef typename CGAL::Coercion_traits< NT, NTX > CTi;
//    typedef typename CT::Type Coercion_i;
    typedef CGAL::Polynomial_traits_d < Coercion > PT_dc;

    std::vector< Coercion > replacements;
    replacements.push_back(typename CT::Cast()(x));
    replacements.push_back(typename CT::Cast()(y));
    for (int i = 0; i < d-2; i++) {
        Coercion repl(typename PT_dc::Coefficient_type(0),
                      typename PT_dc::Coefficient_type(1));
        typename PT_dc::Move move;
        repl = move(repl, d-3,i);
        replacements.push_back(repl);
    }

    typename PT_d::Substitute substitute;
    Coercion sub =
        substitute(p, replacements.begin(), replacements.end());

    return sub;
}

} //namespace CGAL

#endif // CGAL_POLYNOMIAL_FUNCTIONS_H
// EOF
