// Copyright (c) 2002-2008 Max-Planck-Institute Saarbruecken (Germany)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
//Author(s) : Michael Hemmer <mhemmer@uni-mainz.de>

/*! \file CGAL/Polynomial/modular_gcd.h
  provides gcd for Polynomials, based on Modular arithmetic.
*/


#ifndef CGAL_POLYNOMIAL_MODULAR_GCD_H
#define CGAL_POLYNOMIAL_MODULAR_GCD_H 1

#include <CGAL/basic.h>
#include <CGAL/Residue.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Scalar_factor_traits.h>
#include <CGAL/Chinese_remainder_traits.h>
#include <CGAL/Polynomial/modular_gcd_utcf_dfai.h>
#include <CGAL/Polynomial/modular_gcd_utcf_algorithm_M.h>

namespace CGAL {
namespace internal {

template <class NT>
Polynomial<NT> modular_gcd_utcf(
        const Polynomial<NT>& FF1 ,
        const Polynomial<NT>& FF2 , Integral_domain_tag){
    return modular_gcd_utcf_dfai(FF1, FF2);
}

template <class NT>
Polynomial<NT> modular_gcd_utcf(
        const Polynomial<NT>& FF1 ,
        const Polynomial<NT>& FF2 , Unique_factorization_domain_tag){
    return modular_gcd_utcf_algorithm_M(FF1, FF2);
}

}// namespace internal
}///namespace CGAL

#endif //#ifndef CGAL_POLYNOMIAL_MODULAR_GCD_H 1

