// Copyright (c) 2002-2008 Max-Planck-Institute Saarbruecken (Germany)
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
 
