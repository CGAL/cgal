// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
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
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Michael Kerber   <mkerber@mpi-inf.mpg.de>
//                 Dominik Huelse   <dominik.huelse@gmx.de>
//                 Michael Hemmer   <hemmer@informatik.uni-mainz.de>
//                 Eric Berberich   <eric.berberich@cgal.org>
// ============================================================================

/*! \file CGAL/Polynomial/polynomial_gcd_ntl.h
 *   \brief special polynomial gcd function via NTL  
 */

#ifndef CGAL_POLYNOMIAL_GCD_NTL_H
#define CGAL_POLYNOMIAL_GCD_NTL_H

#include <CGAL/config.h>

#ifndef CGAL_USE_NTL
#warning This header file needs NTL installed in order to work properly.
#endif

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#endif

#ifdef CGAL_USE_CORE
#include <CGAL/CORE_BigInt.h>
#endif

#include <CGAL/Polynomial.h>
#include <CGAL/polynomial_utils.h>
#include <CGAL/Polynomial_traits_d.h>

#include <sstream>

#include <NTL/ZZX.h>

namespace CGAL{
template <class A> class Polynomial; // fwd 

template <typename Polynomial_d> class Polynomial_traits_d;

} // namespace CGAL


// This part forms the bridge to NTL to use the modular gcd algorithm. If 
// NTL is not available, the usual strategy is applied. 

namespace CGAL {
namespace internal {

// Forward
template <class NT> 
Polynomial<NT> gcd_utcf(
        const Polynomial<NT>& FF1 ,
        const Polynomial<NT>& FF2 );
 
template<typename PolyInt>
inline
void polynomial_to_ntl(const PolyInt& p, NTL::ZZX& q) {
    std::stringstream ss;
    ss << "[ ";
    for(int i=0;i<=p.degree();i++) {
        ss << p[i] << " ";
    }
    ss << "]";
    ss >> q;
}
      
template<typename PolyInt>
inline
void ntl_to_polynomial(const NTL::ZZX& q,PolyInt& p) {
    int d = NTL::deg(q);
    if(d==-1) {
        p=PolyInt(1);
        return;
    }
    std::stringstream ss;
    ss << "P[";
    ss << d;
    for(int i=0;i<=d;i++) {
        ss << "(" << i << "," << NTL::coeff(q,i) << ")";
    }
    ss << "]";
    p=PolyInt::input_ascii(ss);
}

template<typename NT> Polynomial<NT>
inline
modular_NTL_gcd_for_univariate_integer_polynomials
  (Polynomial<NT> p1, Polynomial<NT> p2) {     
  //    std::cout<<" NTL GCD"<<std::endl;
    
    NTL::ZZX q1,q2,h;   
    Polynomial<NT> g;
    internal::polynomial_to_ntl(p1,q1);
    internal::polynomial_to_ntl(p2,q2);
#ifdef CGAL_MODULAR_GCD_TIMER
    timer_ntl2.start(); 
#endif
    NTL::GCD(h,q1,q2);
#ifdef CGAL_MODULAR_GCD_TIMER
    timer_ntl2.stop(); 
#endif
    internal::ntl_to_polynomial(h,g);
    return g;
}

template<typename NT> Polynomial<NT>
inline
canonical_modular_NTL_gcd_for_univariate_integer_polynomials
  (Polynomial<NT> p1, Polynomial<NT> p2) {     
  //    std::cout<<" NTL canonical GCD"<<std::endl;
  return CGAL::canonicalize(modular_NTL_gcd_for_univariate_integer_polynomials(p1,p2));
}


#ifdef CGAL_USE_LEDA

template <> 
inline
CGAL::Polynomial<leda::integer>
gcd_utcf_(const CGAL::Polynomial<leda::integer>& p1,
          const CGAL::Polynomial<leda::integer>& p2) {
    CGAL::Polynomial<leda::integer> gcd = 
        internal::canonical_modular_NTL_gcd_for_univariate_integer_polynomials(p1,p2);
    return gcd;
}

template <> 
inline
CGAL::Polynomial<leda::integer>
gcd_(const CGAL::Polynomial<leda::integer>& p1,
     const CGAL::Polynomial<leda::integer>& p2) {
    return internal::modular_NTL_gcd_for_univariate_integer_polynomials(p1,p2);
}

#endif // CGAL_USE_LEDA

#ifdef CGAL_USE_CORE

template <> 
inline
Polynomial<CORE::BigInt>
gcd_utcf_(const Polynomial<CORE::BigInt>& p1,
          const Polynomial<CORE::BigInt>& p2) {
    Polynomial<CORE::BigInt> gcd = canonical_modular_NTL_gcd_for_univariate_integer_polynomials(p1,p2);
    return gcd;
}

template <> 
inline
Polynomial<CORE::BigInt>
gcd_(const Polynomial<CORE::BigInt>& p1,
    const Polynomial<CORE::BigInt>& p2) {
    return modular_NTL_gcd_for_univariate_integer_polynomials(p1,p2);
}

#endif //CGAL_USE_CORE

  } // namespace internal

} // namespace CGAL

#endif // CGAL_POLYNOMIAL_GCD_NTL_H

// EOF
