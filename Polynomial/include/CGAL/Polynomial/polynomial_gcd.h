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
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Tobias Reithmann <treith@mpi-inf.mpg.de>
//                 Michael Hemmer   <hemmer@informatik.uni-mainz.de>
//                 Michael Kerber   <mkerber@mpi-inf.mpg.de>
//                 Dominik Huelse   <dominik.huelse@gmx.de>
// ============================================================================

/*! \file CGAL/Polynomial/polynomial_gcd.h
 *   \brief Greatest common divisors and related operations on polynomials. 
 */

#ifndef CGAL_POLYNOMIAL_GCD_H
#define CGAL_POLYNOMIAL_GCD_H

#include <CGAL/config.h>

#ifndef CGAL_USE_INTERNAL_MODULAR_GCD
#define CGAL_USE_INTERNAL_MODULAR_GCD 1
#endif
 
#include <CGAL/basic.h>
#include <CGAL/Residue.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Scalar_factor_traits.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Polynomial/Polynomial_type.h>
#include <CGAL/Polynomial/misc.h>
#include <CGAL/Polynomial/polynomial_gcd_implementations.h>
#include <CGAL/polynomial_utils.h>

#ifdef CGAL_USE_NTL
#include <CGAL/Polynomial/polynomial_gcd_ntl.h>
#endif

#if CGAL_USE_INTERNAL_MODULAR_GCD 
#include <CGAL/Polynomial/modular_gcd.h>
#endif


// 1) gcd (basic form without cofactors)
//     uses three level of dispatch on tag types:
//     a) if the algebra type of the innermost coefficient is a field,
//         ask for decomposability. for UFDs compute the gcd directly
//     b) if NT supports integralization, the gcd is computed on
//         integralized polynomials
//     c) over a field (unless integralized), use the Euclidean algorithm;
//         over a UFD, use the subresultant algorithm
//
// NOTICE: For better performance, especially in AlciX, there exist special
// modular implementations for the polynmials with coefficient type 
// leda::integer and the CORE::BigInt type which use, when the 
// NTL library is available.
// see CGAL/Polynomial/polynomial_gcd_ntl.h

namespace CGAL { 
namespace internal {

template <class NT> 
inline
Polynomial<NT> gcd_(
        const Polynomial<NT>& p1, 
        const Polynomial<NT>& p2, 
        Field_tag)
{ 
    return CGAL::internal::gcd_utcf_(p1,p2);
}

template <class NT> 
inline
Polynomial<NT> gcd_(
        const Polynomial<NT>& p1, 
        const Polynomial<NT>& p2, 
        Unique_factorization_domain_tag)
{ 
	typedef Polynomial<NT> POLY;
    typedef Polynomial_traits_d<POLY> PT;  
    typedef typename PT::Innermost_coefficient_type IC; 
        
    typename PT::Multivariate_content mcont; 
    IC mcont_p1 = mcont(p1);
    IC mcont_p2 = mcont(p2);
    
    typename CGAL::Coercion_traits<POLY,IC>::Cast ictp; 
    POLY p1_ = CGAL::integral_division(p1,ictp(mcont_p1));
    POLY p2_ = CGAL::integral_division(p2,ictp(mcont_p2));
            
    return CGAL::internal::gcd_utcf_(p1_, p2_) * ictp(CGAL::gcd(mcont_p1, mcont_p2)); 
}



// name gcd() forwarded to the internal::gcd_() dispatch function
/*! \ingroup CGAL_Polynomial
 *  \relates CGAL::Polynomial
 *  \brief return the greatest common divisor of \c p1 and \c p2
 *
 *  \pre Requires \c Innermost_coefficient_type to be a \c Field or a \c UFDomain.
 */
template <class NT> 
inline
Polynomial<NT> gcd_(const Polynomial<NT>& p1, const Polynomial<NT>& p2)
{ 
    typedef typename internal::Innermost_coefficient_type<Polynomial<NT> >::Type IC;
    typedef typename Algebraic_structure_traits<IC>::Algebraic_category Algebraic_category;

    // Filter for zero-polynomials
    if( p1 == Polynomial<NT>(0) )
    	return p2;
    if( p2 == Polynomial<NT>(0) )
    	return p1;
    return internal::gcd_(p1,p2,Algebraic_category());
}

} // namespace internal

// 2) gcd_utcf computation
//     (gcd up to scalar factors, for non-UFD non-field coefficients)
//     a) first try to decompose the coefficients
//     b) second dispatch depends on the algebra type of NT

namespace internal {

template <class NT> Polynomial<NT> inline
gcd_utcf_(const Polynomial<NT>& p1, const Polynomial<NT>& p2){
    typedef CGAL::Fraction_traits< Polynomial<NT> > FT;
    typedef typename FT::Is_fraction Is_fraction;
    return gcd_utcf_is_fraction_(p1, p2, Is_fraction());
}

// is fraction ? 
template <class NT> Polynomial<NT> inline
gcd_utcf_is_fraction_( 
        const Polynomial<NT>& p1, 
        const Polynomial<NT>& p2, 
        ::CGAL::Tag_true)
{
    typedef Polynomial<NT> POLY;
    typedef Polynomial_traits_d<POLY> PT;
    typedef Fraction_traits<POLY> FT;

    typename FT::Denominator_type dummy;
    typename FT::Numerator_type p1i, p2i; 
    
    typename FT::Decompose()(p1,p1i, dummy);
    typename FT::Decompose()(p2,p2i, dummy);

    typename Coercion_traits<POLY,typename FT::Numerator_type>::Cast cast;
    return typename PT::Canonicalize()(cast(internal::gcd_utcf_(p1i, p2i)));
}

template <class NT> Polynomial<NT> inline
gcd_utcf_is_fraction_( 
        const Polynomial<NT>& p1, 
        const Polynomial<NT>& p2, 
        ::CGAL::Tag_false)
{
    typedef Algebraic_structure_traits< Polynomial<NT> > NTT;
    typedef CGAL::Modular_traits<Polynomial<NT> > MT;
    
    return gcd_utcf_modularizable_algebra_(
            p1,p2,typename MT::Is_modularizable(),typename NTT::Algebraic_category());
}

// is type modularizable 
template <class NT> Polynomial<NT> inline
gcd_utcf_modularizable_algebra_( 
        const Polynomial<NT>& p1,  
        const Polynomial<NT>& p2, 
        ::CGAL::Tag_false, 
        Integral_domain_tag){
    return internal::gcd_utcf_Integral_domain(p1, p2);   
}
template <class NT> Polynomial<NT> inline
gcd_utcf_modularizable_algebra_( 
        const Polynomial<NT>& p1, 
        const Polynomial<NT>& p2, 
        ::CGAL::Tag_false, 
        Unique_factorization_domain_tag){
    return internal::gcd_utcf_UFD(p1, p2);   
}
template <class NT> Polynomial<NT> inline
gcd_utcf_modularizable_algebra_( 
        const Polynomial<NT>& p1, 
        const Polynomial<NT>& p2, 
        ::CGAL::Tag_false, 
        Euclidean_ring_tag){
    return internal::gcd_Euclidean_ring(p1, p2);
}

#if CGAL_USE_INTERNAL_MODULAR_GCD 
template <class NT> Polynomial<NT> inline
gcd_utcf_modularizable_algebra_( 
        const Polynomial<NT>& p1, 
        const Polynomial<NT>& p2, 
        ::CGAL::Tag_true, 
        Integral_domain_tag tag){
    return modular_gcd_utcf(p1, p2, tag);
}
template <class NT> Polynomial<NT> inline
gcd_utcf_modularizable_algebra_( 
        const Polynomial<NT>& p1, 
        const Polynomial<NT>& p2, 
        ::CGAL::Tag_true, 
        Unique_factorization_domain_tag tag){
    return modular_gcd_utcf(p1, p2, tag);
//    return modular_gcd_utcf_algorithm_M(p1, p2);
}
#else
template <class NT> Polynomial<NT> inline
gcd_utcf_modularizable_algebra_( 
        const Polynomial<NT>& p1, 
        const Polynomial<NT>& p2, 
        ::CGAL::Tag_true, 
        Integral_domain_tag){
    return internal::gcd_utcf_Integral_domain(p1, p2);
}
template <class NT> Polynomial<NT> inline
gcd_utcf_modularizable_algebra_( 
        const Polynomial<NT>& p1, 
        const Polynomial<NT>& p2, 
        ::CGAL::Tag_true, 
        Unique_factorization_domain_tag){
    return internal::gcd_utcf_UFD(p1, p2);   
}
#endif

template <class NT> Polynomial<NT> inline
gcd_utcf_modularizable_algebra_( 
        const Polynomial<NT>& p1, 
        const Polynomial<NT>& p2, 
        ::CGAL::Tag_true, 
        Euclidean_ring_tag){
    // No modular algorithm available
    return internal::gcd_Euclidean_ring(p1, p2);
}

template <class NT> Polynomial<NT> inline
gcd_utcf(const Polynomial<NT>& p1, const Polynomial<NT>& p2){
    return internal::gcd_utcf_(p1,p2);
}

} // namespace internal 


// 3) extended gcd computation (with cofactors)
//     with dispatch similar to gcd

namespace internal {

template <class NT> 
inline
Polynomial<NT> gcdex_(
        Polynomial<NT> x, Polynomial<NT> y,
        Polynomial<NT>& xf, Polynomial<NT>& yf,
        ::CGAL::Tag_false
) {
    typedef typename Algebraic_structure_traits<NT>::Algebraic_category Algebraic_category;
    return gcdex_(x, y, xf, yf, Algebraic_category());
}

template <class NT>
inline
Polynomial<NT> gcdex_(
        Polynomial<NT> x, Polynomial<NT> y,
        Polynomial<NT>& xf, Polynomial<NT>& yf,
        Field_tag
) {
    /* The extended Euclidean algorithm for univariate polynomials.
     * See [Cohen, 1993], algorithm 3.2.2
     */
    typedef Polynomial<NT> POLY;
    typename Algebraic_structure_traits<NT>::Integral_div idiv;

    // handle trivial cases
    if (x.is_zero()) {
        if (y.is_zero()) CGAL_error_msg("gcdex(0,0) is undefined");
        xf = NT(0); yf = idiv(NT(1), y.unit_part());
        return yf * y;
    }
    if (y.is_zero()) {
        yf = NT(0); xf = idiv(NT(1), x.unit_part());
        return xf * x;
    }
    bool swapped = x.degree() < y.degree();
    if (swapped) { POLY t = x; x = y; y = t; }

    // main loop
    POLY u = x, v = y, q, r, m11(1), m21(0), m21old;
    for (;;) {
        /* invariant: (i) There exist m12 and m22 such that
         *   u = m11*x + m12*y
         *   v = m21*x + m22*y
         * (ii) and we have
         *   gcd(u,v) == gcd(x,y)
         */

        // compute next element of remainder sequence
        POLY::euclidean_division(u, v, q, r);  //  u == qv + r
        if (r.is_zero()) break;

        // update u and v while preserving invariant
        u = v; v = r;
        /* Since r = u - qv, this preserves invariant (part ii)
         * and corresponds to the matrix assignment
         *   (u) = (0  1) (u)
         *   (v)   (1 -q) (v)
         */
        m21old = m21; m21 = m11 - q*m21; m11 = m21old;
        /* This simulates the matching matrix assignment
         *   (m11 m12) = (0  1) (m11 m12)
         *   (m21 m22)   (1 -q) (m21 m22)
         * which preserves the invariant (part i)
         */
        if (r.degree() == 0) break;
    }
    /* postcondition:  invariant holds  and  v divides u */

    // make gcd unit-normal
    m21 /= v.unit_part(); v /= v.unit_part();

    // obtain m22 such that  v == m21*x + m22*y
    POLY m22;
    POLY::euclidean_division(v - m21*x, y, m22, r);
    CGAL_assertion(r.is_zero());

    // check computation
    CGAL_assertion(v == m21*x + m22*y);

    // return results
    if (swapped) {
        xf = m22; yf = m21;
    } else {
        xf = m21; yf = m22;
    }
    return v;
}

template <class NT> 
inline
Polynomial<NT> gcdex_(
        Polynomial<NT> x, Polynomial<NT> y,
        Polynomial<NT>& xf, Polynomial<NT>& yf,
        ::CGAL::Tag_true
) {
    typedef Polynomial<NT> POLY;
    typedef typename CGAL::Fraction_traits<POLY>::Numerator_type INTPOLY;
    typedef typename CGAL::Fraction_traits<POLY>::Denominator_type DENOM;
    typedef typename INTPOLY::NT INTNT;
    
    typename CGAL::Fraction_traits<POLY>::Decompose decompose;
    typename CGAL::Fraction_traits<POLY>::Compose   compose;
    

    // rewrite  x as xi/xd  and  y as yi/yd  with integral polynomials xi, yi
    DENOM xd, yd;
    x.simplify_coefficients();
    y.simplify_coefficients();
    INTPOLY xi ,yi; 
    decompose(x,xi,xd);
    decompose(y,yi,yd);

    // compute the integral gcd with cofactors:
    // vi = gcd(xi, yi);  vfi*vi == xfi*xi + yfi*yi
    INTPOLY xfi, yfi; INTNT vfi;
    INTPOLY vi = pseudo_gcdex(xi, yi, xfi, yfi, vfi);

    // proceed to vfi*v == xfi*x + yfi*y  with v = gcd(x,y) (unit-normal)
    POLY v = compose(vi, vi.lcoeff());
    v.simplify_coefficients();
    CGAL_assertion(v.unit_part() == NT(1));
    vfi *= vi.lcoeff(); xfi *= xd; yfi *= yd;

    // compute xf, yf such that gcd(x,y) == v == xf*x + yf*y
    xf = compose(xfi, vfi);
    yf = compose(yfi, vfi);
    xf.simplify_coefficients();
    yf.simplify_coefficients();
    return v;
}

} // namespace internal

/*! \ingroup CGAL_Polynomial
 *  \relates CGAL::Polynomial
 *  \brief compute gcd with cofactors
 *
 *  This function computes the gcd of polynomials \c p1 and \c p2
 *  along with two other polynomials \c f1 and \c f2 such that
 *  gcd(\e p1, \e p2) = <I>f1*p1 + f2*p2</I>. This is called
 *  <I>extended</I> gcd computation, and <I>f1, f2</I> are called
 *  <I>B&eacute;zout factors</I> or <I>cofactors</I>.
 *
 *  CGALially, computation is performed ``denominator-free'' if
 *  supported by the coefficient type via \c CGAL::Fraction_traits
 *  (using \c pseudo_gcdex() ), otherwise the euclidean remainder
 *  sequence is used.
 *
 *  \pre \c NT must be a \c Field.
 *
 *  The result <I>d</I> is unit-normal,
 *  i.e. <I>d</I><TT>.lcoeff() == NT(1)</TT>.
 *
 */
template <class NT> 
inline
Polynomial<NT> gcdex(
        Polynomial<NT> p1, Polynomial<NT> p2,
        Polynomial<NT>& f1, Polynomial<NT>& f2
) {
    typedef typename CGAL::Fraction_traits< Polynomial<NT> >
        ::Is_fraction Is_fraction;
    return internal::gcdex_(p1, p2, f1, f2, Is_fraction());
}


/*! \ingroup CGAL_Polynomial
 *  \relates CGAL::Polynomial
 *  \brief compute gcd with ``almost'' cofactors
 *
 *  This is a variant of \c exgcd() for use over non-field \c NT.
 *  It computes the gcd of polynomials \c p1 and \c p2
 *  along with two other polynomials \c f1 and \c f2 and a scalar \c v
 *  such that \e v * gcd(\e p1, \e p2) = <I>f1*p1 + f2*p2</I>,
 *  using the subresultant remainder sequence. That \c NT is not a field
 *  implies that one cannot achieve \e v = 1 for all inputs.
 *
 *  \pre \c NT must be a \c UFDomain of scalars (not polynomials).
 *
 *  The result is unit-normal.
 *
 */
template <class NT>
inline
Polynomial<NT> pseudo_gcdex(
#ifdef DOXYGEN_RUNNING
        Polynomial<NT> p1, Polynomial<NT> p2,
        Polynomial<NT>& f2, Polynomial<NT>& f2, NT& v
#else
        Polynomial<NT> x, Polynomial<NT> y,
        Polynomial<NT>& xf, Polynomial<NT>& yf, NT& vf
#endif // DOXYGEN_RUNNING
) {
    /* implemented using the extended subresultant algorithm
     * for gcd computation with Bezout factors
     *
     * To understand this, you need to understand the computation of
     * cofactors as in the basic extended Euclidean algorithm (see
     * the code above of gcdex_(..., Field_tag)), and the subresultant
     * gcd algorithm, see gcd_(..., Unique_factorization_domain_tag).
     *
     * The crucial point of the combination of both is the observation
     * that the subresultant factor (called rho here) divided out of the
     * new remainder in each step can also be divided out of the
     * cofactors.
     */

    typedef Polynomial<NT> POLY;
    typename Algebraic_structure_traits<NT>::Integral_division idiv;
    typename Algebraic_structure_traits<NT>::Gcd          gcd;

    // handle trivial cases
    if (x.is_zero()) {
        if (y.is_zero()) CGAL_error_msg("gcdex(0,0) is undefined");
        xf = POLY(0); yf = POLY(1); vf = y.unit_part();
        return y / vf;
    }
    if (y.is_zero()) {
        xf = POLY(1); yf = POLY(0); vf = x.unit_part();
        return x / vf;
    }
    bool swapped = x.degree() < y.degree();
    if (swapped) { POLY t = x; x = y; y = t; }

    // compute gcd of content
    NT xcont = x.content(); NT ycont = y.content();
    NT gcdcont = gcd(xcont, ycont);

    // compute gcd of primitive parts
    POLY xprim = x / xcont; POLY yprim = y / ycont;
    POLY u = xprim, v = yprim, q, r;
    POLY m11(1), m21(0), m21old;
    NT g(1), h(1), d, rho;
    for (;;) {
        int delta = u.degree() - v.degree();
        POLY::pseudo_division(u, v, q, r, d);
        CGAL_assertion(d == ipower(v.lcoeff(), delta+1));
        if (r.is_zero()) break;
        rho = g * ipower(h, delta);
        u = v; v = r / rho;
        m21old = m21; m21 = (d*m11 - q*m21) / rho; m11 = m21old;
        /* The transition from (u, v) to (v, r/rho) corresponds
         * to multiplication with the matrix
         *   __1__ (0 rho)
         *    rho  (d  -q)
         * The comments and correctness arguments from
         * gcdex(..., Field_tag) apply analogously.
         */
        g = u.lcoeff();
        CGAL::internal::hgdelta_update(h, g, delta);
        if (r.degree() == 0) break;
    }

    // obtain v == m21*xprim + m22*yprim
    // the correct m21 was already computed above
    POLY m22;
    POLY::euclidean_division(v - m21*xprim, yprim, m22, r);
    CGAL_assertion(r.is_zero());

    // now obtain gcd(x,y) == gcdcont * v/v.content() == (m21*x + m22*y)/denom
    NT vcont = v.content(), vup = v.unit_part();
    v /= vup * vcont; v *= gcdcont;
    m21 *= ycont; m22 *= xcont;
    vf = idiv(xcont, gcdcont) * ycont * (vup * vcont);
    CGAL_assertion(vf * v == m21*x + m22*y);

    // return results
    if (swapped) {
        xf = m22; yf = m21;
    } else {
        xf = m21; yf = m22;
    }
    return v;
}





} // namespace CGAL

#endif // CGAL_POLYNOMIAL_GCD_H

// EOF
