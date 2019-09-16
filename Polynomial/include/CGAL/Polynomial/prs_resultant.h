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

// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

/*! \file CGAL/prs_resultant.h
 *  \brief Resultant computation via polynomial remainder sequences (PRS)
 *
 */

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>

#include <CGAL/ipower.h>
#include <CGAL/Polynomial/hgdelta_update.h>

#ifndef CGAL_POLYNOMIAL_PRS_RESULTANT_H
#define CGAL_POLYNOMIAL_PRS_RESULTANT_H

namespace CGAL {

template <class NT> inline
NT prs_resultant_integral_domain(Polynomial<NT> A, Polynomial<NT> B) {
    // implemented using the subresultant algorithm for resultant computation
    // see [Cohen, 1993], algorithm 3.3.7

    if (A.is_zero() || B.is_zero()) return NT(0);

    int signflip;
    if (A.degree() < B.degree()) {
        Polynomial<NT> T = A; A = B; B = T;
        signflip = (A.degree() & B.degree() & 1);
    } else {
        signflip = 0;
    }
    
    typedef CGAL::Scalar_factor_traits<Polynomial<NT> > SFT;
    typedef typename SFT::Scalar Scalar; 
    typename SFT::Scalar_factor scalar_factor;
    typename CGAL::Coercion_traits<Scalar, NT>::Cast cast_scalar_nt; 
    
    Scalar a = scalar_factor(A), b = scalar_factor(B);
    NT g(1), h(1);
    NT t = cast_scalar_nt (CGAL::ipower(a, B.degree()) * CGAL::ipower(b, A.degree()));

    Polynomial<NT> Q, R; NT d;
    int delta;

    A /= cast_scalar_nt(a); B /= cast_scalar_nt(b);
    do {
        signflip ^= (A.degree() & B.degree() & 1);
        Polynomial<NT>::pseudo_division(A, B, Q, R, d);
        delta = A.degree() - B.degree();
        CGAL_expensive_assertion_code
          (typedef typename CGAL::Algebraic_structure_traits<NT>::Is_exact
           Is_exact;)
        CGAL_expensive_assertion(CGAL::check_tag(Is_exact()) == false
          || d == CGAL::ipower(B.lcoeff(), delta + 1) );
        A = B;
        B = R / (g * CGAL::ipower(h, delta));
        g = A.lcoeff();
        // h = h^(1-delta) * g^delta
        internal::hgdelta_update(h, g, delta);
    } while (B.degree() > 0);
    // h = h^(1-deg(A)) * lcoeff(B)^deg(A)
    delta = A.degree();
    g = B.lcoeff();
    internal::hgdelta_update(h, g, delta);
    h = signflip ? -(t*h) : t*h;
    typename Algebraic_structure_traits<NT>::Simplify simplify;
    simplify(h);
    return h;
}



template <class NT> inline
NT prs_resultant_ufd(Polynomial<NT> A, Polynomial<NT> B) {
    // implemented using the subresultant algorithm for resultant computation
    // see [Cohen, 1993], algorithm 3.3.7

    if (A.is_zero() || B.is_zero()) return NT(0);

    int signflip;
    if (A.degree() < B.degree()) {
        Polynomial<NT> T = A; A = B; B = T;
        signflip = (A.degree() & B.degree() & 1);
    } else {
        signflip = 0;
    }

    NT a = A.content(), b = B.content();
    NT g(1), h(1), t = CGAL::ipower(a, B.degree()) * CGAL::ipower(b, A.degree());
    Polynomial<NT> Q, R; NT d;
    int delta;

    A /= a; B /= b;
    do {
        signflip ^= (A.degree() & B.degree() & 1);
        Polynomial<NT>::pseudo_division(A, B, Q, R, d);
        delta = A.degree() - B.degree();
        CGAL_expensive_assertion_code
          (typedef typename CGAL::Algebraic_structure_traits<NT>::Is_exact
           Is_exact;)
        CGAL_expensive_assertion(CGAL::check_tag(Is_exact()) == false
          || d == CGAL::ipower(B.lcoeff(), delta + 1) );
        A = B;
        B = R / (g * CGAL::ipower(h, delta));
        g = A.lcoeff();
        // h = h^(1-delta) * g^delta
        internal::hgdelta_update(h, g, delta);
    } while (B.degree() > 0);
    // h = h^(1-deg(A)) * lcoeff(B)^deg(A)
    delta = A.degree();
    g = B.lcoeff();
    internal::hgdelta_update(h, g, delta);
    if (signflip)
      h = -(t*h);
    else
      h = t*h;
    typename Algebraic_structure_traits<NT>::Simplify simplify;
    simplify(h);
    return h;
}

template <class NT> inline
NT prs_resultant_field(Polynomial<NT> A, Polynomial<NT> B) {
    // implemented using the Euclidean algorithm for resultant computation
    // compare [Cox et al, 1997], p.157

    if (A.is_zero() || B.is_zero()) return NT(0);

    int signflip;
    if (A.degree() < B.degree()) {
        Polynomial<NT> T = A; A = B; B = T;
        signflip = (A.degree() & B.degree() & 1);
    } else {
        signflip = 0;
    }

    NT res(1);
    Polynomial<NT> Q, R;
    while (B.degree() > 0) {
        signflip ^= (A.degree() & B.degree() & 1);
        Polynomial<NT>::euclidean_division(A, B, Q, R);
        res *= CGAL::ipower(B.lcoeff(), A.degree() - R.degree());
        A = B;
        B = R;
    }
    res = CGAL::ipower(B.lcoeff(), A.degree()) * (signflip ? -res : res);
    typename Algebraic_structure_traits<NT>::Simplify simplify;
    simplify(res);
    return res;
}

// definition follows below
template <class NT> inline
NT prs_resultant_decompose(Polynomial<NT> A, Polynomial<NT> B);

namespace INTERN_PRS_RESULTANT {
    template <class NT> inline
    NT prs_resultant_(Polynomial<NT> A, Polynomial<NT> B, ::CGAL::Tag_false) {
        return prs_resultant_field(A, B);
    }

    template <class NT> inline
    NT prs_resultant_(Polynomial<NT> A, Polynomial<NT> B, ::CGAL::Tag_true) {
        return prs_resultant_decompose(A, B);
    }

    template <class NT> inline
    NT prs_resultant_(Polynomial<NT> A, Polynomial<NT> B, Field_tag) {
        typedef typename Fraction_traits<NT>::Is_fraction Is_decomposable;
        return prs_resultant_(A, B, Is_decomposable());     
    }

    template <class NT> inline
    NT prs_resultant_(Polynomial<NT> A, Polynomial<NT> B, Unique_factorization_domain_tag) {
        return prs_resultant_ufd(A, B);
    }
} // namespace internal

template <class NT> inline
NT prs_resultant_decompose(Polynomial<NT> A, Polynomial<NT> B){
    typedef Polynomial<NT> POLY;
    typedef typename Fraction_traits<POLY>::Numerator_type INTPOLY;
    typedef typename Fraction_traits<POLY>::Denominator_type DENOM;
    typename Fraction_traits<POLY>::Decompose decompose;
    typedef typename INTPOLY::NT RES;
    
    DENOM a, b;
    A.simplify_coefficients();
    B.simplify_coefficients();
    INTPOLY A0; decompose(A,A0,a);
    INTPOLY B0; decompose(B,B0,b);
    DENOM c = CGAL::ipower(a, B.degree()) * CGAL::ipower(b, A.degree());
    typedef typename Algebraic_structure_traits<RES>::Algebraic_category Algebraic_category;
    RES res0 = INTERN_PRS_RESULTANT::prs_resultant_(A0, B0, Algebraic_category());
    typename Fraction_traits<NT>::Compose comp_frac;
    NT res = comp_frac(res0, c);
    typename Algebraic_structure_traits<NT>::Simplify simplify;
    simplify(res);
    return res;
}


/*! \ingroup CGAL_Polynomial
 *  \relates CGAL::Polynomial
 *  \brief compute the resultant of polynomials \c A and \c B
 *
 *  The resultant of two polynomials is computed from their
 *  polynomial remainder sequence (PRS), in the Euclidean or
 *  subresultant version. This depends on the coefficient type:
 *  If \c NT is a \c UFDomain , the subresultant PRS is formed.
 *  If \c NT is a \c Field that is not decomposable (see
 *  \c CGAL::Fraction_traits ), then a Euclidean PRS is formed.
 *  If \c NT is a \c Field that is decomposable, then the
 *  \c Numerator must be a \c UFDomain, and the subresultant
 *  PRS is formed for the decomposed polynomials.
 *
 *  Using \c CGAL::hybrid_bezout_subresultant() may be faster in some cases
 *  and works for non-UFDomains, too.
 *  Using \c CGAL::resultant() from \c CGAL/resultant.h
 *  chooses automatically among these alternative methods of resultant
 *  computation for you.
 *
 *  For the benefit of those who want to do their own template
 *  metaprogramming to choose the method of resultant computation,
 *  the three variants of resultant computation from a PRS
 *  can be called directly as \c prs_resultant_field() ,
 *  \c prs_resultant_ufd() and \c prs_resultant_decompose() .
 *  <b>Do not use them directly unless you know what you are doing!</b>
 *
 */
template <class NT> inline
NT prs_resultant(Polynomial<NT> A, Polynomial<NT> B) {
    typedef typename Algebraic_structure_traits<NT>::Algebraic_category
                                                                   Algebraic_category;
    return INTERN_PRS_RESULTANT::prs_resultant_(A, B, Algebraic_category());     
}

} //namespace CGAL

#endif // CGAL_POLYNOMIAL_PRS_RESULTANT_H

// EOF
