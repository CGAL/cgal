// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     :  
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_POLYNOMIAL_UNIVARIATE_POLYNOMIAL_UTILS_H
#define CGAL_POLYNOMIAL_UNIVARIATE_POLYNOMIAL_UTILS_H

#include <CGAL/Polynomial.h>

CGAL_BEGIN_NAMESPACE

//! return an upper bound on the absolute value of all real roots of \c P.
/*! The upper bound is a power of two. Only works for univariate polynomials.
 *  \pre \c NT must be \c RealComparable.
 *  \relates NiX::Polynomial
 */
template <class NT>
NT weak_upper_root_bound(const Polynomial<NT>& P) { 
    // code comes from Kurt Mehlhorn
    // see [Mignotte, 1992], p.144 for a proof
    CGAL_precondition(Polynomial_traits_d<NT>::d == 0);
    typename Real_embeddable_traits<NT>::Abs abs;
    const int n = P.degree();
    NT x(1); 
    NT val;
    for (;;) {
        val = -abs(P[n]);
        for (int i = n-1; i >= 0; i--) {
            val = val*x + abs(P[i]);
        }
        if (val < NT(0)) return x;
        x *= NT(2);
    }
}

//! return the number of sign variations in the coefficient sequence of \c P.
/*! This is the number of sign changes (+ to - or - to +) in the
 *  coefficient sequence of the polynomial, ignoring zeroes.
 *  Only meaningful for univariate polynomials.
 *  \pre \c NT must be \c RealComparable.
 *  \relates NiX::Polynomial
 */
template <class NT>
int sign_variations(const Polynomial<NT>& P) { 
    typename Real_embeddable_traits<NT>::Sign sign;
    const int n = P.degree();
    int variations = 0;
    int old_sign = sign(P[n]); // never zero unless P is zero
    for (int i = n-1; i >= 0; i--) {
        int s = sign(P[i]);
        if (s == 0) continue;
        if (old_sign != s) {
            old_sign = s;
            variations++;
        }
    }
    return variations;
}

CGAL_END_NAMESPACE


#endif // CGAL_POLYNOMIAL_UNIVARIATE_POLYNOMIAL_UTILS_H
