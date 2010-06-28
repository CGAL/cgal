// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Polynomial/include/CGAL/Polynomial.h $
// $Id: Polynomial.h 47254 2008-12-06 21:18:27Z afabri $
// 
//
// Author(s)     :  Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_POLYNOMIAL_UNIVARIATE_POLYNOMIAL_UTILS_H
#define CGAL_POLYNOMIAL_UNIVARIATE_POLYNOMIAL_UTILS_H

#include <CGAL/Polynomial.h>

namespace CGAL {

namespace internal {
    //! return an upper bound on the absolute value of all real roots of \c P.
    /*! The upper bound is a power of two. Only works for univariate polynomials.
     *  \pre \c NT must be \c RealComparable.
     *  \relates CGAL::Polynomial
     */
    template <class NT>
    NT weak_upper_root_bound(const Polynomial<NT>& P) { 
        // code comes from Kurt Mehlhorn
        // see [Mignotte, 1992], p.144 for a proof
        CGAL_precondition(Polynomial_traits_d<NT>::d == 0);
        typename Real_embeddable_traits<NT>::Abs abs;
        const int n = CGAL::degree(P);
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
     *  \relates CGAL::Polynomial
     */
    template <class NT>
    int sign_variations(const Polynomial<NT>& P) { 
        const int n = CGAL::degree(P);
        int variations = 0;
        int old_sign = CGAL::sign(P[n]); // never zero unless P is zero
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

    /*! \ingroup CGAL_polynomial_utils
     *  \brief checks whether a univariate polynomial is square-free
     */    
    
    /*template < class NT >
    bool is_square_free(const Polynomial<NT>& p) {
        if( may_have_multiple_factor(p) ) {
            Polynomial<NT> d = p; d.diff();
            return CGAL::degree(gcd_utcf(p, d)) == 0;
        } else {
            return true;
        }
    }
    
    template< class NT >
    bool is_square_free( const Polynomail< Polynomial< NT > >& ) { 
        
        
        return true;
    } */   

} // namespace internal

} //namespace CGAL


#endif // CGAL_POLYNOMIAL_UNIVARIATE_POLYNOMIAL_UTILS_H
