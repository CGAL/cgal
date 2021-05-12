// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_REFINE_ZERO_AGAINST_H
#define CGAL_REFINE_ZERO_AGAINST_H

#include <CGAL/basic.h>

#include <CGAL/Polynomial.h>

namespace CGAL {

namespace internal {
/* computes an upper bound on the number of roots ]low,high[ using Descartes'
 * Sign Rule
*/
template <class Polynomial, class Field>
int descartes(Polynomial& p, const Field& low,const Field& high){
// decompose interval length and upper bound
    CGAL_precondition(low<high);
    typedef typename Polynomial::NT Coefficient;
    typedef typename Fraction_traits<Field>::Numerator_type Numerator;
    typedef typename Fraction_traits<Field>::Denominator_type Denominator;

    typename Fraction_traits<Field>::Decompose decomp;
    //typename Algebraic_structure_traits<Field>::Simplify simplify;

    //simplify(low);
    //simplify(high);

    Numerator num_high, num_low_sub_high;
    Denominator den_high, den_low_sub_high;

    decomp(high, num_high, den_high);
    decomp(low - high, num_low_sub_high, den_low_sub_high);

    Coefficient tmp(num_high);
    // apply Descartes' rule to count zeros of q in ]low,high[
    Polynomial transf = // q(high + (low-high)/(1+x))
      CGAL::translate_by_one(
          CGAL::reversal(
              CGAL::scale_homogeneous(
                  CGAL::translate_homogeneous(p
                      ,Coefficient(num_high)
                      ,Coefficient(den_high))
                  ,Coefficient(num_low_sub_high)
                  ,Coefficient(den_low_sub_high)
              )
          )
      );
    return sign_variations(transf);
}

/*! \ingroup \NiX_univariate_polynomial_utils
 *  \brief refine isolating interval for \c p w.r.t \c q
 *
 *  This function refines the interval ]<TT>low</TT>, <TT>high</TT>[
 *  such that it does not contain any zero of \c q different from the
 *  unique zero of \c p in ]<TT>low</TT>, <TT>high</TT>[.  It is returned
 *  whether \c q has a zero in ]<TT>low</TT>, <TT>high</TT>[ equal to
 *  that of \c p . Note that zeroes of \c q at the boundaries are
 *  ignored.
 *
 *  This function is implemented using bisection and Descartes' Rule.
 *  If the interval boundaries have denominators 2<SUP>k</SUP>, then
 *  this property will still hold after refinement.  Although this
 *  function works similar to \c NiX::Algebraic_real<>.compare() ,
 *  it is different insofar that it always maintains an open interval
 *  and never simplifies.
 *
 *  \pre Both polynomials must be square-free. \c p must not vanish at the
 *  interval boundaries \c low and \c high .
 *
 *  \todo Provide a means to let an \c NiX::Algebraic_real benefit
 *  from the interval refinement if it is the origin of the respective
 *  input data.
 */
template <class Polynomial, class Field>
bool refine_zero_against(Field& low, Field& high, Polynomial p, Polynomial q) {
    typename Algebraic_structure_traits<Field>::Simplify simplify;

    CGAL_precondition(low < high);
    CGAL_precondition(CGAL::degree(p) > 0);
    CGAL_precondition((CGAL::degree(q) >= 0) && !q.is_zero());

    if (CGAL::degree(q) == 0) return false;

    CGAL::Sign sign_p_low  = p.sign_at(low);
    CGAL::Sign sign_p_high = p.sign_at(high);
    CGAL_precondition(sign_p_low  != CGAL::ZERO);
    CGAL_precondition(sign_p_high != CGAL::ZERO);
    CGAL_precondition(sign_p_high != sign_p_low);

    Polynomial gcd_pq; // computed below if necessary

    for (;;) {
        int sv = internal::descartes(q,low,high);
        CGAL_assertion(sv >= 0);

        if (sv == 0) {
            // q has no zero in ]low,high[
            return false;
        } else if (sv == 1) {
            if (CGAL::degree(gcd_pq) < 0) {
                if (may_have_common_factor(p, q)) {
                    gcd_pq = gcd_utcf(p, q);
                } else {
                    gcd_pq = Polynomial(1);
                }
            }
            std::cout << CGAL::to_double(low) << " "
                      << CGAL::to_double(high) << " "
                      << CGAL::degree(gcd_pq) << " "
                      << gcd_pq
                      << std::endl;
            if (CGAL::degree(gcd_pq) > 0 // constant poly cannot change sign
            && gcd_pq.sign_at(low) != gcd_pq.sign_at(high)) {
                // q has exactly one zero in ]low,high[
                // and it's equal to that of p
                return true;
            }
        }
        // q may have a zero in ]low,high[ not equal to that of p
        Field mid = (low+high)/Field(2);
        CGAL::Sign s = p.sign_at(mid);
        if (s == CGAL::ZERO) {
            mid = (low+mid)/Field(2);
            simplify(mid);
            s = p.sign_at(mid);
        }
        CGAL_postcondition(s != CGAL::ZERO);

        if (s == sign_p_low) {
            low = mid;
            sign_p_low = s;
        } else {
            CGAL_postcondition(s == sign_p_high);
            high = mid;
            sign_p_high = s;
        }
    }
}


// Uses refine_zero_against first and refines the interval further, if any
//  of the interval boarders has sign zero.
template < class Polynomial, class Field >
static bool strong_refine_zero_against(Field& low, Field& high,
                                       Polynomial p, Polynomial q){
    std::cout << "comp has_common_root" << std::endl;

    bool has_common_root = refine_zero_against(low,high,p,q);

    std::cout << "done, " << has_common_root << std::endl;

    CGAL::Sign sign_p_low = p.sign_at(low);
    CGAL::Sign sign_p_high = p.sign_at(high);

    Field mid;
    CGAL::Sign s;

    while ((q.sign_at(low)==CGAL::ZERO)||(q.sign_at(high)==CGAL::ZERO)) {
        mid = (low+high)/Field(2);
        simplify(mid);
        s = p.sign_at(mid);
        if (s == CGAL::ZERO) {
            mid = (low+mid)/Field(2);
            simplify(mid);
            s = p.sign_at(mid);
        }
        CGAL_postcondition(s != CGAL::ZERO);

        if (s == sign_p_low) {
            low = mid;
            sign_p_low = s; //bogus?
        }
        else {
            CGAL_assertion(s == sign_p_high);
            high = mid;
            sign_p_high = s; //bogus?
        }
    }

    return has_common_root;
}

} //namespace internal

} //namespace CGAL

#endif //CGAL_REFINE_ZERO_AGAINST_H
