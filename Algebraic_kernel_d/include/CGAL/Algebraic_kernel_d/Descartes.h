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

/*! \file NiX/Descartes.h
  \brief Defines class NiX::Descartes.

  Isolate real roots of polynomials.

  This file provides a class to isolate real roots of polynomials,
  using the algorithm based on the method of Descartes.

  The polynomial has to be a univariat polynomial over any number
  type which is contained in the real numbers.

*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_DESCARTES_H
#define CGAL_ALGEBRAIC_KERNEL_D_DESCARTES_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>

#include <CGAL/Algebraic_kernel_d/univariate_polynomial_utils.h>
#include <CGAL/Algebraic_kernel_d/construct_binary.h>

#define POLYNOMIAL_REBIND( coeff ) \
    typename CGAL::Polynomial_traits_d<Polynomial>::template \
    Rebind<coeff,1>::Other::Type


namespace CGAL {

namespace internal {

/*! \ingroup NiX_Algebraic_real
 *  \brief A model of concept RealRootIsolator.
 */
template <class Polynomial_, class Rational_>
class Descartes {
    typedef CGAL::Fraction_traits<Polynomial_> FT_poly;
    typedef Fraction_traits<Rational_> FT_rat;
public:
    //! First template parameter
    typedef Polynomial_ Polynomial;
    //! Second template parameter
    typedef Rational_ Rational;
    //! Bound type of the isolating intervals
    typedef Rational_ Bound;
    // Integer or Numerator/Denominator type of bound.
    typedef typename CGAL::Fraction_traits<Rational>::Numerator_type Integer;
private:
    typedef typename Polynomial::NT Coeff;
    typedef Integer IT;

    Polynomial poly_;
    int number_of_real_roots_;
    IT* numerator;
    IT* denominator_exponent;
    bool* is_exact;
    IT LEFT,SCALE,DENOM;
    bool is_strong_;
    int k;
    bool interval_given;

public:
    /*! \brief Constructor from univariate square free polynomial.

    The RealRootIsolator provides isolating intervals for the real
    roots of the polynomial.
    \pre the polynomial is square free
    */
    Descartes(const Polynomial& P = Polynomial(Coeff(0)),
            bool is_strong = false,
            int kk = 2)
        : poly_(P) ,
          is_strong_(is_strong),
          k(kk),
          interval_given(false) {

        numerator = new IT[CGAL::degree(P)];
        denominator_exponent = new IT[CGAL::degree(P)];
        is_exact = new bool[CGAL::degree(P)];
        number_of_real_roots_ = 0;
        if(CGAL::degree(P) == 0)
            {
                if(P.is_zero()) number_of_real_roots_ = -1;
                return;
            }

        intern_decompose(poly_,typename FT_poly::Is_fraction());
    }

    // constructor for coefficient types \c Coeff with given interval
    // (experimental)
    Descartes(const Polynomial& P,
            const Rational& left,
            const Rational& right,
            bool is_strong = false,
            int kk = 2)
        : poly_(P) ,
          is_strong_(is_strong),
          k(kk),
          interval_given(true) {

        numerator = new IT[CGAL::degree(P)];
        denominator_exponent = new IT[CGAL::degree(P)];
        is_exact = new bool[CGAL::degree(P)];
        number_of_real_roots_ = 0;
        if(CGAL::degree(P) == 0)
            {
                if(P.is_zero()) number_of_real_roots_ = -1;
                return;
            }
        typename FT_rat::Decompose decompose;
        typedef typename FT_rat::Numerator Numerator;
        typedef typename FT_rat::Denominator Denominator;
        Numerator numleft, numright;
        Denominator denleft, denright;

        decompose(left,numleft,denleft);
        decompose(right,numright,denright);

        LEFT = numleft * denright;
        SCALE = numright * denleft - LEFT;
        DENOM = denleft * denright;
        poly_.scale_down(denleft*denright);

        intern_decompose(poly_,typename FT_poly::Is_decomposable());
    }

    //! copy constructor
    Descartes(const Descartes& D)
        : poly_(D.poly_),
          number_of_real_roots_(D.number_of_real_roots_),
          LEFT(D.LEFT),
          SCALE(D.SCALE),
          DENOM(D.DENOM),
          is_strong_(D.is_strong_),
          k(D.k),
          interval_given(D.interval_given) {

        numerator = new IT[CGAL::degree(poly_)];
        denominator_exponent = new IT[CGAL::degree(poly_)];
        is_exact = new bool[CGAL::degree(poly_)];
        for(int i=0; i<number_of_real_roots(); i++)
            {
                numerator[i] = D.numerator[i];
                denominator_exponent[i] = D.denominator_exponent[i];
                is_exact[i] = D.is_exact[i];
            }
    }

    // destructor
    ~Descartes() {
        delete[] numerator;
        delete[] denominator_exponent;
        delete[] is_exact;
    }

public: // functions

    /*! \brief returns the defining polynomial*/
    Polynomial polynomial() const { return poly_; }

    //! returns the number of real roots
    int number_of_real_roots() const { return number_of_real_roots_; }

    /*! \brief returns true if the isolating interval is degenerated to a
      single point.

      If is_exact_root(i) is true,
      then left_bound(int i) equals  \f$root_i\f$. \n
      If is_exact_root(i) is true,
      then right_bound(int i) equals  \f$root_i\f$. \n
    */
    bool is_exact_root(int i) const { return is_exact[i]; }



public:


    void left_bound(int i, IT& numerator_, IT& denominator_) const {
        CGAL_assertion(i >= 0 && i < number_of_real_roots_);
        construct_binary(denominator_exponent[i], denominator_);
        numerator_= SCALE * numerator[i] + LEFT * denominator_;
        denominator_ = denominator_ * DENOM;
    }


    void right_bound(int i,IT& numerator_, IT& denominator_) const {
        CGAL_assertion(i >= 0 && i < number_of_real_roots_);
        if(is_exact[i]){
            return left_bound(i,numerator_,denominator_);
        }
        else{
            construct_binary(denominator_exponent[i],denominator_);
            numerator_= SCALE * (numerator[i]+1) + LEFT * denominator_;
            denominator_ = denominator_ * DENOM;
        }
    }
public:

    /*! \brief returns  \f${l_i}\f$ the left bound of the isolating interval
      for root  \f$root_{i}\f$.

      In case is_exact_root(i) is true,  \f$l_i = root_{i}\f$,\n
      otherwise:  \f$l_i < root_{i}\f$.

      If  \f$i-1>=0\f$, then  \f$l_i > root_{i-1}\f$. \n
      If  \f$i-1>=0\f$, then  \f$l_i >= r_{i-1}\f$,
      the right bound of  \f$root_{i-1}\f$\n

      \pre 0 <= i < number_of_real_roots()
    */
    Rational left_bound(int i) const {
        IT numerator_, denominator_;
        left_bound(i,numerator_,denominator_);
        return Rational(numerator_) / Rational(denominator_);
    }

    /*! \brief returns  \f${r_i}\f$ the right bound of the isolating interval
      for root  \f$root_{i}\f$.

      In case is_exact_root(i) is true,  \f$r_i = root_{i}\f$,\n
      otherwise:  \f$r_i > root_{i}\f$.

      If  \f$i+1< n \f$, then  \f$r_i < root_{i+1}\f$,
      where \f$n\f$ is number of real roots.\n
      If  \f$i+1< n \f$, then  \f$r_i <= l_{i+1}\f$,
      the left bound of  \f$root_{i+1}\f$\n


      \pre 0 <= i < number_of_real_roots()
    */
    Rational right_bound(int i) const {
        IT numerator_, denominator_;
        right_bound(i,numerator_,denominator_);
        return Rational(numerator_) / Rational(denominator_);
    }

private:
    void intern_decompose( Polynomial P_, ::CGAL::Tag_true){
        typename FT_poly::Decompose decompose;
        typename FT_poly::Numerator_type NumP;
        typename FT_poly::Denominator_type dummy;

        decompose(P_,NumP,dummy);
        init_with(NumP);
    }

    void intern_decompose( Polynomial P, ::CGAL::Tag_false){
        init_with(P);
    }


    template<class Polynomial__>
    void init_with(const Polynomial__& P){
        typedef typename  Polynomial__::NT Coeff;
        if(!interval_given)
            {
                LEFT = -weak_upper_root_bound<Coeff>(P);
                SCALE = - LEFT * IT(2);
                DENOM = IT(1);
            }
        Polynomial__ R = ::CGAL::translate(P,Coeff(LEFT));
        Polynomial__ Q = ::CGAL::scale_up(R,Coeff(SCALE));
        zero_one_descartes<Coeff>(Q,0,0);
    }


    //! returns the polynomial $(1 + x)^n P(1/(1 + x))$.
    template <class Coeff__>
    /*
    typename
    CGAL::Polynomial_traits_d<Polynomial>
    ::template Rebind<Coeff__,1>::Other::Type
    */
    POLYNOMIAL_REBIND(Coeff__)
        variation_transformation(const POLYNOMIAL_REBIND(Coeff__)& P) {
        POLYNOMIAL_REBIND(Coeff__) R = reversal(P);
        return translate_by_one(R);
    }

    //! Returns an upper bound on the absolute value of all roots of $P$.
    /*! The upper bound is a power of two. Only works for univariate
     * polynomials.
     */
    template <class Coeff__>
    IT weak_upper_root_bound(const POLYNOMIAL_REBIND(Coeff__)& P) {

        typename Real_embeddable_traits<Coeff__>::Abs abs;
        const int n = CGAL::degree(P);
        IT r(1);  // return value
        Coeff__ x(1);  // needed to "evaluate" the polynomial
        Coeff__ val;
        for (;;) {
            val = -abs(P[n]);
            for (int i = n-1; i >= 0; i--) {
                val = val*x + abs(P[i]);
            }
            if (val < Coeff__(0)) return r;
            r *= IT(2);
            x = Coeff__(r);
        }
    }

    //! tests if the polynomial has no root in the interval.
    template <class Coeff__>
    bool not_zero_in_interval(const POLYNOMIAL_REBIND(Coeff__)& P)
    {
        if(CGAL::degree(P) == 0) return true;
        if(internal::sign_variations(variation_transformation<Coeff__>(P)) != 0)
            return false;
        return (P[0] != Coeff__(0) && P.evaluate(Coeff__(1)) != Coeff__(0));
    }
    //! Descartes algoritm to determine isolating intervals for the roots
    //! lying in the interval (0,1).
    // The parameters $(i,D)$ describe the interval $(i/2^D, (i+1)/2^D)$.
    // Here $0\leq i < 2^D$.
    template <class Coeff__>
    void zero_one_descartes(const POLYNOMIAL_REBIND(Coeff__)& P,
            IT i, IT D) {
        // Determine the number of sign variations of the transformed
        // polynomial $(1+x)^nP(1/(1+x))$. This gives the number of
        // roots of $P$ in $(0,1)$.

        POLYNOMIAL_REBIND(Coeff__) R = variation_transformation<Coeff__>(P);
        int descarte = sign_variations(R);

        // no root
        if ( descarte == 0 ) return;

        // exactly one root
        // Note the termination criterion $P(0)\neq 0$ and $P(1)\neq 0$.
        // This ensures that the given interval is an isolating interval.
        if ( descarte == 1
                && P[0] != Coeff__(0)
                && P.evaluate(Coeff__(1)) != Coeff__(0) ) {
            if(is_strong_) {
                strong_zero_one_descartes<Coeff__>(P,i,D);
                return;
            }
            else {
                numerator[number_of_real_roots_] = i;
                denominator_exponent[number_of_real_roots_] = D;
                is_exact[number_of_real_roots_] = false;
                number_of_real_roots_++;
                return;
            }
        }

        // more than one root
        // Refine the interval.
        i = 2*i; D = D+1;

        // Transform the polynomial such that the first half of the interval
        // is mapped to the unit interval.
        POLYNOMIAL_REBIND(Coeff__) Q = scale_down(P,Coeff__(2));

        // Consider the first half of the interval.
        zero_one_descartes<Coeff__>(Q,i,D);

        // Test if the polynomial is zero at the midpoint of the interval
        POLYNOMIAL_REBIND(Coeff__)  S = translate_by_one(Q);
        if ( S[0] == Coeff__(0) ) {
            numerator[number_of_real_roots_] = i + 1;
            denominator_exponent[number_of_real_roots_] = D;
            is_exact[number_of_real_roots_] = true;
            number_of_real_roots_++;
        }

        // Consider the second half of the interval.
        zero_one_descartes<Coeff__>(S,i+1,D);
    }


    //! Strong Descartes algoritm to determine isolating intervals for the
    //! roots lying in the interval (0,1), where the first
    //! derivative have no sign change. \pre $P$ has only one root in the
    //! interval given by $(i,D)$.
    // The parameters $(i,D)$ describe the interval $(i/2^D, (i+1)/2^D)$.
    // Here $0\leq i < D$.
    template <class Coeff__>
    void strong_zero_one_descartes(const POLYNOMIAL_REBIND(Coeff__)& P,
            IT i, IT D) {

        // Test if the polynomial P' has no roots in the
        // interval. For further use in Newton, the interval should be not
        // too large.

        // test if isolating interval is smaller than epsilon
        // [l,r]  ->  r-l < epsilon
        // l = (r-l) * i/2^D + l
        // r = (r-l) * (i+1)/2^D + l
        // r-l = (r-l) * 1/2^D
        // r-l < epsilon = 2^(-k)
        // <=> (r-l) * 1/2^D < 2^(-k)
        // <=> 2^D > (r-l) / 2^(-k)
        // <=> 2^D > (r-l) * 2^k

      POLYNOMIAL_REBIND(Coeff__) PP = CGAL::differentiate(P);
        if(not_zero_in_interval<Coeff__>(PP)) { // P'
            IT tmp;
            construct_binary(D-k, tmp);  // tmp = 2^{D-k}
            if(tmp * DENOM > SCALE ) {
                numerator[number_of_real_roots_] = i;
                denominator_exponent[number_of_real_roots_] = D;
                is_exact[number_of_real_roots_] = false;
                number_of_real_roots_++;
                return;
            }
        }

        // either $P'$ fails the test,
        // or the interval is too large
        // Refine the interval.
        i = 2*i; D = D+1;

        // Transform the polynomial such that the first half of the interval
        // is mapped to the unit interval.
        POLYNOMIAL_REBIND(Coeff__) Q = scale_down(P,Coeff__(2));

        // Test if the polynomial is zero at the midpoint of the interval
        POLYNOMIAL_REBIND(Coeff__)  S = translate_by_one(Q);
        if ( S[0] == Coeff__(0) ) {
            numerator[number_of_real_roots_] = i + 1;
            denominator_exponent[number_of_real_roots_] = D;
            is_exact[number_of_real_roots_] = true;
            number_of_real_roots_++;
            return;
        }

        // Consider the first half of the interval.
        if(sign_variations(variation_transformation<Coeff__>(Q)) == 1) {
            strong_zero_one_descartes<Coeff__>(Q,i,D);
            return;
        }

        // Consider the second half of the interval.
        strong_zero_one_descartes<Coeff__>(S,i+1,D);
        return;
    }
};

} // namespace internal

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_D_DESCARTES_H
