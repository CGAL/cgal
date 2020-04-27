// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Michael Hemmer <hemmer@informatik.uni-mainz.de>
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.


#ifndef CGAL_POLYNOMIAL_FRACTION_TRAITS_H
#define CGAL_POLYNOMIAL_FRACTION_TRAITS_H

#include <CGAL/basic.h>

namespace CGAL {

// We need to play a similar game to provide Fraction_traits

template <class POLY, class TAG>
class Poly_Ftr_base;

// Use this if the coefficients cannot be decomposed
// into numerator and denominator
template <class NT_>
class Poly_Ftr_base< Polynomial<NT_>, CGAL::Tag_false > {
public:
    typedef Polynomial<NT_> Type;
    typedef CGAL::Tag_false Is_fraction;
    typedef CGAL::Null_tag Numerator;
    typedef CGAL::Null_tag Denominator_type;
    typedef CGAL::Null_functor Common_factor;
    typedef CGAL::Null_functor Decompose;
    typedef CGAL::Null_functor Compose;
};

// If they can, use this
template <class NT_>
class Poly_Ftr_base< Polynomial<NT_>, CGAL::Tag_true > {
    typedef Polynomial<NT_> Poly;
    typedef NT_ Coefficient_type;
public:
    typedef Polynomial<NT_> Type;
    typedef CGAL::Tag_true Is_fraction;
    typedef Polynomial<typename Fraction_traits<NT_>::Numerator_type>
        Numerator_type;
    typedef typename Fraction_traits<NT_>::Denominator_type Denominator_type;
    typedef typename Fraction_traits<NT_>::Common_factor Common_factor;
    class Decompose {
    public:
        typedef Type first_argument_type;
        typedef Numerator_type& second_argument_type;
        typedef Denominator_type& third_argument_type;
        inline void operator () (
                const Type& p,
                Numerator_type& num,
                Denominator_type& den){

            typedef Numerator_type INTPOLY;
            typedef Denominator_type DENOM;

            typedef Fraction_traits<Coefficient_type> CFTRAITS;
            typedef typename CFTRAITS::Numerator_type INTCOEFF;

            const int d = p.degree();
            std::vector<INTCOEFF> integ(d+1);
            std::vector<DENOM> denom(d+1);

            int i;

            // decompose each coefficient into integral part and denominator
            typename CFTRAITS::Decompose decomp_coeff;
            for (i = 0; i <= d; i++) {
                decomp_coeff(p[i], integ[i], denom[i]);
            }

            // c = lcm(denom[0], ..., denom[d])
            typename Algebraic_structure_traits<DENOM>::Integral_division idiv;
            typename CFTRAITS::Common_factor  gcd;  // not really `greatest'

            den = denom[0];
            for (i = 1; i <= d; i++) {
                den *= idiv(denom[i], gcd(den, denom[i]));
            }

            // expand each (integ, denom) pair to common denominator
            for (i = 0; i <= d; i++) {
                integ[i] *= INTCOEFF(idiv(den, denom[i]));
            }
            num =  INTPOLY(integ.begin(), integ.end());
        }
    };

    class Compose {
    public:
        typedef Numerator_type first_argument_type;
        typedef Denominator_type second_argument_type;
        typedef Type result_type;
        inline Type operator () (const Numerator_type& n,
                                 const Denominator_type& d){
            typename Fraction_traits<NT_>::Compose comp_coeff;
            (void)comp_coeff;

            std::vector< NT_> coeffs(n.degree()+1);

            for (int i = 0; i <= n.degree(); i++) {
                coeffs[i] = comp_coeff(n[i], d);
            }

            return Type(coeffs.begin(), coeffs.end());
        };
    };
};


// Select the right alternative as Fraction_traits
/*! \ingroup CGAL_Polynomial
    \brief \c CGAL::Fraction_traits < \c CGAL::Polynomial<NT> >
 *
 *  Polynomials provide suitable specializations of \c CGAL::Fraction_traits.
 *  They are decomposable iff their coefficient type is.
 *  The denominator \e d of a polynomial \e p is a low common multiple
 *  (see \c CGAL::Fraction_traits::Common_factor for details) of the
 *  denominators of its coefficients.  The numerator is the polynomial
 *  \e d*p with a fraction-free coefficient type.
 *
 *  This works for nested polynomials, too.
 */
template <class NT_>
class Fraction_traits< Polynomial<NT_> >
    : public Poly_Ftr_base< Polynomial<NT_>,
                 typename Fraction_traits<NT_>::Is_fraction >
{
    // nothing new
};

} //namespace CGAL
#endif // CGAL_POLYNOMIAL_FRACTION_TRAITS_H
