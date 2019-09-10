// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// =============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_SCALAR_FACTOR_TRAITS_H
#define CGAL_SCALAR_FACTOR_TRAITS_H

#include <CGAL/Algebraic_structure_traits.h>

/*! \file NiX/Scalar_factor_traits.h
 *  \brief Definition of traits class NiX::Scalar_factor_traits.
 */
namespace CGAL {

/*! \nosubgrouping
 *  \ingroup NiX_Scalar_factor_traits_spec
 *  \brief traits class for scalar factor extraction
 *  
 *  A model of the concept ScalarFactorTraits. 
 * 
 *  This is the default implementation of NiX::Scalar_factor_traits. 
 *  It is valid for most basic number types as Integer, Rational etc. 
 */
template< class NT >
class Scalar_factor_traits {
public:
    // the supported number type
    typedef NT Type;
    // NT is also 
    typedef NT Scalar;

    // functor computing the scalar factor of a basic number type. 
    class Scalar_factor {
        static Scalar scalar_factor(const NT& a, Integral_domain_tag) {
            typename Algebraic_structure_traits<NT>::Unit_part upart;
            typename Algebraic_structure_traits<NT>::Integral_division idiv;
            return idiv(a, upart(a));
        }
        static Scalar scalar_factor(const NT& a, Field_tag) {
            return (a == NT(0)) ? NT(0) : NT(1);
        }
        static Scalar scalar_factor(const NT& /* a */, 
                                    const Scalar& /* d */,
                                    Integral_domain_tag) {
            return Scalar(1);
        }
        static Scalar scalar_factor(const NT& a, 
                                    const Scalar& d,
                                    Field_tag) {
            return (a == NT(0) && d == Scalar(0)) ? NT(0) : NT(1);
        }
        static Scalar scalar_factor(const NT& a, 
                                    const Scalar& d,
                                    Unique_factorization_domain_tag) { 
            if( d == Scalar(0)) {
                return a;
            } else {
                typename Algebraic_structure_traits<NT>::Gcd gcd;
                return gcd(a,d);
            }
        }
    public: 
        // argument type
        typedef NT argument_type;
        // first argument type
        typedef NT first_argument_type;
        // second argument type
        typedef Scalar second_argument_type;
        // result type
        typedef Scalar result_type;
  
        // determine extractable scalar factor
        Scalar operator () (const NT& a) {
            CGAL_static_assertion(( ::boost::is_same< NT,Scalar >::value));  
            typedef typename Algebraic_structure_traits<NT>::Algebraic_category SAT;
            return scalar_factor(a, SAT());
        }
        // determine extractable scalar factor
        Scalar operator () (const NT& a, const Scalar& d) {     
            CGAL_static_assertion(( ::boost::is_same< NT,Scalar >::value));  
            typedef typename Algebraic_structure_traits<NT>::Algebraic_category SAT;
            return scalar_factor(a,d,SAT());
        }
    };

    class Scalar_div {
    public:
        // first argument type
        typedef NT& first_argument_type;
        // second argument type
        typedef Scalar second_argument_type;
        // result type
        typedef void result_type;
        // divide \a a in place by scalar factor \a b
        void operator () (NT& a, const Scalar& b) {
            CGAL_precondition(b != Scalar(0));
            typename Algebraic_structure_traits<NT>::Integral_division idiv;
            a = idiv(a, b);
        }
    };
};

/*! \ingroup NiX_Scalar_factor_traits
 *  \brief removes the scalar factor
 * 
 *   Removes the scalar factor provided by NiX::Scalar_factor_traits
 */
template <class NT>
void remove_scalar_factor(NT& x){
    typedef CGAL::Scalar_factor_traits<NT> SFT;
    typedef typename SFT::Scalar Scalar;
    typename SFT::Scalar_factor sfac;
    typename SFT::Scalar_div    sdiv;
    
    Scalar scalar=sfac(x);
    if(scalar != Scalar(0))
        sdiv(x,scalar);
}

template <class NT>
inline
typename CGAL::Scalar_factor_traits<NT>::Scalar
scalar_factor(const NT& x){
    typename CGAL::Scalar_factor_traits<NT>::Scalar_factor scalar_factor;
	return scalar_factor(x);
}

template <class NT>
inline
typename CGAL::Scalar_factor_traits<NT>::Scalar 
scalar_factor(const NT& x,
        const typename CGAL::Scalar_factor_traits<NT>::Scalar& d)
{
    typename CGAL::Scalar_factor_traits<NT>::Scalar_factor scalar_factor;
    return scalar_factor(x,d);
}


} //namespace CGAL

#endif // NiX_SCALAR_FACTOR_TRAITS_H
// EOF
