// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
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
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Michael Hemmer <hemmer@informatik.uni-mainz.de> 
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.


#ifndef CGAL_POLYNOMIAL_SCALAR_FACTOR_TRAITS_H
#define CGAL_POLYNOMIAL_SCALAR_FACTOR_TRAITS_H

#include <CGAL/basic.h>

namespace CGAL {

// see <CGAL/Scalar_factor_traits.h>
/*! \ingroup CGAL_Polynomial
 *  \ingroup CGAL_Scalar_factor_traits_spec
    \brief \c CGAL::Scalar_factor_traits < \c CGAL::Polynomial<NT> >
 *
 *  The \c CGAL::Scalar_factor_traits::Scalar_factor of a polynomial,
 *  even of a nested polynomial, is essentially the gcd of its scalar
 *  coefficients. (This differs from the \c content() in that the content
 *  is the gcd of the coefficient sequence, which is again a polynomial,
 *  albeit with one variable less, for nested polynomials.) The scalar
 *  factor of the zero polynomial is 0.
 *
 *  This currently only works if the scalar type is a model of \c UFDomain.
 */
template <class Coeff>
class Scalar_factor_traits< Polynomial<Coeff> > {
public:
    typedef Polynomial<Coeff> NT;
    typedef typename Scalar_factor_traits<Coeff>::Scalar Scalar;
    class Scalar_factor {
    public:
        //! argument type
        typedef NT argument_type;
        //! first argument type
        typedef NT first_argument_type;
        //! second argument type
        typedef Scalar second_argument_type;
        //! result type
        typedef Scalar result_type;
        //! returns the gcd of a polynomials coefficients
        Scalar operator()(const NT& p, const Scalar& d_=Scalar(0)) const {
            
            typename Scalar_factor_traits<Coeff>::Scalar_factor sfac;
            const Scalar unity(1);
            
            Scalar d(d_);
            if (p.is_zero()) return d;
            
            int i = p.degree();
            while((d != unity) && (i >= 0)) {
                d = sfac(p[i--],d);
            }
            return d;
        }
    };
    class Scalar_div {
    public:
        //! first argument type
        typedef NT first_argument_type;
         //! second argument type
        typedef Scalar second_argument_type;
        //! divides a polynomial \c p by a scalar factor \c b
        void operator () (NT& p, const Scalar& b) const { 
            CGAL_precondition(b != Scalar(0));
            p.scalar_div(b); 
        }
    };
};

} //namespace CGAL
#endif // CGAL_POLYNOMIAL_SCALAR_FACTOR_TRAITS_H
