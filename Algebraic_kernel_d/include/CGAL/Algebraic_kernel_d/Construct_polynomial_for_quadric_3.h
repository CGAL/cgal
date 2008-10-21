// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $URL$
// $Id$
//
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_D_CONSTRUCT_POLYNOMIAL_FOR_QUADRIC_3
#define CGAL_ALGEBRAIC_KERNEL_D_CONSTRUCT_POLYNOMIAL_FOR_QUADRIC_3 1

/*!\file include/CGAL/Algebraic_kernel_d/Construct_polynomial_for_quadric_3.h
 * \brief Class that defines functor
 * Construct_polynomial_for_quadric_3
 */

#include <CGAL/config.h>

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>

CGAL_BEGIN_NAMESPACE

/*!\brief Offers different possibilities to construct the defining
 * polynomial for a quadric (canonicalized version)
 */
template < typename Coefficient_ >
class Construct_polynomial_for_quadric_3 {

public:
    //! type of coefficient
    typedef Coefficient_ Coefficient;

    //! type of univariate polynonomial
    typedef CGAL::Polynomial< Coefficient > Poly_coeff_1;
    //! type of univariate polynonomial
    typedef CGAL::Polynomial< Poly_coeff_1 > Poly_coeff_2;
    //! type of univariate polynonomial
    typedef CGAL::Polynomial< Poly_coeff_2 > Poly_coeff_3;

    //! type of result
    typedef Poly_coeff_3 result_type;

private:

    /*! \brief
     * constructs the trivariate polynomial of a quadric from ten coefficients
     *
     * \f[ f = Ax^2 + Bxy + Cxz + Dy^2 + Eyz + Fz^2 + Gx + Hy + Kz + L \f]
     */
    Poly_coeff_3 _polynomial_from_coefficients(
            const Coefficient &a, 
            const Coefficient &b, 
            const Coefficient &c, 
            const Coefficient &d, 
            const Coefficient &e, 
            const Coefficient &f, 
            const Coefficient &g, 
            const Coefficient &h, 
            const Coefficient &k, 
            const Coefficient &l) {
        
        Poly_coeff_1 p00(l, g, a);
        Poly_coeff_1 p01(h, b);
        Poly_coeff_1 p02(d);
        
        Poly_coeff_1 p10(k,c);
        Poly_coeff_1 p11(e);
        
        Poly_coeff_1 p20(f);
        
        Poly_coeff_2 p0(p00,p01,p02);
        Poly_coeff_2 p1(p10,p11);
        Poly_coeff_2 p2(p20);
        
        Poly_coeff_3 ret(p0,p1,p2);
        
        return ret;
    }
    
public:
    
    //!\brief return a_0z^0 + a_1z^1 + a_2^2
    result_type operator()(const Poly_coeff_2& a_0, 
                           const Poly_coeff_2& a_1, 
                           const Poly_coeff_2& a_2) {
        CGAL_precondition(CGAL::total_degree(a_0) <= 2);
        CGAL_precondition(CGAL::total_degree(a_0) <= 1);
        CGAL_precondition(CGAL::total_degree(a_0) == 0);
        Poly_coeff_3 poly(a_0,a_1,a_2);
        // TODO replace with non-CGALi-version
        Poly_coeff_3 ret = CGAL::CGALi::canonicalize_polynomial(poly);
        return ret;
    }

    /*!\brief constructs from 10 coefficients.
     *
     * \f[ f = Ax^2 + Bxy + Cxz + Dy^2 + Eyz + Fz^2 + Gx + Hy + Kz + L \f]
     */
    result_type operator()(const Coefficient& a, 
                           const Coefficient& b, 
                           const Coefficient& c,
                           const Coefficient& d, 
                           const Coefficient& e, 
                           const Coefficient& f,
                           const Coefficient& g, 
                           const Coefficient& h, 
                           const Coefficient& k,
                           const Coefficient& l) {
        Poly_coeff_3 poly = 
            _polynomial_from_coefficients(a,b,c,d,e,f,g,h,k,l);
        // TODO replace with non-CGALi-version
        Poly_coeff_3 ret = CGAL::CGALi::canonicalize_polynomial(poly);
        return ret;
    }
    
    // TODO from string
    
    // TODO from stream
    
};

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_D_CONSTRUCT_POLYNOMIAL_FOR_QUADRIC_3
// EOF
