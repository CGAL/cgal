// Copyright (c) 2005  Stanford University (USA).
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
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_INTERNAL_QUOTIENT_REMAINDER_H
#define CGAL_POLYNOMIAL_INTERNAL_QUOTIENT_REMAINDER_H

#include <CGAL/Polynomial/basic.h>

/*!
  \file Quotient_remainder.h A class that computes quotient and remainder.
*/

namespace CGAL { namespace POLYNOMIAL { namespace internal {

//! Compute the quotient and remainder of two polynomials.
template<class Polynomial>
struct Quotient_remainder
{
    typedef typename Polynomial::NT   NT;
    typedef std::pair<Polynomial,Polynomial>  result_type;
    typedef Polynomial       argument_type;
    typedef Polynomial       argument_type1;
    typedef Polynomial       argument_type2;

    void write(std::ostream &out) const
    {
        out << "quo_rem";
    }

//! compute the quotient and remainder
    result_type
        operator()(const Polynomial& t, const Polynomial& v) const
    {
        int n = v.degree();
        int m = t.degree();
        int divdeg = m - n;

        std::vector<NT> rem_coef(m + 1);
        std::vector<NT> div_coef(divdeg + 1);

        for (int i = 0; i <= m; i++) {
            rem_coef[i] = t[i];
        }

        for (int k = divdeg; k >= 0; k--) {
            div_coef[k] = rem_coef[n + k] / v[n];
            for (int j = n + k - 1; j >= k; j--) {
                rem_coef[j] = rem_coef[j] - div_coef[k] * v[j - k];
            }
        }

        rem_coef.resize(n);
#if 0
        while ( rem.degree() > n - 1 ) {
            rem.set_leading_coef_to_zero();
        }
#endif
        Polynomial div(div_coef.begin(), div_coef.end());
        Polynomial rem(rem_coef.begin(), rem_coef.end());

        return result_type(div,rem);
    }
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif                                            // CGAL_POLYNOMIAL_INTERNAL_QUOTIENT_REMAINDER_H
