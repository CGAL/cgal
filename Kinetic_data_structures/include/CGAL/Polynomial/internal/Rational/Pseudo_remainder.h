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

#ifndef CGAL_POLYNOMIAL_INTERNAL_PSEUDO_REMAINDER_H
#define CGAL_POLYNOMIAL_INTERNAL_PSEUDO_REMAINDER_H

#include <CGAL/Polynomial/basic.h>

/*!
  \file Pseudo_remainder.h A class to compute pseudo remainders.
*/

namespace CGAL { namespace POLYNOMIAL { namespace internal {

//! Compute the Pseudo remainder of two polynomials.
/*!
  I pulled this out of Polynomial because I did not think polynomial should have such complicated methods.
*/
template<class Polynomial>
struct Pseudo_remainder
{
    typedef typename Polynomial::NT   NT;
    typedef Polynomial       result_type;
    typedef Polynomial       argument_type;
    typedef Polynomial       argument_type1;
    typedef Polynomial       argument_type2;

    void write(std::ostream &out) const
    {
        out << "prem";
    }

//! compute the pseudo remainder
/*!
  \todo Why do we need the !STRIP_ZEROS line?
*/
    Polynomial
        operator()(const Polynomial& t, const Polynomial& v) const
    {
        CGAL_Polynomial_precondition( t.degree() >= v.degree() );

        int n = v.degree();
        int m = t.degree();
        int divdeg = m - n;
        CGAL_Polynomial_assertion( divdeg >= 0 );

        std::vector<NT> r_coef(m+1);
        std::vector<NT> q_coef(divdeg+1);

        for (int i = 0; i <= m; i++) {
            r_coef[i] = t[i];
        }

//    Polynomial r(t);
//    Polynomial q;

//    q.set_nominal_degree(divdeg);

        for (int k = divdeg; k >= 0; k--) {
            q_coef[k] = r_coef[n + k];
            for (int j = n + k - 1; j >= k; j--) {
                r_coef[j] = v[n] * r_coef[j] - q_coef[k] * v[j - k];
            }
            for (int j = k - 1; j >= 0; j--) {
                r_coef[j] = r_coef[j] * v[n];
            }
            q_coef[k] = q_coef[k] * v[n];
        }

        r_coef.resize(n);

        Polynomial r(r_coef.begin(), r_coef.end());

/*!
  \todo Why did we need this? Negation should be safe.
*/
/*if ( !STRIP_ZEROS ) {
  if ( divdeg % 2 == 0 ) {
r = r * v[n];
  }
  return r;
  }*/

        CGAL::Sign s_vn = CGAL::sign(v[n]);

        if ( (divdeg % 2 == 0) && s_vn == CGAL::NEGATIVE ) {
            r = -r;
        }

        return r;
    }
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
