// Copyright (c) 2002  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// Author(s)     : Andreas Fabri, Susan Hert, Sylvain Pion

#ifndef CGAL_TO_RATIONAL_H
#define CGAL_TO_RATIONAL_H

#include <CGAL/number_type_basic.h>
#include <CGAL/Fraction_traits.h>

#include <cmath>

namespace CGAL {

template <class Rational>
Rational
to_rational(double x)
{
    typedef Fraction_traits<Rational> FT;
    typedef typename FT::Is_fraction Is_fraction;
    typedef typename FT::Numerator_type Numerator_type;
    typedef typename FT::Denominator_type Denominator_type;
    typename FT::Compose compose;

    CGAL_static_assertion((::boost::is_same<Is_fraction,Tag_true>::value));
    CGAL_USE_TYPE(Is_fraction);
    CGAL_static_assertion((::boost::is_same<Numerator_type,Denominator_type>::value));
    CGAL_USE_TYPE(Denominator_type);

    Numerator_type num(0),den(1);

    if (x != 0.0)
        { bool neg = (x < 0);
            if (neg) x = -x;

            const unsigned shift = 15;   // a safe shift per step
            const int shift_pow = 32768; // = 2^shift
            const double width = 32768;  // = 2^shift
            const int maxiter = 20;      // ought not be necessary, but just in
            // case, max 300 bits of precision
            int expt;
            double mantissa = std::frexp(x, &expt);
            long exponent = expt;
            double intpart;
            int k = 0;

            while (mantissa != 0.0 && k++ < maxiter)
                {
                    mantissa *= width; // shift double mantissa
                    mantissa = std::modf(mantissa, &intpart);
                    num *= shift_pow;
                    num += (int)intpart;
                    exponent -= shift;
                }
            int expsign = (exponent>0 ? +1 : (exponent<0 ? -1 : 0));
            exponent *= expsign;
            Numerator_type twopot(2);
            Numerator_type exppot(1);
            while (exponent!=0) {
                if (exponent & 1)
                    exppot *= twopot;
                exponent >>= 1;
                twopot *= twopot;
            }
            if (expsign > 0)
                num *= exppot;
            else if (expsign < 0)
                den *= exppot;
            if (neg)
                num = -num;
        }
    return compose(num,den);
}

} //namespace CGAL

#endif // CGAL_TO_RATIONAL_H
