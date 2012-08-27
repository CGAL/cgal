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

#ifndef CGAL_POLYNOMIAL_INTERNAL_SHIFT_POWER_H
#define CGAL_POLYNOMIAL_INTERNAL_SHIFT_POWER_H

#include <CGAL/Polynomial/basic.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

//------------------------------------------------------------------

template <class Polynomial>
class Shift_power
{
    public:
        typedef typename Polynomial::NT NT;
        typedef unsigned int argument_type;
        typedef Polynomial result_type;
        Shift_power(){}
        Shift_power(unsigned int p): p(p){}

        void write(std::ostream &out) const
        {
            out << "shift_power";
        }

        Polynomial operator()(const Polynomial& f) const
        {
#if 0
// this is the old buggy code: the for-loop below was accessing
// an array element which was out of the array bounds
            int deg = f.degree();
            unsigned int sdeg = deg + p + 1;

            std::vector<NT> ret_coef(sdeg);

            for (unsigned int i = 0; i < p; i++) {
                ret_coef[i] = NT(0);
            }
            for (unsigned int i = p; i <= sdeg; i++) {
                ret_coef[i] = f[i-p];
            }

            Polynomial ret(ret_coef.begin(), ret_coef.end());

            return ret;
#endif
            if ( p == 0 ) { return f; }

            int new_deg = f.degree() + p;

            if ( new_deg < 0 ) { return Polynomial(); }

            unsigned int unew_deg = new_deg;

            std::vector<NT> ret_coef(unew_deg+1);

            for (unsigned int i = 0; i < p; i++) {
                ret_coef[i] = NT(0);
            }
            for (unsigned int i = p; i <= unew_deg; i++) {
                ret_coef[i] = f[i-p];
            }

            return Polynomial(ret_coef.begin(), ret_coef.end());
        }
    protected:
        unsigned int p;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
