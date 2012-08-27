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

#ifndef CGAL_POLYNOMIAL_INTERNAL_REMAINDER_H
#define CGAL_POLYNOMIAL_INTERNAL_REMAINDER_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Rational/Quotient_remainder.h>
/*!
  \file Remainder.h A class to compute remainders.
*/

namespace CGAL { namespace POLYNOMIAL { namespace internal {

//! Compute the remainder of two polynomials.
/*!
  I pulled this out of Polynomial because I did not think polynomial should have such complicated methods.
*/
template<class Polynomial>
class Remainder : private Quotient_remainder<Polynomial>
{
    private:
        typedef Quotient_remainder<Polynomial>  Base;

    public:
        typedef typename Polynomial::NT   NT;
        typedef Polynomial       result_type;
        typedef Polynomial       argument_type;
        typedef Polynomial       argument_type1;
        typedef Polynomial       argument_type2;

        void write(std::ostream &out) const
        {
            out << "rem";
        }
//! compute the remainder
        result_type
            operator()(const Polynomial& t, const Polynomial& v) const
        {
            return Base::operator()(t,v).second;
        }
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif                                            // CGAL_POLYNOMIAL_INTERNAL_REMAINDER_H
