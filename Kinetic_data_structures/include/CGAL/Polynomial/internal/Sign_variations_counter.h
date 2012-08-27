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

#ifndef CGAL_POLYNOMIAL_SIGN_VARIATIONS_COUNTER_H
#define CGAL_POLYNOMIAL_SIGN_VARIATIONS_COUNTER_H

#include <CGAL/Polynomial/basic.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

struct Sign_variations_counter
{
    template<class Iterator>
        static
        unsigned int sign_variations(const Iterator& first,
        const Iterator& beyond,
    bool stop_if_more_than_one = false) {
        Iterator it1 = first;
        Iterator last = beyond;
        last--;

        int count(0);

        while ( it1 != beyond && (int)(*it1) == (int)CGAL::ZERO ) { it1++; }

        while ( it1 != last && it1 != beyond ) {
            Iterator it2 = it1;
            it2++;
            while ( it2 != beyond && (int)(*it2) == (int)CGAL::ZERO ) { it2++; }

            if ( it2 != beyond && (int)(*it1) != (int)(*it2) ) {
                if ( stop_if_more_than_one && count > 1 ) { break; }
                count++;
            }
            it1 = it2;
        }

        return count;
    }

    template<class Iterator>
        unsigned int operator()(const Iterator& first,
        const Iterator& beyond,
        bool stop_if_more_than_one = false) const
    {
        return sign_variations(first, beyond, stop_if_more_than_one);
    }

};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif                                            // CGAL_POLYNOMIAL_SIGN_VARIATIONS_COUNTER_H
