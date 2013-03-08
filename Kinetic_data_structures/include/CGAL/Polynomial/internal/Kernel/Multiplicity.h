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

#ifndef CGAL_POLYNOMIAL_INTERNAL_MULTIPLICITY_H
#define CGAL_POLYNOMIAL_INTERNAL_MULTIPLICITY_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Rational/Rational_multiplicity.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

//! Compute the sign after a root.
/*!
  This has specializations for Explicit_roots.
*/
template <class K>
class Multiplicity: public Rational_multiplicity<K>
{
    typedef Rational_multiplicity<K> P;
    public:
        Multiplicity(const typename K::Function &p, K k=K()): P(p,k){  }
        Multiplicity(){}

        using P::operator();
        typename P::result_type operator()(const typename K::Root &) const
        {
            CGAL_Polynomial_precondition(0);
            return 1;
        }
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
