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

#ifndef CGAL_POLYNOMIAL_INTERNAL_UPPER_BOUND_ENUMERATOR_TRAITS_BASE_H
#define CGAL_POLYNOMIAL_INTERNAL_UPPER_BOUND_ENUMERATOR_TRAITS_BASE_H

#include <CGAL/Polynomial/basic.h>

#include <CGAL/Polynomial/internal/Isolating_interval.h>
#include <CGAL/Polynomial/internal/Rational/Rational_traits_base.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <class Poly>
class Root_stack_traits_base: public Rational_traits_base<Poly>
{
    private:
        typedef Root_stack_traits_base<Poly> This;
        typedef Rational_traits_base<Poly> P;

    public:
  typedef CGAL_POLYNOMIAL_NS::internal::Isolating_interval<typename P::FT> Isolating_interval;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
