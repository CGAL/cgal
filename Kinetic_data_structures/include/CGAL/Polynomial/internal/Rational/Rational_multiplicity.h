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

#ifndef CGAL_POLYNOMIAL_INTERNAL_ROOT_MULTIPLICITY_H
#define CGAL_POLYNOMIAL_INTERNAL_ROOT_MULTIPLICITY_H

#include <CGAL/Polynomial/basic.h>
#include <vector>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <class K>
class Rational_multiplicity
{
public:
  Rational_multiplicity(){}

  Rational_multiplicity(const K &k):d_(k.differentiate_object()),
				   k_(k) {
  }

  typedef unsigned int result_type;
  //typedef Bound_type argument_type;
  typedef typename K::Function first_argument_type;
  typedef typename K::Function::NT second_argument_type;

  result_type operator()(const first_argument_type &f, const second_argument_type &t) const
  {
    CGAL_precondition(f.degree() > -1);
    typename K::Function h=f;
    //CGAL_Polynomial_exactness_assertion(k_.sign_at_object(h)( t)== CGAL::ZERO);
    //POLYNOMIAL_NS::Sign sn;
    //if ( k.sign_at_object(fh)(t) != POLYNOMIAL_NS::ZERO ) return 0;
    unsigned int deg=0;
    
    typename K::Sign_at sa= k_.sign_at_object();
    
    //unsigned int mdegree= h.degree();
    while (sa(h, t) ==CGAL::ZERO) {
      ++deg;
      h= d_(h);
    }
    return deg;
  }

protected:

  typename K::Differentiate d_;
  K k_;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
