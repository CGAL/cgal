// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_INTERNAL_BELOW_AT_ROOT_H
#define CGAL_POLYNOMIAL_INTERNAL_BELOW_AT_ROOT_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Rational/Sign_below_rational.h>


CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

//! Compute the sign after a root.
/*!
  This has specializations for Explicit_roots. 
*/
template <class R, class K>
class Sign_below: public Sign_below_rational<K>{
  typedef Sign_below_rational<K> P;
public:
  Sign_below(const typename K::Function &p, K k): P(p,k) {
  }
  Sign_below(){}
  using P::operator();
  typename P::result_type operator()(const typename K::Root &v) const {
    CGAL_Polynomial_precondition(0);
    return ZERO;
  }

};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
