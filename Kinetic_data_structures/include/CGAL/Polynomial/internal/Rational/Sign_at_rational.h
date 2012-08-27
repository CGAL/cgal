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

#ifndef CGAL_POLYNOMIAL_INTERNAL_SIGN_AT_H
#define CGAL_POLYNOMIAL_INTERNAL_SIGN_AT_H

#include <CGAL/Polynomial/basic.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {
template<class Fn>
class Sign_at_rational
{
protected:
  typedef CGAL::Sign  Sign;

public:
  typedef typename Fn::NT   NT;

  typedef Sign                 result_type;
  typedef Fn                   first_argument_type;
  typedef NT                   second_argument_type;

  Sign_at_rational(){}

  //! Evaluate the sign of the value of the polynomial at x
  template <class TNT>
  Sign operator()(const first_argument_type &p, const TNT& t) const
  {
    return CGAL::sign( p(t) );
  }

  //! Evaluate the sign of the value of the polynomial
  //  at n/d without any divisions
  Sign operator()(const first_argument_type &p, const NT &n, const NT& d) const
  {
    Sign s_n = CGAL::sign(n);
    Sign s_d = CGAL::sign(d);

    if ( s_d == CGAL::ZERO ) {
      CGAL_precondition( s_n != CGAL::ZERO );

      Sign s = CGAL::sign( p[p.degree()] );
      if ( CGAL::sign(n) == CGAL::POSITIVE ) {
	return s;
      }
      else {
	return opposite(s);
      }
    }

    Fn copy(p);
    if (p.is_zero()) {
      return CGAL::ZERO;
    }
    else if ( p.is_constant() == 0 || s_n == CGAL::ZERO ) {
      return CGAL::sign( copy[0] ) * s_d;
    }

    NT result = copy[copy.degree()];
    for (int i = copy.degree() - 1; i >= 0; i--) {
      result *= n;
      result += copy[i] * d;
      for (int j = i - 1; j >= 0; j--) {
	copy.set_coef(j, d * copy[j]);
      }
    }

    return CGAL::sign(result) * s_d;
  }
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
