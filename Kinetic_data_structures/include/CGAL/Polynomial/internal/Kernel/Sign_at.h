// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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

#ifndef CGAL_POLYNOMIAL_INTERNAL_SIGN_AT_ROOT_H
#define CGAL_POLYNOMIAL_INTERNAL_SIGN_AT_ROOT_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Explicit_root.h>
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/polynomial_converters.h>
#include <CGAL/Polynomial/internal/Rational/Sign_at_rational.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE
template <class Root, class K>
class Sign_at
{
  typedef typename K::Function Poly;
public:
  Sign_at(const Poly &p, K k=K()): p_(p), k_(k) {
  }
  Sign_at(){}
  typedef typename K::Root argument_type;
  typedef CGAL_POLYNOMIAL_NS::Sign result_type;

  template <class T>
  result_type operator()(const T &v) const
  {
    return eval(v);
  }
  /*result_type operator()(const typename K::NT &nt) const {
    return eval(nt);
    }*/

protected:

  template <class R>
  CGAL_POLYNOMIAL_NS::Sign eval(const R &r) const
  {
    typedef typename K::Root_stack_traits::Sign_at SA;
    SA sa= k_.root_stack_traits_object().sign_at_object(p_);
    return sa(r);
  }

  CGAL_POLYNOMIAL_NS::Sign eval(const typename K::Root &r) const {
    typename K::Is_rational ir= k_.is_rational_object();

    //std::pair<double, double> i= to_interval(r);
    if (ir(r)) {
      typename K::To_rational tr= k_.to_rational_object();
      typename K::NT nt= tr(r);
      return eval(nt);
    }
    else {
      typename K::To_isolating_interval tii= k_.to_isolating_interval_object();
      std::pair<typename K::NT, typename K::NT> ii= tii(r);
      typename K::Root_stack s= k_.root_stack_object(p_,
						     typename K::Root(ii.first),
						     typename K::Root(ii.second));
      if (s.empty()) {
	// there are no roots
	typename K::NT mid= (ii.first + ii.second)*typename K::NT(.5);
	return eval(mid);
      }
      else {
	while (!s.empty() && s.top() < r) {
	  s.pop();
	}
	if (!s.empty()) {
	  if (s.top()==r) {
	    return CGAL_POLYNOMIAL_NS::ZERO;
	  }
	  // now we know it is not a root

	  typename K::Sign_between_roots sbr= k_.sign_between_roots_object(r, s.top());

	  return sbr(p_);

	}
	else {
	  // There were roots below r.
	  typename K::Sign_between_roots sbr= k_.sign_between_roots_object(r, typename K::Root(ii.second));
	  return sbr(p_);

	}
      }
      //}
      //return sb;
    }
    CGAL_postcondition(false);
    return CGAL_POLYNOMIAL_NS::ZERO;
  }

  template <class RT>
  CGAL_POLYNOMIAL_NS::Sign eval(const CGAL_POLYNOMIAL_NS::internal::Explicit_root<RT> &r) {
    typedef  internal::Explicit_root<RT> R;
    typename R::Representation rep= r.representation();
    typedef  typename CGAL_POLYNOMIAL_NS::Polynomial<typename R::Representation> Rep_poly;
    typename CGAL_POLYNOMIAL_NS::Polynomial_converter<typename K::Polynomial, Rep_poly> pc;
    return CGAL_POLYNOMIAL_NS::sign(pc(p_)(rep));
  }

  /* CGAL_POLYNOMIAL_NS::Sign eval(const typename Poly::NT &nt) const
  {
   
  }*/

  Poly p_;
  K k_;
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
