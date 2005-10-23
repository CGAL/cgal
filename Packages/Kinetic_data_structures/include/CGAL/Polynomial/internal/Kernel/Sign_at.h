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

#ifndef CGAL_POLYNOMIAL_INTERNAL_SIGN_AT_ROOT_H
#define CGAL_POLYNOMIAL_INTERNAL_SIGN_AT_ROOT_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Explicit_root.h>
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/polynomial_converters.h>
#include <CGAL/Polynomial/internal/Rational/Sign_at_rational.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE
template <class Root, class K>
class Sign_at{
  typedef typename K::Function Poly;
public:
  Sign_at(const Poly &p, K k=K()): p_(p), k_(k) {
  }
  Sign_at(){}
  typedef typename K::Root argument_type;
  typedef CGAL_POLYNOMIAL_NS::Sign result_type;

  template <class T> 
  result_type operator()(const T &v) const {
    return eval(v);
  }
  /*result_type operator()(const typename K::NT &nt) const {
    return eval(nt);
    }*/
  
protected:

  template <class R>
  CGAL_POLYNOMIAL_NS::Sign eval(const R &r) const {
    std::pair<double, double> i= to_interval(r);
    if (i.first==i.second){
      double d= i.second;
      return eval(typename Poly::NT(d));
    } else {
      typename K::Root_stack s= k_.root_stack_object(p_, 
						     typename K::Root(i.first), 
						     typename K::Root(i.second));
      if (s.empty()){
	// there are no roots
	typename Poly::NT mid= (typename Poly::NT(i.first) 
				+ typename Poly::NT(i.second))*typename Poly::NT(.5);
	return eval(mid);
      } else {
	while (!s.empty() && s.top() < r){
	  s.pop();
	}
	if (!s.empty()){
	  if (s.top()==r){
	    return CGAL_POLYNOMIAL_NS::ZERO;
	  }
	  // now we know it is not a root
	  
	  typename K::Sign_between_roots sbr= k_.sign_between_roots_object(r, s.top());
	  
	  return sbr(p_);

	} else {
	  // There were roots below r.
	  typename K::Sign_between_roots sbr= k_.sign_between_roots_object(r, R(i.second));
	  return sbr(p_);

	}
	CGAL_postcondition(0);
	return CGAL_POLYNOMIAL_NS::ZERO;
      }
      //}
      //return sb;
    }
  }

  template <class RT>
  CGAL_POLYNOMIAL_NS::Sign eval(const CGAL_POLYNOMIAL_NS::Explicit_root<RT> &r){
    typedef  Explicit_root<RT> R;
    typename R::Representation rep= r.representation();
    typedef  typename CGAL_POLYNOMIAL_NS::Polynomial<typename R::Representation> Rep_poly;
    typename CGAL_POLYNOMIAL_NS::Polynomial_converter<typename K::Polynomial, Rep_poly> pc;
    return CGAL_POLYNOMIAL_NS::sign(pc(p_)(rep));
  }

  CGAL_POLYNOMIAL_NS::Sign eval(const typename Poly::NT &nt) const {
    typedef typename K::Root_stack_traits::Sign_at SA;
    SA sa= k_.root_stack_traits_object().sign_at_object(p_);
    return sa(nt);
  }

  Poly p_;
  K k_;
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
