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

#ifndef CGAL_POLYNOMIAL_INTERNAL_AFTER_AT_ROOT_H
#define CGAL_POLYNOMIAL_INTERNAL_AFTER_AT_ROOT_H

#include <CGAL/Polynomial/basic.h>
//#include <CGAL/Polynomial/internal/Explicit_root.h>
#include <CGAL/Polynomial/internal/Rational/Evaluate_polynomial.h>
#include <CGAL/Polynomial/internal/Rational/Sign_above_rational.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

//! Compute the sign after a root.
/*!
  This has specializations for Explicit_roots.
*/
template <class RNT, class K>
class Sign_above: Sign_above_rational<K>
{
  typedef typename K::Function Poly;
  typedef Sign_above_rational<K> P;
public:
  Sign_above(const Poly &p, K k):P(p,k) {
  }
  Sign_above(){}

  using P::operator();

  typename P::result_type operator()(const typename K::Root &v) const
  {
    CGAL_Polynomial_expensive_precondition(k_.sign_at_root_object(p_)(v)==CGAL::ZERO);
    return eval(v);
  }
protected:

  template <class R>
  CGAL::Sign eval(const R &r) const
  {
    double ub= CGAL::to_interval(r).second+.00001;
    if (ub== CGAL::to_interval(r).second) ub= ub*2;
    typename K::Root_stack rc=  k_.root_stack_object(p_, r, ub);
    if (rc.empty()) {
      return CGAL::sign(p_(typename Poly::NT(ub)));
    } else {
      typename K::Root rr=rc.top();
      typename K::Sign_between_roots sb= k_.sign_between_roots_object(r, rr);
      return sb(p_);
    }
  }

#if 0
  template <class RT>
  CGAL::Sign eval(const Explicit_root<RT> &r) {
    typedef  Explicit_root<RT> R;
    typename R::Representation rep= r.representation();
    CGAL_POLYNOMIAL_NS::Polynomial_converter<typename K::Polynomial, Polynomial<typename R::Representation> > pc;
    CGAL::Sign sn=  CGAL::sign(pc(p_)(rep));
    if (sn==CGAL::ZERO) {
      typename K::Root rr= k_.solver_object(p_, r).next_root();
      typename R::Representation mr= .5*(rr+r);
      return CGAL::sign(pc(p_)(mr));
    }
    else {
      return sn;
    }
  }
#endif
  Poly p_;
  K k_;
};

//! An explicit instantiation with Explicit_root
/*!
  This evaluates the polynomials at the root and takes the sign. Much faster for doubles at least.
*/
template <class RNT, class Kernel>
class Sign_above<CGAL_POLYNOMIAL_NS::internal::Explicit_root<RNT>, Kernel>
{
  typedef Evaluate_polynomial<Kernel, RNT> Eval;
public:
  Sign_above(const typename Kernel::Function &fn, Kernel ) {
    fns_.push_back(Eval(fn));
    last_= fn;
  }
  Sign_above(){}

  typedef CGAL_POLYNOMIAL_NS::internal::Explicit_root<RNT> argument_type;
  typedef CGAL::Sign result_type;

  result_type operator()(const argument_type &v) const
  {
    RNT rep= v.representation();

    CGAL::Sign sn;
    unsigned int i=0;
    while ((sn=tsign(i, rep)) == CGAL::ZERO) {
      ++i;
    }
    return sn;
  }

protected:
  CGAL::Sign tsign(unsigned int i, const RNT &rnt) const
  {
    while  (i >= fns_.size()) {
      last_= last_.derivative();
      fns_.push_back(Eval(last_));;
    }
    return CGAL::sign(fns_[i](rnt));
  }

  mutable std::vector<Eval> fns_;
  mutable typename Kernel::Function last_;
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
