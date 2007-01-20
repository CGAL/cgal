// Copyright(c) 2005  Stanford University (USA).
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

#ifndef CGAL_POLYNOMIAL_KERNEL_RATIONAL_BETWEEN_ROOTS_H
#define CGAL_POLYNOMIAL_KERNEL_RATIONAL_BETWEEN_ROOTS_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/interval_arithmetic.h>
#include <CGAL/Polynomial/internal/nt_converters.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template <class K>
struct Rational_between_roots
{
  typedef typename K::FT result_type;
  typedef typename K::Root first_argument_type;
  typedef typename K::Root second_argument_type;

  Rational_between_roots(const K& k): tii_(k.to_isolating_interval_object()){}

  template <class T>
  static bool is_minf(T r){
    if (std::numeric_limits<result_type>::has_infinity) return r== -std::numeric_limits<result_type>::infinity();
    else return false;
  }


protected:
  template <class Rt, class Result> 
  Result compute(const Rt &r0, const Rt &r1)const  {
  CGAL_exactness_precondition(r0<r1);
    if (r0 > r1) return operator()(r1, r0);
    typedef std::pair<Result, Result> Ival;
    //Ival i0= to_interval(r0);
    //Ival i1= to_interval(r1);
    Ival i0= tii_(r0);
    if (is_minf(i0.second)) {
      result_type ret(-33554432);
      while (r1 <= second_argument_type(ret)) {
	ret= ret*1024;
      }
      //std::cout << "Returning "<< ret << " between " << r0 << " and " << r1 << std::endl;
      return ret;
    }

    if (i0.second == i0.first) {
      result_type ret= i0.second;
      //2^-24
      result_type incr=.0000000596046447753906250000000;
      CGAL_Polynomial_postcondition(incr!=0);
      while (second_argument_type(ret+incr) >= r1) {
	incr=incr/2;
      }
      //std::cout << "Returning "<< ret+incr << " between " << r0 << " and " << r1 << std::endl;
      return ret+incr;
    }
    else if (second_argument_type(i0.second) < r1) {
      //std::cout << "Returning "<< i0.second << " between " << r0 << " and " << r1 << std::endl;
      return Result(i0.second);
    }
    else {
      result_type lb = i0.first;
      result_type ub = i0.second;
      result_type mid;
      while (true) {
	mid = result_type(.5)*(lb+ub);
	if (RT(mid) < r0) {
	  lb= mid;
	}
	else if (RT(mid) > r1) {
	  ub= mid;
	}
	else {
	  //std::cout << "Returning "<< mid << " between " << r0 << " and " << r1 << std::endl;
	  return mid;
	}
      }
    }
    CGAL_Polynomial_postcondition(0);
    return Result(0);
  }

  template <class OT>
  OT compute(const OT &r0, const OT &r1) const {
    if (std::numeric_limits<OT>::has_infinity 
	&& r1 == std::numeric_limits<OT>::infinity()) {
      return 2*r0;
    } else {
      return (r0+r1)/2.0;
    }
  }
  
  typename K::To_isolating_interval tii_;

public:

  result_type operator()(const first_argument_type &r0, const second_argument_type &r1) const
  {
    return compute<result_type>(r0, r1);
  }

};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
