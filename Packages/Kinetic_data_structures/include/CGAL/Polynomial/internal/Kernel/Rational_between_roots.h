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

#ifndef CGAL_POLYNOMIAL_KERNEL_RATIONAL_BETWEEN_ROOTS_H
#define CGAL_POLYNOMIAL_KERNEL_RATIONAL_BETWEEN_ROOTS_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/interval_arithmetic.h>
#include <CGAL/Polynomial/internal/nt_converters.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template <class K>
struct Rational_between_roots {
  typedef typename K::NT result_type;
  typedef typename K::Root first_argument_type;
  typedef typename K::Root second_argument_type;
  
  result_type operator()(const first_argument_type &r0, const first_argument_type &r1) const {
    CGAL_precondition(r0<r1);
    typedef std::pair<double, double> Ival;
    Ival i0= to_interval(r0);
    Ival i1= to_interval(r1);
    if (r0== -infinity<first_argument_type>()) {
      if (r1== infinity<second_argument_type>()){
	return result_type(0);
      } else {
	return result_type(i1.first-1);
      }
    } else if (r1== infinity<second_argument_type>()){
      return i0.second+1;
    }
   
    // the intervals are disjoint
    if (i0.second < i1.first) {
      return result_type(.5*(i0.second+ i1.first));
    } else {
      // try to shrink the itnervals
      double d0= to_double(r0);
      double d1= to_double(r1);
      if (0) {
	++d0; ++d1;
      }
      Ival i0= to_interval(r0);
      Ival i1= to_interval(r1);

      if (i0.second < i1.first) {
	// they are now disjoint
	return result_type(.5*(i0.second+ i1.first));
      } else {
	result_type lb = i0.first;
	result_type ub = i1.second;
	result_type mid;
	while (true){
	  mid = result_type(.5)*(lb+ub);
	  if (first_argument_type(mid) < r0) {
	    lb= mid;
	  } else if (second_argument_type(mid) > r1){
	    ub= mid;
	  } else {
	    return mid;
	  }
	}
      }
    }
    CGAL_postcondition(0);
    return result_type(0);
  }
  
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

#endif
