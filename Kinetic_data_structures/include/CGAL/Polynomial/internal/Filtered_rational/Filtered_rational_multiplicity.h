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

#ifndef CGAL_POLYNOMIAL_INTERNAL_FILTERED_ROOT_MULTIPLICITY_H
#define CGAL_POLYNOMIAL_INTERNAL_FILTERED_ROOT_MULTIPLICITY_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/interval_arithmetic.h>
#include <vector>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <class Kernel>
class Filtered_rational_multiplicity
{
  typedef typename Kernel::Function::NT NT;
  typedef typename Kernel::Function Function;
public:

  Filtered_rational_multiplicity(){}

  Filtered_rational_multiplicity(const Function &fh, const Kernel& k= Kernel()):kernel_(k) {
    h_.push_back(fh);
    h_.push_back(k.differentiate_object()(fh));
  }

  typedef int result_type;
  //typedef Bound_type argument_type;
  typedef NT argument_type;

  template <class NTT>
  result_type operator()(const NTT &t) const
  {
    //typename Kernel::Sign_at sa=;
    CGAL_Polynomial_exactness_assertion( kernel_.sign_at_object(h_[0])(t)==CGAL_POLYNOMIAL_NS::ZERO);
    //if (kernel_.sign_at_object(fh)(t)!= ::CGAL::ZERO) return 0;
    // need to check if it is an even root
    int deg=1;

    typename Kernel::Interval_traits::NT i;
    {
      CGAL_POLYNOMIAL_NS::Interval_arithmetic_guard gd;
      if (h_[0].interval_function().is_zero()) return -1;
      i = kernel_.exact_to_interval_converter_object().nt_converter()(t);
    }
    //typename POLYNOMIAL_NS::To_interval<NTT> ti;

    typename Kernel::Exact_traits::NT e(t);
    do {
      if (h_.size() == static_cast<unsigned int>(deg)) {
	h_.push_back(kernel_.differentiate_object()(h_.back()));
      }
      typename Kernel::Interval_traits::NT vali;
      {
	CGAL_POLYNOMIAL_NS::Interval_arithmetic_guard gd;
	//gd.set_enabled(true); // turn on filtering
	vali = h_[deg].interval_function()(i);
	//gd.set_enabled(false); // turn off filtering
      }

      if (vali.sup() < 0 || vali.inf() >0) {
	return deg;
      }
      else if (vali.sup()==0 && vali.inf()==0) {

      }
      else {
	typename Kernel::Exact_traits::NT vale= h_[deg].exact_function()(e);
	// catch up exact and evaluate

	if (CGAL::sign(vale) == CGAL::ZERO) {

	}
	else {
	  // filtering is off
	  return deg;
	}
      }
      ++deg;
    } while (true);

    /*CGAL_postcondition(0);
      return 1;*/
  }

protected:
  Kernel kernel_;
  mutable std::vector<Function> h_;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
