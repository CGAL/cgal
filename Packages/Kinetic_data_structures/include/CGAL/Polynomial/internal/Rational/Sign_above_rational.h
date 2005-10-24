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

#ifndef CGAL_POLYNOMIAL_INTERNAL_SIGN_ABOVE_RATIONAL_H
#define CGAL_POLYNOMIAL_INTERNAL_SIGN_ABOVE_RATIONAL_H

#include <CGAL/Polynomial/basic.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE
template <class Kernel>
CGAL_POLYNOMIAL_NS::Sign sign_above(const typename Kernel::Function &p,
		      const typename Kernel::NT &nt, const Kernel &k){
  // to make sure this is not called in vain
  CGAL_precondition( false );
  // to avoid warning
  return 0;
}



template <class Kernel>
class Sign_above_rational {
public:
  Sign_above_rational(){}
  Sign_above_rational(const typename Kernel::Function &p, Kernel k= Kernel()):p_(p), k_(k){}
  typedef typename Kernel::NT argument_type;
  //  typedef typename POLYNOMIAL_NS::Sign result_type;
  // g++ 3.4 does not like the above declaration
  typedef CGAL_POLYNOMIAL_NS::Sign result_type;
  result_type operator()(const argument_type &nt) const {
    //CGAL_exactness_precondition(k.sign_at_object(p)(nt)==CGAL::ZERO);
    CGAL_POLYNOMIAL_NS::Sign sn= k_.sign_at_object(p_)(nt);
    if (sn != CGAL_POLYNOMIAL_NS::ZERO) return sn;
    typename Kernel::Differentiate d= k_.differentiate_object();
    typename Kernel::Function pcur= d(p_);
    do {
      CGAL_POLYNOMIAL_NS::Sign sn= k_.sign_at_object(pcur)(nt);
      if (sn != CGAL_POLYNOMIAL_NS::ZERO) return sn;
      else {
	pcur=d(pcur);
      }
    } while (!pcur.is_zero());
    return CGAL_POLYNOMIAL_NS::ZERO;
  }
protected:
  typename Kernel::Function p_;
  Kernel k_;
};
CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
