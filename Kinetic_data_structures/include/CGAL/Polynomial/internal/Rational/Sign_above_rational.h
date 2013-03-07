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

#ifndef CGAL_POLYNOMIAL_INTERNAL_SIGN_ABOVE_RATIONAL_H
#define CGAL_POLYNOMIAL_INTERNAL_SIGN_ABOVE_RATIONAL_H

#include <CGAL/Polynomial/basic.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {
template <class Kernel>
CGAL::Sign sign_above(const typename Kernel::Function &,
		      const typename Kernel::NT &, const Kernel &)
{
  // to make sure this is not called in vain
  CGAL_precondition( false );
  // to avoid warning
  return CGAL::ZERO;
}


template <class Kernel>
class Sign_above_rational
{
public:
  Sign_above_rational( Kernel k= Kernel()):k_(k){}
  typedef typename Kernel::FT second_argument_type;
  typedef typename Kernel::Function first_argument_type;
  //  typedef typename POLYNOMIAL_NS::Sign result_type;
  // g++ 3.4 does not like the above declaration
  typedef CGAL::Sign result_type;
  result_type operator()(const first_argument_type &p, const second_argument_type &nt) const
  {
    //CGAL_exactness_precondition(k.sign_at_object(p)(nt)==CGAL::ZERO);
    CGAL::Sign sn= k_.sign_at_object()(p, nt);
    if (sn != CGAL::ZERO) return sn;
    typename Kernel::Differentiate d= k_.differentiate_object();
    typename Kernel::Function pcur= d(p);
    do {
      CGAL::Sign sn= k_.sign_at_object()(pcur, nt);
      if (sn != CGAL::ZERO) return sn;
      else {
	pcur=d(pcur);
      }
    } while (pcur.degree()>=0);
    return CGAL::ZERO;
  }
protected:
  Kernel k_;
};
} } } //namespace CGAL::POLYNOMIAL::internal
#endif
