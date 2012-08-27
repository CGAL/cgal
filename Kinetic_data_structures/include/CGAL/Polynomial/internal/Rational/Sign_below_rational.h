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

#ifndef CGAL_POLYNOMIAL_INTERNAL_SIGN_BELOW_RATIONAL_H
#define CGAL_POLYNOMIAL_INTERNAL_SIGN_BELOW_RATIONAL_H

#include <CGAL/Polynomial/basic.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <class Kernel>
CGAL::Sign sign_below(const typename Kernel::Function &p,
const typename Kernel::NT &nt, const Kernel &k)
{
//CGAL_exactness_precondition(k.sign_at_object(p)(nt)==CGAL::ZERO);
    CGAL::Sign sn= k.sign_at_object(p)(nt);
    if (sn != CGAL::ZERO) return sn;
    typename Kernel::Function pcur=k.differentiate_object()(p);
    int count=1;
    do {
        CGAL::Sign sn= k.sign_at_object(pcur)(nt);
        if (sn != CGAL::ZERO) {
            if ((count %2==1 && sn== CGAL::POSITIVE)
            || (count%2==0 && sn==CGAL::NEGATIVE)) {
                return CGAL::NEGATIVE;
            }
            else {
                return CGAL::POSITIVE;
            }
        }
	pcur= k.differentiate_object()(pcur);
	if (pcur.degree() <0) return  CGAL::ZERO;
    } while (1);
    //return CGAL_POLYNOMIAL_NS::ZERO;
}


template <class Kernel>
class Sign_below_rational
{
public:
  Sign_below_rational(Kernel k= Kernel()):k_(k){}
  typedef typename Kernel::Function first_argument_type;
  typedef typename Kernel::NT second_argument_type;
  //  typedef typename POLYNOMIAL_NS::Sign result_type;
  // g++ 3.4 does not like the above typedef
  typedef CGAL::Sign result_type;
  result_type operator()(const first_argument_type &f, const second_argument_type &n) const {
    return sign_below(f,n, k_);
  }
protected:
  Kernel k_;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
