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

#ifndef CGAL_POLYNOMIAL_KERNEL_SIGN_BETWEEN_ROOTS_H
#define CGAL_POLYNOMIAL_KERNEL_SIGN_BETWEEN_ROOTS_H

#include <CGAL/Polynomial/basic.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <class K>
struct Sign_between_roots
{
  typedef CGAL::Sign result_type;
  typedef typename K::Root second_argument_type;
  typedef typename K::Root third_argument_type;
  typedef typename K::Function first_argument_type;

  Sign_between_roots(){}
  Sign_between_roots(const K &k): k_(k) {
  
  }

  result_type operator()(const first_argument_type &f,
			 const second_argument_type &r0,
			 const third_argument_type &r1) const
  {
    typename K::Rational_between_roots rbr= k_.rational_between_roots_object();
    typename K::FT rat= rbr(r0,r1);
    CGAL_postcondition(second_argument_type(rat) > r0 && second_argument_type(rat) < r1);
    return k_.sign_at_object()(f, rat);
  }
protected:
  K k_;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
