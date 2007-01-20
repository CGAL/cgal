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

#ifndef CGAL_POLYNOMIAL_KERNEL_SIGN_BETWEEN_ROOTS_H
#define CGAL_POLYNOMIAL_KERNEL_SIGN_BETWEEN_ROOTS_H

#include <CGAL/Polynomial/basic.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template <class K>
struct Sign_between_roots
{
  typedef CGAL::Sign result_type;
  typedef typename K::Function argument_type;

  Sign_between_roots(){}
  Sign_between_roots(const typename K::Root &r0,
		     const typename K::Root &r1,
		     const K &k): k_(k) {
    typename K::Rational_between_roots rbr= k.rational_between_roots_object();
    rat_= rbr(r0,r1);
  }

  result_type operator()(const argument_type &f) const
  {
    typename K::Sign_at sa= k_.sign_at_object(f);
    return sa(rat_);
  }
protected:
  typename K::FT rat_;
  K k_;
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
