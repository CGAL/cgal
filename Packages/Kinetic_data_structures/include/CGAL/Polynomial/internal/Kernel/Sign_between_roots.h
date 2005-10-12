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

#ifndef CGAL_POLYNOMIAL_KERNEL_SIGN_BETWEEN_ROOTS_H
#define CGAL_POLYNOMIAL_KERNEL_SIGN_BETWEEN_ROOTS_H

#include <CGAL/Polynomial/basic.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template <class K>
struct Sign_between_roots {
  typedef CGAL_POLYNOMIAL_NS::Sign result_type;
  typedef typename K::Function argument_type;
    
  Sign_between_roots(){}
  Sign_between_roots(const typename K::Root &r0,
		     const typename K::Root &r1,
		     const K &k): k_(k){
    typename K::Rational_between_roots rbr= k.rational_between_roots_object();
    rat_= rbr(r0,r1);
  }

  result_type operator()(const argument_type &f) const {
    typename K::Sign_at sa= k_.sign_at_object(f);
    return sa(rat_);
  }
protected:
  typename K::NT rat_;
  K k_;
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

#endif
