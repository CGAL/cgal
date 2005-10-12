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

#ifndef CGAL_COMPARE_ISOLATED_ROOTS_IN_INTERVAL_H
#define CGAL_COMPARE_ISOLATED_ROOTS_IN_INTERVAL_H

#include <CGAL/Polynomial/basic.h>


CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template<class Kernel, class M = CGAL::Ring_tag>
class Compare_isolated_roots_in_interval
{
protected:
  typedef M                        Method_tag;
  typedef typename Kernel::NT  NT;

  typedef typename 
  Kernel::Sign_Sturm_sequence Sign_Sturm_sequence;

public:
  Compare_isolated_roots_in_interval(const typename Kernel::Function& p,
				     const typename Kernel::Function& q,
				     Kernel k): kernel_(k),q(q), seq(kernel_.sign_Sturm_sequence_object(p, q)) {}

  Compare_isolated_roots_in_interval(){}


  typedef CGAL_POLYNOMIAL_NS::Comparison_result  result_type;
  typedef NT                       argument_type;
  typedef NT                       first_argument_type;
  typedef NT                       second_argument_type;
  
  result_type operator()(const NT& a, const NT& b) const
  {
    int sgn = seq.sum_of_signs(a, b);

    if ( sgn == CGAL_POLYNOMIAL_NS::ZERO ) { return CGAL_POLYNOMIAL_NS::EQUAL; }

    int s_a = static_cast<int>(  kernel_.sign_at_object(q)(a) );

    if ( sgn == s_a ) { return CGAL_POLYNOMIAL_NS::SMALLER; }
    return CGAL_POLYNOMIAL_NS::LARGER;
  }


protected:
  Kernel kernel_;
  typename Kernel::Function          q;
  Sign_Sturm_sequence seq;
};






CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE


#endif // CGAL_COMPARE_ISOLATED_ROOTS_IN_INTERVAL_H
