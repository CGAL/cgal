// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra
//                 Ron Wein        <wein@post.tau.ac.il>
 

#ifndef CGAL__TEST_FCT_POINT_CONVERSION_H
#define CGAL__TEST_FCT_POINT_CONVERSION_H

#include <CGAL/Homogeneous.h>
#include <cassert>
#include <CGAL/Cartesian.h>
#include <CGAL/cartesian_homogeneous_conversion.h>

// Test for number-types that support division.
template <class NT>
bool
_test_fct_point_conversion(const NT&, CGAL::Field_tag)
{
  std::cout << "Testing Point Conversion Functions" ;

  typedef CGAL::Homogeneous<NT>               H;
  typedef CGAL::Cartesian<NT>                 C;
  typedef CGAL::Cartesian<CGAL::Quotient<NT> > QC;

  typedef CGAL::Point_2< H >                  PtH2;
  typedef CGAL::Point_2< C >                  PtC2;
  typedef CGAL::Point_2< QC>                  PtQC2;
  typedef CGAL::Point_3< H >                  PtH3;
  typedef CGAL::Point_3< C >                  PtC3;
  typedef CGAL::Point_3< QC>                  PtQC3;

  // NT n0  = NT(0);
  // NT n1  = NT(1);
  NT n3  = NT(3);
  NT n5  = NT(5);
  NT n10 = NT(10);
  NT n15 = NT(15);
  NT n25 = NT(25);
  NT n50 = NT(50);
  CGAL::Quotient<NT> q3(n3);
  CGAL::Quotient<NT> q5(n5);
  CGAL::Quotient<NT> q10(n10);

  PtH2  ph2(n25, n15, n5);
  PtC2  pc2(n25, n15, n5);
  PtQC2 pq2(q5, q3);
  PtH3  ph3(n50, n25, n15, n5);
  PtC3  pc3(n50, n25, n15, n5);
  PtQC3 pq3(q10, q5, q3);

  std::cout << '.';

  assert( CGAL::homogeneous_to_cartesian(ph2) == pc2 );
  assert( CGAL::cartesian_to_homogeneous(pc2) == ph2 );
  assert( CGAL::homogeneous_to_quotient_cartesian(ph2) == pq2);
  assert( CGAL::quotient_cartesian_to_homogeneous(pq2) == ph2);
  std::cout << '.';

  assert( CGAL::homogeneous_to_cartesian(ph3) == pc3 );
  assert( CGAL::cartesian_to_homogeneous(pc3) == ph3 );
  assert( CGAL::homogeneous_to_quotient_cartesian(ph3) == pq3);
  assert( CGAL::quotient_cartesian_to_homogeneous(pq3) == ph3);

  std::cout << "done" << std::endl;
  return true;
}

// Test for number-types that do not support division.
template <class NT>
bool
_test_fct_point_conversion(const NT&, CGAL::Integral_domain_without_division_tag)
{
  std::cout << "Testing Point Conversion Functions" ;

  typedef CGAL::Homogeneous<NT>               H;
  typedef CGAL::Cartesian<CGAL::Quotient<NT> > QC;

  typedef CGAL::Point_2< H >                  PtH2;
  typedef CGAL::Point_2< QC>                  PtQC2;
  typedef CGAL::Point_3< H >                  PtH3;
  typedef CGAL::Point_3< QC>                  PtQC3;

  NT n3  = NT(3);
  NT n5  = NT(5);
  NT n10 = NT(10);
  NT n15 = NT(15);
  NT n25 = NT(25);
  NT n50 = NT(50);
  CGAL::Quotient<NT> q3(n3);
  CGAL::Quotient<NT> q5(n5);
  CGAL::Quotient<NT> q10(n10);

  PtH2  ph2(n25, n15, n5);
  PtQC2 pq2(q5, q3);
  PtH3  ph3(n50, n25, n15, n5);
  PtQC3 pq3(q10, q5, q3);

  std::cout << '.';

  assert( CGAL::homogeneous_to_quotient_cartesian(ph2) == pq2);
  assert( CGAL::quotient_cartesian_to_homogeneous(pq2) == ph2);
  std::cout << '.';

  assert( CGAL::homogeneous_to_quotient_cartesian(ph3) == pq3);
  assert( CGAL::quotient_cartesian_to_homogeneous(pq3) == ph3);

  std::cout << "done" << std::endl;
  return true;
}

template <class NT>
bool
_test_fct_point_conversion (const NT& x)
{
    typedef CGAL::Algebraic_structure_traits<NT> AST;
    return _test_fct_point_conversion (x, typename AST::Algebraic_category());
}

#endif // CGAL__TEST_FCT_POINT_CONVERSION_H
