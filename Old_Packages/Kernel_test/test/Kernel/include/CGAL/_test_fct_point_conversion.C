// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : _test_fct_point_conversion.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_FCT_POINT_CONVERSION_C
#define CGAL__TEST_FCT_POINT_CONVERSION_C

#ifndef CGAL__TEST_FCT_POINT_CONVERSION_H
#include <CGAL/_test_fct_point_conversion.h>
#endif // CGAL__TEST_FCT_POINT_CONVERSION_H

template <class NT>
bool
_test_fct_point_conversion(const NT& )
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
#endif // CGAL__TEST_FCT_POINT_CONVERSION_C
