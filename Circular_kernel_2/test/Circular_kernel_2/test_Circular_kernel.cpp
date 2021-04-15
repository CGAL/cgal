// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/_test_circles_predicates.h>
#include <CGAL/_test_circles_constructions.h>
#include <CGAL/_test_circles_extention.h>

int main()
{
  typedef CGAL::Quotient<CGAL::MP_Float>                       NT1;
  typedef CGAL::Cartesian<NT1>                                 Linear_k1;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT1>          Algebraic_k1;
  typedef CGAL::Circular_kernel_2<Linear_k1,Algebraic_k1>      CK1;
  CK1 ck1;
  std::cout << "Testing predicates..." << std::endl;
  _test_circle_predicat(ck1);
  std::cout << "Testing constructions..." << std::endl;
  _test_circle_construct(ck1);
  std::cout << "Testing bboxes..." << std::endl;
  _test_circle_bbox(ck1);
  std::cout << "Testing circular_arc_bboxes..." << std::endl;
  _test_circular_arc_bbox(ck1);
  std::cout << "Testing circular_arc_point_bboxes..." << std::endl;
  _test_circular_arc_point_bbox(ck1);
  std::cout << "Testing has_on..." << std::endl;
  _test_has_on(ck1);
  return 0;
}
