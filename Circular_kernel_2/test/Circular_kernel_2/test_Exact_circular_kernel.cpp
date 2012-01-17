// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#define CGAL_NO_DEPRECATION_WARNINGS
#include <CGAL/basic.h>
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Circular_kernel_intersections.h>
#include <iostream>

typedef CGAL::Exact_circular_kernel_2 CK;

  CK ck;

#include <CGAL/_test_circles_predicates.h>
#include <CGAL/_test_circles_constructions.h>
#include <CGAL/_test_circles_extention.h>
  
int main() {

  _test_circle_predicat(ck);
  _test_circle_construct(ck);
  _test_circle_bbox(ck);
  _test_circular_arc_bbox(ck);
  _test_has_on(ck);

  return 0;
}
