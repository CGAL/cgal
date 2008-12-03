// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#include <CGAL/Exact_spherical_kernel_3.h>

typedef CGAL::Exact_spherical_kernel_3::FT FT_Q;

#include <CGAL/_test_sphere_predicates.h>
#include <CGAL/_test_sphere_constructions.h>
#include <CGAL/_test_sphere_compute.h>

int main()
{ 
  CGAL::Exact_spherical_kernel_3  sk1;
  _test_spherical_kernel_predicates(sk1);
  _test_spherical_kernel_construct(sk1); 
  _test_spherical_kernel_compute(sk1);
  return 0;
}

