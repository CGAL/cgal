// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://pmachado@scm.gforge.inria.fr/svn/cgal/trunk/Circular_kernel_3/test/Circular_kernel_3/test_Spherical_kernel.cpp $
// $Id: test_Spherical_kernel.cpp 44381 2008-07-24 09:26:28Z pmachado $
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifdef CGAL_USE_CORE
#include <CGAL/Cartesian.h>
#include <CGAL/Spherical_kernel_3.h>
#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Gmpq.h>
#include <CGAL/CORE_BigRat.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/_test_sphere_predicates.h>
#include <CGAL/_test_sphere_constructions.h>
#include <CGAL/_test_sphere_compute.h>
#include <CGAL/Polynomials_1_3.h>
#include <CGAL/Polynomials_2_3.h>
#include <CGAL/Polynomials_for_line_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/_test_functionalities_on_sphere.h>
#endif

int main()
{
#ifdef CGAL_USE_CORE	
  typedef CORE::Expr                                          FT;
  typedef CGAL::Cartesian<FT>                                 Linear_k1;
  typedef CGAL::Algebraic_kernel_for_spheres_2_3<FT>          Algebraic_k1;
  typedef CGAL::Spherical_kernel_3<Linear_k1,Algebraic_k1>    SK1;
  SK1 sk1;
  _test_spherical_kernel_predicates(sk1);
  _test_spherical_kernel_construct(sk1); 
  _test_spherical_kernel_compute(sk1);
  test_functionalities_on_a_reference_sphere<SK1>(SK1::Point_3(4.5,1.2,0.7));
#endif
  return 0;
}
