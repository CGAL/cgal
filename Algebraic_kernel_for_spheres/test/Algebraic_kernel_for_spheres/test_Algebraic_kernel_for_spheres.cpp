// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
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
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)
//
// $URL$
// $Id$
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion
//             Pedro Machado

#include <CGAL/Cartesian.h>
#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/_test_predicates.h>
#include <CGAL/_test_constructor.h>

int main()
{
  typedef CGAL::Quotient<CGAL::MP_Float>                       NT1;
  typedef CGAL::Algebraic_kernel_for_spheres_2_3<NT1>          Algebraic_k1;
  Algebraic_k1 ak1;
  _test_solve(ak1);
  _test_sign_at(ak1);
  _test_critical_points(ak1);
  _test_compare_Root_for_spheres(ak1);
  _test_constuctor(ak1);

  return 0;
}
