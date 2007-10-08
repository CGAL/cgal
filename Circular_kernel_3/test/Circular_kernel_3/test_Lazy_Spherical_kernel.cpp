// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
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
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)
//
// $URL: svn+ssh://sloriot@scm.gforge.inria.fr/svn/cgal/trunk/Circular_kernel_3/test/Circular_kernel_3/test_Spherical_kernel.cpp $
// $Id: test_Spherical_kernel.cpp 36659 2007-02-28 11:37:19Z afabri $
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//             Pedro Machado    <tashimir@gmail.com>


#include <CGAL/Cartesian.h>
#include <CGAL/Spherical_kernel_3.h>
#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Gmpq.h>
typedef CGAL::Gmpq FT_Q;
#include <CGAL/_test_sphere_predicates.h>
#include <CGAL/_test_sphere_constructions.h>
#include <CGAL/_test_sphere_compute.h>
#include <CGAL/Polynomials_1_3.h>
#include <CGAL/Polynomials_2_3.h>
#include <CGAL/Polynomials_for_line_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

int pipo(int r){return r;}

int main()
{
  typedef CGAL::Exact_predicates_exact_constructions_kernel Linear_k1;
  typedef Linear_k1::FT FT;
  typedef CGAL::Algebraic_kernel_for_spheres_2_3<FT>          Algebraic_k1;
  typedef CGAL::Spherical_kernel_3<Linear_k1,Algebraic_k1>    SK1;
  SK1 sk1;
  _test_spherical_kernel_predicates(sk1);
  _test_spherical_kernel_construct(sk1); 
  _test_spherical_kernel_compute(sk1);
  return 0;
}
