// Copyright (c) 2011  GeometryFactory Sarl (France)
// 
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//
// Author(s)     : Laurent Rineau


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Filtered_kernel.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Cartesian<double> K2;
typedef CGAL::Simple_cartesian<CGAL::Gmpq> K3;
typedef K3::RT FT3;

typedef CGAL::Algebraic_kernel_for_circles_2_2<double>          AK;
typedef CGAL::Circular_kernel_2<K, AK> CK;

typedef CGAL::Filtered_kernel<K> FK;

int main() {
  K::Point_2 default_p;
  K::Point_2 p(-1./3, 2.);
  K::Vector_2 v = p - CGAL::ORIGIN;
  K::Circle_2 c(p, 10);

  K2::Point_2 p2(-1./3, 2.);
  K2::Vector_2 v2 = p2 - CGAL::ORIGIN;
  K2::Circle_2 c2(p2, 10);

  // no correct pretty-printer for CGAL::Gmpq
  K3::Point_2 p3(-3, 10);
  K3::Vector_2 v3 = p3 - CGAL::ORIGIN; 

  CK::Point_2 p4;
  CK::Vector_2 v4;

  FK::Point_2 default_p5;
  FK::Vector_2 default_v5;

  FK::Point_2 p5(1, 2);
  FK::Vector_2 v5(3, 4);

  return 0;
}
