// Copyright (c) 2001, 2004  Max-Planck-Institute Saarbruecken (Germany).
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
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#include <iostream>
#include <string>
#include <cstdlib>

#if !defined(_MSC_VER)

#include <CGAL/Homogeneous.h>
#include <CGAL/Convex_hull_d_traits_3.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/random_selection.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Convex_hull_d_to_polyhedron_3.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer RT;
#else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz RT;
#else
#include <CGAL/double.h>
typedef double RT;
#endif
#endif

typedef CGAL::Homogeneous<RT> Kernel_3;
typedef CGAL::Convex_hull_d_traits_3<Kernel_3> Kernel_d_3;
typedef CGAL::Convex_hull_d<Kernel_d_3> Convex_hull_d;

typedef  CGAL::Point_3<Kernel_3> Point_3;
typedef  CGAL::Polyhedron_3<Kernel_3> Polyhedron;

typedef CGAL::Creator_uniform_3<RT,Point_3> Creator;
typedef CGAL::Random_points_in_cube_3<Point_3,Creator> Point_source;

int main(int argc, char* argv[]) {
  int dimension = 3;
  int n = 100;
  if (argc > 1 && std::string(argv[1])=="-h") {
    std::cout << "usage: ch5-demo [#points]\n";
    std::exit(1);
  }
  if (argc > 1) n = std::atoi(argv[1]);

  int r = 2*n;
  CGAL::Geomview_stream gv(CGAL::Bbox_3(-r, -r, -r, r, r, r));
  gv.clear();


  Convex_hull_d T(dimension);
  Point_source g(n);
  for(int i=0; i<n; ++i) {
    T.insert(*g++);
    if (i%10==0) std::cout << i << " points inserted" << std::endl;
  }
  T.is_valid(true);

  Polyhedron P;
  CGAL::convex_hull_d_to_polyhedron_3(T,P);
  gv << P;
  std::cout << "Enter a key to finish" << std::endl;
  char ch;
  std::cin >> ch;

  return 0;
}

#else // on windows:

int main(int argc, char* argv[]) {
  std::cerr <<
  "This demo requires geomview, that is is not present on windows\n";
  return 0;
}

#endif
