// Copyright (c) 2002-2011  Max Planck Institut fuer Informatik (Germany).
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
// $URL: svn+ssh://sloriot@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/Faster_static_convex_hull_3/demo/Convex_hull_3/incremental_hull_3_demo.cpp $
// $Id: incremental_hull_3_demo.cpp 44910 2008-08-12 12:58:18Z spion $
//
//
// Author(s)     : Susan Hert
//


#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/Convex_hull_d_traits_3.h>
#include <CGAL/Convex_hull_d_to_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/algorithm.h>
#include <vector>
#include <cassert>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpq.h>
typedef CGAL::Gmpq RT;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::Quotient<CGAL::MP_Float> RT;
#endif

typedef CGAL::Cartesian<RT>                        K;
typedef K::Point_3                                 Point_3;
typedef CGAL::Polyhedron_3< K>                     Polyhedron_3;

typedef CGAL::Convex_hull_d_traits_3<K>            Hull_traits_3;
typedef CGAL::Convex_hull_d< Hull_traits_3 >       Convex_hull_3;
typedef CGAL::Creator_uniform_3<double, Point_3>   Creator;

int main ()
{
  Convex_hull_3 CH(3);  // create instance of the class with dimension == 3

  // generate 250 points randomly on a sphere of radius 100
  // and insert them into the convex hull
  CGAL::Random_points_in_sphere_3<Point_3, Creator> gen(100);

  for (int i = 0; i < 250 ; i++, ++gen)
     CH.insert(*gen);

  assert(CH.is_valid());

  // define polyhedron to hold convex hull and create it
  Polyhedron_3 P;
  CGAL::convex_hull_d_to_polyhedron_3(CH,P);

  std::cout << "The convex hull has " << P.size_of_vertices() 
            << " vertices" << std::endl;
  return 0;
}
