// Copyright (c) 2002  Max Planck Institut fuer Informatik (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Susan Hert

#include <CGAL/config.h>

#ifdef CGAL_USE_GEOMVIEW

#include <CGAL/Homogeneous.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/predicates_on_points_3.h>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer RT;
#else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz RT;
#else
// NOTE: the choice of double here for a number type may cause problems
//       for degenerate point sets
#include <CGAL/double.h>
typedef double RT;
#endif
#endif

#include <vector>
#include <cstdlib>

// NOTE: the choice of double here for a number type may cause problems
//       for degenerate point sets
typedef CGAL::Homogeneous<RT>                     K;
typedef CGAL::Convex_hull_traits_3<K>             Traits;
typedef Traits::Polyhedron_3                      Polyhedron_3;

// define point creator
typedef K::Point_3                                Point_3;
typedef CGAL::Creator_uniform_3<double, Point_3>  PointCreator;
typedef CGAL::Random_points_in_sphere_3<Point_3, PointCreator> Generator;

void draw_points_and_hull(const std::vector<Point_3>& points,
                          const CGAL::Object& object)
{
   std::vector<Point_3>::const_iterator p_it;

   CGAL::Geomview_stream geomview;
   geomview << CGAL::red();
   for (p_it = points.begin(); p_it != points.end(); p_it++)
   {
      geomview << *p_it;
   }

   K::Segment_3     segment;
   K::Triangle_3    triangle;
   Point_3          point;
   Polyhedron_3     polyhedron;

   geomview << CGAL::blue();
   if ( CGAL::assign(point, object) )
      geomview << point;
   else if ( CGAL::assign(segment, object) )
      geomview << segment;
   else if ( CGAL::assign(triangle, object) )
      geomview << triangle;
   else if ( CGAL::assign(polyhedron, object) )
      geomview << polyhedron;


   std::cout << "Press any key to end the program: ";
   char wait;
   std::cin >> wait;
}



int main(int argc, char* argv[])
{

  if (argc != 2)
  {
      std::cerr << "Usage: " << argv[0] << " #points " << std::endl;
      std::exit(0);
  }

  int num = std::atoi(argv[1]);
  if (num < 0)
  {
     std::cerr << "Usage: " << argv[0] << " #points " << std::endl;
     std::cerr << " #points must be >= 0" << std::endl;
     std::exit(0);
  }

  std::vector<Point_3> points;
  Generator gen(100.0);

  // generate num points and copy them to a vector
  std::copy_n( gen, num, std::back_inserter(points) );

  // define object to hold convex hull
  CGAL::Object ch_object;

  // compute convex hull
  CGAL::convex_hull_3(points.begin(), points.end(), ch_object);
  draw_points_and_hull(points, ch_object);
  return 0;
}

#else

#include <iostream>

int main() {
  std::cerr <<
  "This demo requires geomview, which is not present on this platform\n";
  return 0;
}

#endif
