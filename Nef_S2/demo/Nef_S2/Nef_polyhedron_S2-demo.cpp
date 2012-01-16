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
// Author(s)     : Michael Seel       <seel@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#include <CGAL/basic.h>

#ifdef CGAL_USE_LEDA

#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/leda_integer.h>
#include <CGAL/random_selection.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/IO/Nef_polyhedron_S2_OGLUT_stream.h>
#if CGAL_LEDA_VERSION < 500
#include <LEDA/param_handler.h>
#else
#include <LEDA/system/param_handler.h>
#endif
#include <CGAL/Nef_S2/SM_items.h>

typedef leda_integer NT;
typedef CGAL::Homogeneous<NT> Kernel;
typedef Kernel::Point_3       Point_3;
typedef Kernel::Plane_3       Plane_3;

typedef CGAL::Nef_polyhedron_S2<Kernel>           Nef_polyhedron_S2;
typedef Nef_polyhedron_S2::Sphere_point           Sphere_point;
typedef Nef_polyhedron_S2::Sphere_segment         Sphere_segment;
typedef Nef_polyhedron_S2::Sphere_circle          Sphere_circle;

typedef CGAL::Creator_uniform_3<NT,Point_3>  Creator;
typedef CGAL::Random_points_in_cube_3<Point_3,Creator> Point_source;


int main(int argc, char **argv)
{
  CGAL::set_pretty_mode ( std::cerr );
  CGAL_NEF_SETDTHREAD(911);
  // Sphere_geometry 11
  // Sphere_geometry_OGL 13
  // Segment_overlay 23
  // SM_overlayer 53
  Point_3 p(0,0,0);

  int seed;
  int n;
  leda_string input_file;
  leda_param_handler H(argc,argv,".nd",false);
  H.add_parameter("random_seed:-r:int:0");
  H.add_parameter("number_of_lines:-n:int:10");
  H.add_parameter("file_of_segments:-i:string:");
  leda_param_handler::init_all();
  H.get_parameter("-r",seed);
  H.get_parameter("-n",n);
  H.get_parameter("-i",input_file);
  CGAL_assertion_msg(n>0,"-n value must be postive.");
  srand(seed);

  std::list<Sphere_circle> L;
  if ( input_file == "" ) { // create random input:
    Point_source S(5);
    Point_3 ph;
    Point_3 o(0,0,0);
    while ( n-- > 0 ) {
      do { ph = *S++; } while ( ph == o );
      Plane_3 h(o,(ph-CGAL::ORIGIN).direction());
      L.push_back( Sphere_circle(h) );
    }
  } else { // read input from file:
    std::ifstream input(input_file);
    CGAL_assertion_msg(input,"no input log.");
    Sphere_circle c;
    while ( input >> c ) L.push_back(c);
  }

  // output log:
  std::ofstream output("nef-demo.log");
  std::list<Sphere_circle>::iterator it;
  CGAL_forall_iterators(it,L) output << *it << ' ';
  output << std::endl;
  output.close();

  // partition input into two lists
  std::list<Sphere_circle> L1,L2;
  int b=0;
  CGAL_forall_iterators(it,L) {
    if ( b == 0 ) L1.push_back(*it);
    else          L2.push_back(*it);
    b = 1-b;
  }

  Nef_polyhedron_S2 N1(L1.begin(), L1.end(), 0.5);
  Nef_polyhedron_S2 N2(L2.begin(), L2.end(), 0.5);
  Nef_polyhedron_S2 N3 = N1 * N2;
  Nef_polyhedron_S2 N4 = N1 ^ N2;
  //std::cerr << N1 << N2 << N3 << N4 << std::endl;
  CGAL::ogl << N1 << N2 << N3 << N4;
  CGAL::ogl << "Nef Polyhedron 1" << "Nef Polyhedron 2"
            << "Intersection" << "Symmetric Difference";
  CGAL::ogl.display();
  return 0;

}


#else

#include <iostream>

int main()
{
  std::cout << "This demo requires LEDA\n";
  return 0;
}

#endif
