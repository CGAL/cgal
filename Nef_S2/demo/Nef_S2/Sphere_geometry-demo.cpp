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

#include <CGAL/leda_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/algorithm.h>
#include <CGAL/random_selection.h>
#include <CGAL/point_generators_3.h>
#if CGAL_LEDA_VERSION < 500
#include <LEDA/param_handler.h>
#include <LEDA/random.h>
#else
#include <LEDA/system/param_handler.h>
#include <LEDA/core/random.h>
#endif
// #include <CGAL/Nef_S2/leda_sphere_map.h>
#include <CGAL/Nef_S2/Sphere_geometry_OGL.h>

typedef leda_integer                  RT;
typedef CGAL::Homogeneous<RT>         HKernel;
typedef CGAL::Sphere_point<HKernel>   SPoint;
typedef CGAL::Sphere_segment<HKernel> SSegment;
typedef CGAL::Plane_3<HKernel>        Plane;
typedef CGAL::Point_3<HKernel>        Point;
typedef CGAL::Direction_3<HKernel>    Direction;

typedef CGAL::Creator_uniform_3<leda_integer,Point>  Creator;
typedef CGAL::Random_points_in_cube_3<Point,Creator> Point_source;

int main(int argc, char **argv)
{
  CGAL::set_pretty_mode ( std::cerr );
  CGAL_NEF_SETDTHREAD(101); //(11*23*31);
  // Sphere_geometry 11
  // Sphere_geometry_OGL 13
  // PM_segment_overlay 23
  // leda_sphere_map 31

  int n;
  leda_string input_file;
  leda_param_handler H(argc,argv,".sg",false);
  H.add_parameter("number_of_lines:-n:int:10");
  H.add_parameter("file_of_segments:-i:string:");
  leda_param_handler::init_all();
  H.get_parameter("-n",n);
  H.get_parameter("-i",input_file);
  CGAL_assertion_msg(n>0,"-n value must be postive.");

  std::list<SSegment> L;
  if ( input_file == "" ) {
    Point_source S(5);
    Point p1,p2,ph;
    Point o(0,0,0);
    while ( n-->0 ) {
      do { ph = *S++; } while ( ph == o );
      Plane h(o,(ph-CGAL::ORIGIN).direction());
      do { p1 = *S++; }
      while ( p1 == o || h.projection(p1) == o );
      do { p2 = *S++; }
      while ( p2 == o || h.projection(p2) == o );
      SPoint p3(h.projection(p1)),p4(h.projection(p2));
      int which = CGAL::default_random.get_int(0,3);
      if ( p3 == p4 ) which = 3;
      if ( p3 == CGAL::ORIGIN - p4 ) which = 2;
      switch ( which ) {
        case 0: // short
          L.push_back( SSegment(p3,p4,true) ); break;
        case 1: // long
          L.push_back( SSegment(p3,p4,false) ); break;
        case 2: // halfcircle
          L.push_back( SSegment(p3,CGAL::ORIGIN - p3,h) ); break;
        case 3: // trivial
          L.push_back( SSegment(p3,p3,h) ); break;
      }
    }
  } else {
    std::ifstream input(input_file);
    CGAL_assertion_msg(input,"no input log.");
    SSegment s;
    while ( input >> s ) L.push_back(s);
  }
  std::ofstream output("sg-demo.log");
  std::list<SSegment>::iterator it;
  CGAL_forall_iterators(it,L) output << *it;
  output << std::endl;
  output.close();
  leda_sphere_map_overlayer<HKernel> SMO;
  SMO.subdivide(L.begin(),L.end());
  //SMO.dump(std::cerr);

  CGAL::OGL::add_sphere();
  CGAL::OGL::add_sphere();
  CGAL::OGL::Unit_sphere& S1(*CGAL::OGL::spheres_.begin());
  CGAL::OGL::Unit_sphere& S2(CGAL::OGL::spheres_.back());

  CGAL_forall_iterators(it,L) {
    S1.push_back(*it);
    S1.push_back(it->source(),CGAL::BLUE);
    S1.push_back(it->target(),CGAL::BLUE);
  }

  leda_edge e;
  forall_edges(e,SMO.sphere_map()) {
    SSegment s(SMO.sphere_map()[source(e)],
               SMO.sphere_map()[target(e)]);
    S2.push_back(s);
  }
  leda_node v;
  forall_nodes(v,SMO.sphere_map()) {
    SPoint p(SMO.sphere_map()[v]);
    S2.push_back(p);
  }
  CGAL::OGL::titles_.push_back(std::string("Input segments"));
  CGAL::OGL::titles_.push_back(std::string("Output map"));
  CGAL::OGL::start_viewer();
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
