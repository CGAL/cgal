// Copyright (c) 2001, 2004  Max-Planck-Institute Saarbruecken (Germany).
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
//
// Author(s)     : Michael Seel       <seel@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Gmpz.h>
#include <CGAL/random_selection.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/IO/Nef_polyhedron_S2_OGLUT_stream.h>

typedef CGAL::Gmpz NT;
typedef CGAL::Homogeneous<NT> Kernel;
typedef Kernel::Point_3       Point_3;
typedef Kernel::Plane_3       Plane_3;

typedef CGAL::Nef_polyhedron_S2<Kernel> Nef_polyhedron_S2;
typedef Nef_polyhedron_S2::Sphere_point   Sphere_point;
typedef Nef_polyhedron_S2::Sphere_segment Sphere_segment;
typedef Nef_polyhedron_S2::Sphere_circle  Sphere_circle;
typedef Nef_polyhedron_S2::Explorer Explorer;
typedef Nef_polyhedron_S2::Sphere_kernel Sphere_kernel;
typedef Nef_polyhedron_S2::Sphere_map Sphere_map;

typedef CGAL::SM_visualizor<Explorer> SM_visualizor;

typedef CGAL::Creator_uniform_3<NT,Point_3>  Creator;
typedef CGAL::Random_points_in_cube_3<Point_3,Creator> Point_source;

int main(int argc, char **argv)
{
  CGAL::set_pretty_mode ( std::cerr );
  int n(1), r(0);
  if ( argc > 1 ) n = atoi( argv[1] );
  if ( argc > 2 ) r = atoi( argv[2] );
  srand(r);

  Point_source S(5);
  Point_3 ph;
  Point_3 o(0,0,0);
  std::list<Sphere_circle> L;
  while ( n-- > 0 ) {
    do { ph = *S++; } while ( ph == o );
    Plane_3 h(o,(ph-CGAL::ORIGIN).direction());
    L.push_back( Sphere_circle(h) );
  }

  Nef_polyhedron_S2 N(L.begin(),L.end(),0.5);

  CGAL::OGL::add_sphere();
  SM_visualizor V1(&N,CGAL::OGL::spheres_.back());
  V1.draw_map();

  CGAL::OGL::add_sphere();
  SM_visualizor V2(&N,CGAL::OGL::spheres_.back());
  V2.draw_triangulation();

  CGAL::OGL::start_viewer();

  return 0;
}
