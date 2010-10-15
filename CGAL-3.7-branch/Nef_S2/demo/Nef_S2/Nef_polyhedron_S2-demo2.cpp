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

#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Gmpz.h>
#include <CGAL/random_selection.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/Nef_S2/SM_items.h>
#include <CGAL/IO/Qt_widget_Nef_S2.h>

typedef CGAL::Gmpz NT;
typedef CGAL::Homogeneous<NT> Kernel;
typedef Kernel::Point_3       Point_3;
typedef Kernel::Plane_3       Plane_3;

typedef CGAL::Nef_polyhedron_S2<Kernel> Nef_polyhedron_S2;
typedef Nef_polyhedron_S2::Sphere_circle  Sphere_circle;

typedef CGAL::Creator_uniform_3<NT,Point_3>  Creator;
typedef CGAL::Random_points_in_cube_3<Point_3,Creator> Point_source;


int main(int argc, char **argv) {
  int n(5), r(0);
  if ( argc > 1 ) n = atoi( argv[1] );  // reads number of Sphere_segments
  if ( argc > 2 ) r = atoi( argv[2] );  // reads seed for random points
  srand(r);

  std::list<Sphere_circle> L;
  Point_source S(5);
  Point_3 ph;
  Point_3 o(0,0,0);
  while ( n-- > 0 ) {
    do { ph = *S++; } while ( ph == o );
    Plane_3 h(o,(ph-CGAL::ORIGIN).direction());
    L.push_back( Sphere_circle(h) );
  }

  // partition input into two lists
  Nef_polyhedron_S2 Ni, N;
  bool first(true);
  std::list<Sphere_circle>::const_iterator it;
  CGAL_forall_iterators(it,L) {
    if ( first ) {
      N = Nef_polyhedron_S2(*it);
      first = false;
    } else {
      Ni = Nef_polyhedron_S2(*it);
      N = N.symmetric_difference(Ni);
    }
  }

  QApplication a(argc, argv);
  CGAL::Qt_widget_Nef_S2<Nef_polyhedron_S2>* w =
    new CGAL::Qt_widget_Nef_S2<Nef_polyhedron_S2>(N);
  a.setMainWidget(w);
  w->show();
  return a.exec();
}
