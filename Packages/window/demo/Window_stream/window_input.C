// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Stefan Schirra

 
#ifndef CGAL_USE_LEDA
#include <iostream>
int main() {
  std::cout << "This demo requires LEDA Window to compile." << std::endl;
  return 0;
}
#else

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/IO/Window_stream.h>

typedef CGAL::Cartesian< double >             RepCls;
typedef CGAL::Point_2<RepCls>                 Point;
typedef CGAL::Segment_2<RepCls>               Segment;
typedef CGAL::Line_2<RepCls>                  Line;
typedef CGAL::Ray_2<RepCls>                   Ray;
typedef CGAL::Iso_rectangle_2<RepCls>         Iso;
typedef CGAL::Triangle_2<RepCls>              Triangle;
typedef CGAL::Circle_2<RepCls>                Circle;


#if defined(CGAL_USE_CGAL_WINDOW)
#define leda_red     CGAL::red
#define leda_yellow  CGAL::yellow
#endif

int
main()
{
  CGAL::Window_stream W;
  CGAL::cgalize( W);
  W.set_fg_color( leda_red);
  W.set_bg_color( leda_yellow);
  W.display();

  CGAL::set_pretty_mode( std::cout);

  W.acknowledge("Input Points");
  W.clear();
  Point p;
  while ( W >> p) { std::cout << p << std::endl; }

  W.acknowledge("Input Segments");
  W.clear();
  Segment s;
  while ( W >> s) { std::cout << s << std::endl; }

  W.acknowledge("Input Lines");
  W.clear();
  Line l;
  while ( W >> l) { std::cout << l << std::endl; }

  W.acknowledge("Input Rays");
  W.clear();
  Ray ray;
  while ( W >> ray) { std::cout << ray << std::endl; }

  W.acknowledge("Input Triangles");
  W.clear();
  Triangle t;
  while ( W >> t) { std::cout << t << std::endl; }

  W.acknowledge("Input Iso_rectangles");
  W.clear();
  Iso x;
  while ( W >> x) { std::cout << x << std::endl; }

  W.acknowledge("Input Circles");
  W.clear();
  Circle c;
  while ( W >> c) { std::cout << c << std::endl; }

  W.acknowledge("Read Functions");
  W.clear();
  W.message("Point");
  CGAL::read( W, p);
  W.message("Segment");
  CGAL::read( W, s);
  W.message("Line");
  CGAL::read( W, l);
  W.message("Ray");
  CGAL::read( W, ray);
  W.message("Triangle");
  CGAL::read( W, t);
  W.message("Iso_rectangle");
  CGAL::read( W, x);
  W.message("Circle");
  CGAL::read( W, c);


  W.acknowledge("THE END");
  return 0;
}
#endif
