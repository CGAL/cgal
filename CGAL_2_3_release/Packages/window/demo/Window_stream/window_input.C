// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : demo/Window_stream/window_input.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

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

