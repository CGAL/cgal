
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
// source        : webCH2D/ch_demo.fw
// file          : demo/ConvexHull/ch_demo.C
// revision      : 3.3
// revision_date : 03 Aug 2000 
// author(s)     : Stefan Schirra
//
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de> 
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/Cartesian.h>

#ifndef CGAL_USE_LEDA
#error Sorry, this demo needs LEDA for visualisation
#endif // CGAL_USE_LEDA

#include <CGAL/leda_real.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <vector>
#include <CGAL/point_generators_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/copy_n.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/segment_intersection_points_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2.h>
#include <list>
#include <CGAL/IO/Window_stream.h>

typedef CGAL::Cartesian<double>                               R;
typedef CGAL::Cartesian<leda_real>                            Real_R;
typedef CGAL::Point_2<R>                                      Point;
typedef CGAL::Point_2<Real_R>                                 Real_Point;
typedef CGAL::Creator_uniform_2<double,Point>                 Point_creator;
typedef CGAL::Segment_2<R>                                    Segment;
typedef CGAL::Segment_2<Real_R>                               Real_Segment;
typedef CGAL::Random_points_on_circle_2<Point,Point_creator>  Source;
typedef CGAL::Creator_uniform_2< Point, Segment>              Segment_creator;
typedef CGAL::Join_input_iterator_2< Source, Source, Segment_creator>
                                                              Segment_iterator;
typedef CGAL::Polygon_traits_2<R>                             PolygonTraits;
typedef std::list<Point>                                      Container;
typedef CGAL::Polygon_2<PolygonTraits,Container>              Polygon;

int
main()
{
  CGAL::Window_stream W("Convex Hull of Intersection Points");
  W.init(-256.0, 255.0, -256.0);
  CGAL::cgalize(W);
  W.display();
  Point p;

  Source S( 250);
  Segment_iterator g( S, S);
  std::vector<Segment> Vs;
  CGAL::copy_n( g, 30, std::back_inserter( Vs) );
  W << CGAL::GREEN;
  cout << "Generating random segments." << endl;
  std::copy( Vs.begin(), Vs.end(),
             CGAL::Ostream_iterator<Segment,CGAL::Window_stream>( W) );
  cout << "Click in the window to see the intersection points." << endl;
  W.read_mouse();

  std::vector<Point>   Vip;
  CGAL::segment_intersection_points_2( Vs.begin(), Vs.end(),
                                       std::back_inserter( Vip),
                                       R() );
  std::cout << "The " << Vip.size();
  std::cout << " intersection points are now shown in blue in the window";
  std::cout << std::endl;
  W << CGAL::BLUE;
  std::copy( Vip.begin(), Vip.end(),
             CGAL::Ostream_iterator<Point,CGAL::Window_stream>( W) );

  cout << "Click in the window to see the convex hull of the intersection"
       << " points." << endl;
  W.read_mouse();

  Polygon Pol;
  CGAL::convex_hull_points_2( Vip.begin(), Vip.end(),
                              std::back_inserter( Pol) );
  W << CGAL::RED;
  W.set_line_width(2);
  W << Pol;
  cout << "Click in the window to exit." << endl; 
  W.read_mouse();
}
