
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
// file          : demo/ConvexHull/ch_demo.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
// ============================================================================
 

#include <CGAL/Cartesian.h>

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

typedef CGAL::Cartesian<double>                               K;
typedef CGAL::Point_2<K>                                      Point;
typedef CGAL::Creator_uniform_2<double,Point>                 Point_creator;
typedef CGAL::Segment_2<K>                                    Segment;
typedef CGAL::Random_points_on_circle_2<Point,Point_creator>  Source;
typedef CGAL::Creator_uniform_2< Point, Segment>              Segment_creator;
typedef CGAL::Join_input_iterator_2< Source, Source, Segment_creator>
                                                              Segment_iterator;
typedef CGAL::Polygon_traits_2<K>                             PolygonTraits;
typedef std::list<Point>                                      Container;
typedef CGAL::Polygon_2<PolygonTraits,Container>              Polygon;

int main()
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
  std::cout << "Generating random segments." << std::endl;
  std::copy( Vs.begin(), Vs.end(),
             CGAL::Ostream_iterator<Segment,CGAL::Window_stream>( W) );
  std::cout << "Click in the window to see the intersection points." 
            << std::endl;
  W.read_mouse();

  std::vector<Point>   Vip;
  CGAL::segment_intersection_points_2( Vs.begin(), Vs.end(),
                                       std::back_inserter( Vip),
                                       K() );
  std::cout << "The " << Vip.size();
  std::cout << " intersection points are now shown in blue in the window";
  std::cout << std::endl;
  W << CGAL::BLUE;
  std::copy( Vip.begin(), Vip.end(),
             CGAL::Ostream_iterator<Point,CGAL::Window_stream>( W) );

  std::cout << "Click in the window to see the convex hull of the intersection"
            << " points." << std::endl;
  W.read_mouse();

  Polygon Pol;
  CGAL::convex_hull_points_2( Vip.begin(), Vip.end(),
                              std::back_inserter( Pol) );
  W << CGAL::RED;
  W.set_line_width(2);
  W << Pol;
  std::cout << "Click in the window to exit." << std::endl; 
  W.read_mouse();

  return 0;
}
