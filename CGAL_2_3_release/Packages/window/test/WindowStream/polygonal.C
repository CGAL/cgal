
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
// file          : include/CGAL/test/WindowStream/polygonal.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/Cartesian.h>
#include <fstream>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#ifdef CGAL_USE_LEDA
#  include <CGAL/IO/leda_window.h>
#endif // CGAL_USE_LEDA
#include <CGAL/IO/polygonal_2.h>

typedef   CGAL::Point_2<CGAL::Cartesian<double> >        Point2;

int
main(int argc, char *argv[])
{
  if (argc >= 2) { return 0; }

  std::ifstream F("polygon_coordinates");
  CGAL::set_ascii_mode(F);
  CGAL::set_ascii_mode(std::cout);
  std::istream_iterator< Point2>  in_start(F);
  std::istream_iterator< Point2>  in_end;
  std::vector<Point2> V ;
  std::copy( in_start, in_end, std::back_inserter(V) );
  #ifdef CGAL_USE_LEDA
  CGAL::Window_stream W;
  CGAL::Bbox_2 b = V.begin()->bbox();
  for ( std::vector<Point2>::iterator it = V.begin(); it != V.end(); ++it)
  { b = b + (*it).bbox(); }
  double x_span = b.xmax() - b.xmin();
  double y_span = b.ymax() - b.ymin();
  double span = std::max( x_span, y_span);
  span *= 1.1;
  W.init((b.xmin()+b.xmax()-span)/2,
         (b.xmin()+b.xmax()+span)/2,
         (b.ymin()+b.ymax()-span)/2 );
  CGAL::cgalize(W);
  W.display();
  W << CGAL::RED;
  CGAL::send_to_stream_as_polygon( W, V.begin(), V.end());
  leda_wait( 5);
  W.clear();
  W << CGAL::GREEN;
  CGAL::send_to_stream_as_polyline( W, V.begin(), V.end());
  leda_wait( 5);
  #else
  std::cout << std::endl << "Polygon" << std::endl;
  CGAL::send_to_stream_as_polygon( std::cout, V.begin(), V.end());
  std::cout << std::endl << "Polyline" << std::endl;
  CGAL::send_to_stream_as_polyline( std::cout, V.begin(), V.end());
  #endif // LEDA

  return 0;
}
