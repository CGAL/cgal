
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
// file          : examples/ConvexHull/ch_of_polyline.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  
// ============================================================================
 

#include <CGAL/Cartesian.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/ch_melkman.C>

#ifdef CGAL_USE_LEDA
#include <CGAL/IO/leda_ps_file.h>
#include <CGAL/IO/leda_window.h>
#endif // CGAL_USE_LEDA

#include <CGAL/IO/polygonal_2.h>
#include <CGAL/IO/Ostream_iterator.h>

typedef  CGAL::Cartesian<double>      K;
typedef  K::Point_2                   Point_2;

CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(Point_2*)

int main()
{
  CGAL::set_ascii_mode(std::cin);
  CGAL::set_ascii_mode(std::cout);
  std::istream_iterator< Point_2>  in_start( std::cin );
  std::istream_iterator< Point_2>  in_end;
  std::vector<Point_2> V ;
  std::copy( in_start, in_end, std::back_inserter(V) );
#ifdef CGAL_USE_LEDA
  std::vector<Point_2> CH;
#ifdef NO_DISPLAY
  typedef CGAL::ps_file  OutStream;
#else
  typedef CGAL::Window_stream  OutStream;
#endif
  OutStream W;
  CGAL::Bbox_2 b = V.begin()->bbox();
  for ( std::vector<Point_2>::iterator it = V.begin(); it != V.end(); ++it)
  { b = b + (*it).bbox(); }
  double x_span = b.xmax() - b.xmin();
  double y_span = b.ymax() - b.ymin();
  double span = CGAL::max( x_span, y_span);
  span *= 1.1;
  W.init((b.xmin()+b.xmax()-span)/2,
         (b.xmin()+b.xmax()+span)/2,
         (b.ymin()+b.ymax()-span)/2 );
  CGAL::cgalize(W);
  W.display();
  CGAL::ch_melkman( V.begin(), V.end(), std::back_inserter(CH) );
  std::cout << "Convex Hull has size " << CH.size() << std::endl;
  W << CGAL::Color(200,200,200);
  CGAL::send_to_stream_as_polygon( W, CH.begin(), CH.end());
  W << CGAL::GREEN;
  CGAL::send_to_stream_as_polygon( W, V.begin(), V.end());
  W << CGAL::RED;
  std::copy( CH.begin(), CH.end(),
             CGAL::Ostream_iterator<Point_2,OutStream>(W));
  leda_wait( 5);
#else
  CGAL::ch_melkman( V.begin(), V.end(),
                    std::ostream_iterator<Point_2>(std::cout,"\n"));
#endif //   CGAL_USE_LEDA
  return 0;
}

