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
// source        : intersection_and_then.fw
// file          : demo/Robustness/hull_of_intersection_points_2.C
// revision      : 1.5
// revision_date : 20 Sep 2000 
// author(s)     : Stefan Schirra
//
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de> 
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/basic.h>
#ifndef CGAL_USE_LEDA
int main() { std::cout << "\nSorry, this demo needs LEDA\n"; return 0; }
#else
#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#include <CGAL/leda_real.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Circle_2.h>
#include <vector>
#include <fstream>
#include <CGAL/segment_intersection_points_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/copy_n.h>
#include <CGAL/IO/leda_window.h>
#include <CGAL/IO/Ostream_iterator.h>

#include <CGAL/convex_hull_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/IO/polygonal_2.h>
#include <CGAL/kernel_to_kernel.h>


int
main( int argc, char** argv)
{
  typedef CGAL::Cartesian<double>       C_double;
  typedef CGAL::Point_2< C_double>      double_Point;
  typedef CGAL::Segment_2< C_double>    double_Segment;
  typedef CGAL::Cartesian<leda_real>    C_real;
  typedef CGAL::Point_2< C_real>        real_Point;
  typedef CGAL::Segment_2< C_real>      real_Segment;
  typedef CGAL::Creator_uniform_2<double, double_Point>
                                        Point_creator;
  typedef CGAL::Random_points_in_square_2<double_Point, Point_creator>
                                        Source;
  typedef CGAL::Creator_uniform_2<double_Point,  double_Segment>
                                        Segment_creator;
  typedef CGAL::Join_input_iterator_2<Source, Source, Segment_creator>
                                        Segment_iterator;

  Source RS(280);
  Segment_iterator g( RS, RS);

  int N;
  if ( argc == 2)
  { N = CGAL_CLIB_STD::atoi( argv[1]); }
  else
  { std::cout << "How many segments? "; std::cin >> N; }

  std::vector< double_Segment>   double_segments;
  CGAL::copy_n( g, N, std::back_inserter( double_segments) );
  std::vector<real_Segment>  real_segments;
  CGAL::Cartesian_double_to_Cartesian<leda_real> converter;
  std::transform( double_segments.begin(),
                  double_segments.end(),
                  std::back_inserter( real_segments),
                  converter );
  std::vector<double_Point >   double_intersection_points;
  std::vector<real_Point >     real_intersection_points;

  CGAL::segment_intersection_points_2(
          double_segments.begin(),
          double_segments.end(),
          std::back_inserter( double_intersection_points),
          C_double() );
  CGAL::segment_intersection_points_2(
          real_segments.begin(),
          real_segments.end(),
          std::back_inserter( real_intersection_points),
          C_real() );


  std::vector<double_Point >   double_convex_hull;
  std::vector<real_Point >     real_convex_hull;

  typedef leda_window  CGAL_Stream;
  CGAL_Stream W0( 450, 500);
  CGAL_Stream W( 400, 400);
  CGAL::cgalize(W0);
  CGAL::cgalize(W);

  W.init( -300.0, 300.0, -300.0);
  W0.init( 0, 450, 0);
  W0.display();
  W.display(W0,25,25);
  W0.set_fg_color(leda_black);
  W0.draw_ctext(225,485,
     leda_string("Convex hull of intersection points"));
  W0.set_fg_color(leda_blue);
  W0 << double_Point(50,50);
  W0.set_fg_color(leda_red);
  W0 << double_Point(50,25);
  W0 << CGAL::Circle_2< C_double>( double_Point(50,25), 40.0);
  W0.set_fg_color(leda_black);
  W0.draw_text(65,58,
     leda_string("correct extreme point"));
  W0.set_fg_color(leda_black);
  W0.draw_text(65,35,
     leda_string("incorrect extreme point computed with double arithmetic"));

  CGAL::convex_hull_points_2(
          double_intersection_points.begin(),
          double_intersection_points.end(),
          std::back_inserter( double_convex_hull));
  W.set_fg_color( leda_grey2);
  CGAL::send_to_stream_as_polygon( W, double_convex_hull.begin(),
                                      double_convex_hull.end());
  W.set_fg_color( leda_green);
  std::copy( double_segments.begin(),
             double_segments.end(),
             CGAL::Ostream_iterator< double_Segment, CGAL_Stream>( W));
  CGAL::convex_hull_points_2(
          real_intersection_points.begin(),
          real_intersection_points.end(),
          std::back_inserter( real_convex_hull));
  W.set_fg_color( leda_blue);
  std::copy( double_convex_hull.begin(),
             double_convex_hull.end(),
             CGAL::Ostream_iterator< double_Point, CGAL_Stream>( W));

  if ( real_convex_hull.size() != double_convex_hull.size() )
  {
    W.set_fg_color( leda_red);
    std::copy( double_convex_hull.begin(),
               double_convex_hull.end(),
               CGAL::Ostream_iterator< double_Point, CGAL_Stream>( W));
    W.set_fg_color( leda_blue);
    std::copy( real_convex_hull.begin(),
               real_convex_hull.end(),
               CGAL::Ostream_iterator< real_Point, CGAL_Stream>( W));
    std::vector<double_Point >::iterator  dble_it;
    std::vector<real_Point >::iterator    real_it;
    real_it = real_convex_hull.begin();
    W.set_fg_color(leda_red);
    for ( dble_it  = double_convex_hull.begin();
          dble_it != double_convex_hull.end();
          ++dble_it )
    {
      if (   (real_it == real_convex_hull.end())
          || (( CGAL::squared_distance(
                      *dble_it,
                      double_Point( CGAL::to_double(real_it->x()),
                                    CGAL::to_double(real_it->y()) ))
               ) > 0.125 )
         )
      {
        W << CGAL::Circle_2< C_double>( *dble_it, 80.0);
      }
      else
      {
        if ( real_it != real_convex_hull.end() )
        { ++real_it; }
      }
    }
  }
  W.read_mouse();

  return 0;
}
#endif // USE_LEDA
