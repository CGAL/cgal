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
// file          : demo/Robustness/triangulation_of_intersection_points_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 
#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#include <vector>
#include <fstream>
#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_real.h>
typedef leda_real exact_NT;
#else
#  include <CGAL/Lazy_exact_nt.h>
#  include <CGAL/MP_Float.h>
#  include <CGAL/Quotient.h>
typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> >  exact_NT;
#endif
#include <CGAL/segment_intersection_points_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/copy_n.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Ostream_iterator.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/Cartesian_converter.h>

#include <CGAL/Filtered_kernel.h>

#if defined(CGAL_USE_CGAL_WINDOW)
#define leda_window  CGAL::window
#define leda_string  std::string
#define leda_green   CGAL::green
#define leda_blue    CGAL::blue
#define leda_black   CGAL::black
#define leda_red     CGAL::red
#endif


namespace CGAL {

class Assertion_violation
{
  public:
    Assertion_violation( const char*  type,
                         const char*  expression,
                         const char*  file,
                         int          line,
                         const char*  explanation)
      : ty( type), ex(expression), fi(file), wo(line)
    {
      if (explanation)  wh = leda_string(explanation);
      else wh = leda_string();
    }

    leda_string ty;
    leda_string ex;
    leda_string fi;
    leda_string wh;
    int         wo;
};


void
throw_exception_for_assertion_violation( const char*  type,
                                         const char*  expression,
                                         const char*  file,
                                         int          line,
                                         const char*  explanation)
{ throw Assertion_violation(type,expression,file,line,explanation); }

} // namespace CGAL


int
main( int argc, char** argv)
{
  CGAL::set_error_behaviour( CGAL::CONTINUE);
  CGAL::set_error_handler( CGAL::throw_exception_for_assertion_violation);

  typedef CGAL::Cartesian<double>       C_double;
  typedef CGAL::Cartesian<exact_NT>    C_real;
  typedef CGAL::Filtered_kernel<C_double, C_real> C_filtered;
  typedef C_double::Point_2             double_Point;
  typedef C_double::Segment_2           double_Segment;
  typedef C_real::Point_2               real_Point;
  typedef C_real::Segment_2             real_Segment;
  typedef C_filtered::Point_2           C_filtered_Point;
  typedef C_filtered::Segment_2         C_filtered_Segment;
  typedef CGAL::Creator_uniform_2<double, double_Point> Point_creator;
  typedef CGAL::Random_points_in_square_2<double_Point, Point_creator>
                                        Source;
  typedef CGAL::Creator_uniform_2<double_Point,  double_Segment>
                                        Segment_creator;
  typedef CGAL::Join_input_iterator_2<Source, Source, Segment_creator>
                                        Segment_iterator;

  typedef CGAL::Delaunay_triangulation_2<C_double>       DT_double;
  typedef CGAL::Delaunay_triangulation_2<C_filtered>     DT_filtered;
  typedef CGAL::Delaunay_triangulation_2<C_real>         DT_real;

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
  CGAL::Cartesian_converter<C_double, C_real> converter;
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



  std::vector<C_filtered_Segment >  filtered_segments;
  CGAL::Cartesian_converter<C_double, C_filtered>  Fconverter;
  std::transform( double_segments.begin(),
                  double_segments.end(),
                  std::back_inserter( filtered_segments),
                  Fconverter );
  std::vector<C_filtered_Point >  filtered_intersection_points;
  CGAL::segment_intersection_points_2(
          filtered_segments.begin(),
          filtered_segments.end(),
          std::back_inserter( filtered_intersection_points),
          C_filtered() );

  std::transform( double_segments.begin(),
                  double_segments.end(),
                  std::back_inserter( real_segments),
                  converter );


  typedef leda_window  CGAL_Stream;
  CGAL_Stream W0( 450, 550);
  CGAL_Stream W( 400, 400);
  CGAL::cgalize(W0);
  CGAL::cgalize(W);

  W.init( -300.0, 300.0, -300.0);
  W0.init( 0, 450, 0);
  W0.display();
  W.display(W0,25,25);
  W0.set_fg_color(leda_black);
  W0.draw_ctext(225,535,
     leda_string("Delaunay triangulation of intersection points"));

  W.set_fg_color( leda_green);
  std::copy( double_segments.begin(),
             double_segments.end(),
             CGAL::Ostream_iterator< double_Segment, CGAL_Stream>( W));
  W0.draw_text(25,115, leda_string("Cartesian< exact_NT >"));
  std::cout << "First with Cartesian< exact_NT > ";
  std::cout << std::endl;
  DT_real   DTR;
  std::copy( real_intersection_points.begin(),
             real_intersection_points.end(),
             std::back_inserter( DTR ));
  std::cout << std::distance( DTR.faces_begin(), DTR.faces_end() );
  std::cout << " triangles" << std::endl;
  W.set_fg_color(leda_blue);
  W << DTR;
  W.set_fg_color( leda_red);
  std::copy( double_intersection_points.begin(),
             double_intersection_points.end(),
             CGAL::Ostream_iterator< double_Point, CGAL_Stream>( W));
  W.read_mouse();

  W0.draw_text(25,100,
      leda_string("Filtered_kernel<Cartesian<double>, Cartesian<exact_NT> >"));
  std::cout << "Next with Filtered_kernel<Cartesian<double>, Cartesian<exact_NT> > ";
  std::cout << std::endl;
  DT_filtered FTD;
  std::copy( filtered_intersection_points.begin(),
             filtered_intersection_points.end(),
             std::back_inserter( FTD ));
  std::cout << std::distance( FTD.faces_begin(), FTD.faces_end() );
  std::cout << " triangles" << std::endl;
  W.read_mouse();

  W0.draw_text(25,85, leda_string("Cartesian< double >"));
  std::cout << "Next with Cartesian< double > ";
  std::cout << std::endl;
  try
  {
    DT_double DTD;
    std::copy( double_intersection_points.begin(),
               double_intersection_points.end(),
               std::back_inserter( DTD ));
    std::cout << std::distance( DTD.faces_begin(), DTD.faces_end() );
    std::cout << " triangles" << std::endl;
  }
  catch (const CGAL::Assertion_violation& Msg)
  {
    W0.set_fg_color( leda_red);
    leda_string fi = Msg.fi;
#ifdef CGAL_USE_CGAL_WINDOW
    int fistart = fi.find(leda_string("include/CGAL/"),0);
#else
    int fistart = fi.pos(leda_string("include/CGAL/"),0);
#endif
    fistart += 13;
    W0.draw_text(25,85,
        leda_string("Cartesian< double >: ")
      + Msg.ty
      + leda_string(" violation")
      + leda_string(" in line ")
      + leda_string("%d", Msg.wo)
      + leda_string(" of file "));
#ifdef CGAL_USE_CGAL_WINDOW
    W0.draw_text(25,70, fi.substr( fistart, 200 ));
    W0.draw_text(25,55, Msg.ex.substr(0,60) );
    W0.draw_text(25,40, Msg.ex.substr(60,120) );
    W0.draw_text(25,25, Msg.ex.substr(120,180) );
#else
    W0.draw_text(25,70, fi( fistart, 200 ));
    W0.draw_text(25,55, Msg.ex(0,60) );
    W0.draw_text(25,40, Msg.ex(60,120) );
    W0.draw_text(25,25, Msg.ex(120,180) );
#endif
    std::cout << Msg.ex << std::endl;
  }

  W.read_mouse();

  return 0;
}
