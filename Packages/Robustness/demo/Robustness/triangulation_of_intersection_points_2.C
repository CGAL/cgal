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
// file          : demo/Robustness/triangulation_of_intersection_points_2.C
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

#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/kernel_to_kernel.h>



// Workaround for VC++
#ifdef CGAL_CFG_MATCHING_BUG_2
#define CGAL_IA_CT double
#define CGAL_IA_PROTECTED true
#define CGAL_IA_ET leda_real
#define CGAL_IA_CACHE No_Filter_Cache
#endif
#include <CGAL/Arithmetic_filter.h>


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

  typedef CGAL::Cartesian<CGAL::Filtered_exact< double, leda_real> >
                                                           C_filtered;
  typedef CGAL::Point_2< C_filtered>                       C_filtered_Point;
  typedef CGAL::Segment_2< C_filtered>                     C_filtered_Segment;

  typedef CGAL::Triangulation_euclidean_traits_2<C_double> Gtd;
  typedef CGAL::Triangulation_vertex_base_2<Gtd>           Vbd;
  typedef CGAL::Triangulation_face_base_2<Gtd>             Fbd;
  typedef CGAL::Triangulation_default_data_structure_2<Gtd,Vbd,Fbd>
                                                           Tdsd;
  typedef CGAL::Delaunay_triangulation_2<Gtd,Tdsd>         DT_double;

  typedef CGAL::Triangulation_euclidean_traits_2<C_filtered> Gtf;
  typedef CGAL::Triangulation_vertex_base_2<Gtf>           Vbf;
  typedef CGAL::Triangulation_face_base_2<Gtf>             Fbf;
  typedef CGAL::Triangulation_default_data_structure_2<Gtf,Vbf,Fbf>
                                                           Tdsf;
  typedef CGAL::Delaunay_triangulation_2<Gtf,Tdsf>         DT_filtered;

  typedef CGAL::Triangulation_euclidean_traits_2<C_real>   Gtr;
  typedef CGAL::Triangulation_vertex_base_2<Gtr>           Vbr;
  typedef CGAL::Triangulation_face_base_2<Gtr>             Fbr;
  typedef CGAL::Triangulation_default_data_structure_2<Gtr,Vbr,Fbr>
                                                           Tdsr;
  typedef CGAL::Delaunay_triangulation_2<Gtr,Tdsr>         DT_real;
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



  std::vector<C_filtered_Segment >  filtered_segments;
  typedef CGAL::Filtered_exact< double, leda_real>  Filtered;
  CGAL::Cartesian_double_to_Cartesian<Filtered>  Fconverter;
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
  W0.draw_text(25,115, leda_string("Cartesian< leda_real >"));
  std::cout << "First with Cartesian< leda_real > ";
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
      leda_string("Cartesian< Filtered_exact< double, leda_real > >"));
  std::cout << "Next with Cartesian< Filtered_exact< double, leda_real> > ";
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
    int fistart = fi.pos(leda_string("include/CGAL/"),0);
    fistart += 13;
    W0.draw_text(25,85,
        leda_string("Cartesian< double >: ")
      + Msg.ty
      + leda_string(" violation")
      + leda_string(" in line ")
      + leda_string("%d", Msg.wo)
      + leda_string(" of file "));
    W0.draw_text(25,70, fi( fistart, 200 ));
    W0.draw_text(25,55, Msg.ex(0,60) );
    W0.draw_text(25,40, Msg.ex(60,120) );
    W0.draw_text(25,25, Msg.ex(120,180) );
    std::cout << Msg.ex << std::endl;
  }

  W.read_mouse();

  return 0;
}
#endif // USE_LEDA
