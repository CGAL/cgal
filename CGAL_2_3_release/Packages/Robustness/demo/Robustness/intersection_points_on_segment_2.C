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
// file          : demo/Robustness/intersection_points_on_segment_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <cassert>
#include <vector>
#include <algorithm>
#include <CGAL/intersection_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/copy_n.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/kernel_to_kernel.h>

#if defined(CGAL_USE_CGAL_WINDOW)
#define leda_window  CGAL::window
#define leda_string  std::string
#endif

#include <CGAL/intersection_test_statistics.h>
#include <CGAL/further_point_generators_2.h>

using namespace CGAL;

typedef Cartesian<float>                 CartesianFloat;
typedef Cartesian<double>                CartesianDouble;
typedef Homogeneous<float>               HomogeneousFloat;
typedef Homogeneous<double>              HomogeneousDouble;

typedef CartesianDouble::Point_2         Point;
typedef CartesianDouble::Segment_2       Segment;
typedef std::vector<Segment>             Vector;

int
main(int argc, char** argv)
{
  int N = (argc > 1) ? CGAL_CLIB_STD::atoi(argv[1]) : 50;

  typedef Creator_uniform_2< Point, Segment>                Seg_creator;
  typedef Join_input_iterator_2< Point_on_vertical_bar,
                                 Point_on_vertical_bar,
                                 Seg_creator>               Seg1_iterator;
  typedef Join_input_iterator_2< Point_on_horizontal_bar,
                                 Point_on_horizontal_bar,
                                 Seg_creator>               Seg2_iterator;
  Point_on_vertical_bar p1( -100000, 100000, -120000);
  Point_on_vertical_bar p2( -100000, 100000,  120000);
  Point_on_horizontal_bar p3( -100000, 100000, -120000);
  Point_on_horizontal_bar p4( -100000, 100000,  120000);
  Seg1_iterator g1( p1, p2);
  Seg2_iterator g2( p3, p4);
  Vector CD1; // vertical   segments
  Vector CD2; // horizontal segments
  CGAL::copy_n( g1, N, std::back_inserter(CD1));
  CGAL::copy_n( g2, N, std::back_inserter(CD2));

  leda_window W(1000,550);
  leda_window W1(250,250);
  leda_window W2(250,250);
  leda_window W3(250,250);
  leda_window W4(250,250);
  CGAL::cgalize(W);
  CGAL::cgalize(W1);
  CGAL::cgalize(W2);
  CGAL::cgalize(W3);
  CGAL::cgalize(W4);
  W.init(0,1000,0);
  W1.init(-133000,133000,-133000);
  W2.init(-133000,133000,-133000);
  W3.init(-133000,133000,-133000);
  W4.init(-133000,133000,-133000);
  W1.set_node_width(2);
  W2.set_node_width(2);
  W3.set_node_width(2);
  W4.set_node_width(2);
  W.display();
  W1.display(W,   0, 170);
  W2.display(W, 250, 170);
  W3.display(W, 500, 170);
  W4.display(W, 750, 170);
  W1 << CGAL::GREEN;  W2 << CGAL::GREEN;  W3 << CGAL::GREEN; W4 << GREEN;
  Vector::iterator s;
  for( s = CD1.begin(); s != CD1.end(); s++)
  { W1 << *s; W2 << *s; W3 << *s; W4 << *s; }
  for( s = CD2.begin(); s != CD2.end(); s++)
  { W1 << *s; W2 << *s; W3 << *s; W4 << *s; }

  std::vector< CartesianFloat::Segment_2>  CF1;
  std::vector< CartesianFloat::Segment_2>  CF2;
  CGAL::Cartesian_converter<CartesianDouble, CartesianFloat> converter1;
  std::transform( CD1.begin(), CD1.end(), std::back_inserter( CF1), converter1);
  std::transform( CD2.begin(), CD2.end(), std::back_inserter( CF2), converter1);
  std::vector< HomogeneousFloat::Segment_2>  HF1;
  std::vector< HomogeneousFloat::Segment_2>  HF2;
  Cartesian_double_to_Homogeneous< HomogeneousFloat::RT > converter2;
  std::transform( CD1.begin(), CD1.end(), std::back_inserter( HF1), converter2);
  std::transform( CD2.begin(), CD2.end(), std::back_inserter( HF2), converter2);
  std::vector< HomogeneousDouble::Segment_2>  HD1;
  std::vector< HomogeneousDouble::Segment_2>  HD2;
  Cartesian_double_to_Homogeneous< HomogeneousDouble::RT > converter3;
  std::transform( CD1.begin(), CD1.end(), std::back_inserter( HD1), converter3);
  std::transform( CD2.begin(), CD2.end(), std::back_inserter( HD2), converter3);

  leda_string str;

  leda_string    header("Two sets of segments are generated at random.");
  header += leda_string("All endpoint coordinates are integral and have ");
  header += leda_string("absolute value at most 100 000.\n");
  header += leda_string("The endpoints of the first set of segments lie on ");
  header += leda_string("two vertical bars; in the second set, the segments ");
  header += leda_string("lie on two horizontal bars.\n");
  header += leda_string("Each segment in the first set is intersected with ");
  header += leda_string("each segment in the second set. Then, it is checked ");
  header += leda_string("whether the intersection points\n ");
  header += leda_string("lie on the intersecting segments. ");
  header += leda_string("Computation is done with kernels ");
  header += leda_string("differing in representation and float precision.");
  W.draw_text( 5, 540, header);
  Point o_p(300,450);
  W << CGAL::ORANGE;
  W << o_p;
  Point r_p(300,435);
  W << CGAL::RED;
  W << r_p;
  W << CGAL::BLACK;
  leda_string
      orange("INCORRECT intersection point: lies on one of the lines only.");
  leda_string
      red("INCORRECT intersection point: lies on neither of the two lines.");
  W.draw_text( 310, 458, orange);
  W.draw_text( 310, 443, red);

  W.draw_ctext( 125, 400, leda_string("Cartesian<float>"));
  intersection_statistics(CF1.begin(), CF1.end(),
                          CF2.begin(), CF2.end(),
                          str, CartesianFloat());
  draw_errors(W1, CF1.begin(), CF1.end(), CF2.begin(), CF2.end(),
              CartesianFloat());
  W.draw_text( 5, 110, str);

  W.draw_ctext( 375, 400, leda_string("Cartesian<double>"));
  intersection_statistics(CD1.begin(), CD1.end(),
                          CD2.begin(), CD2.end(),
                          str, CartesianDouble());
  draw_errors(W2, CD1.begin(), CD1.end(), CD2.begin(), CD2.end(),
              CartesianDouble());
  W.draw_text( 255, 110, str);

  W.draw_ctext( 625, 400, leda_string("Homogeneous<float>"));
  intersection_statistics(HF1.begin(), HF1.end(),
                          HF2.begin(), HF2.end(),
                          str, HomogeneousFloat());
  draw_errors(W3, HF1.begin(), HF1.end(), HF2.begin(), HF2.end(),
              HomogeneousFloat());
  W.draw_text( 505, 110, str);

  W.draw_ctext( 875, 400, leda_string("Homogeneous<double>"));
  intersection_statistics(HD1.begin(), HD1.end(),
                          HD2.begin(), HD2.end(),
                          str, HomogeneousDouble());
  draw_errors(W4, HD1.begin(), HD1.end(), HD2.begin(), HD2.end(),
              HomogeneousDouble());
  W.draw_text( 755, 110, str);

  W.read_mouse();
  return 0;
}
