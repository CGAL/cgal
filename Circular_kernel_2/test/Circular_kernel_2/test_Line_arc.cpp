// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#include <CGAL/Cartesian.h>
#include <cassert>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Circular_kernel_2/Circular_arc_2.h>
#include <CGAL/Circular_kernel_2/Line_arc_2.h>
#include <CGAL/Exact_circular_kernel_2.h>

#include <CGAL/Random.h>

template <class Expected, class V>
bool assign_variant(Expected& e, const V& v)
{
  if (std::get_if<Expected>(&v) != nullptr)
  {
    e = std::get<Expected>(v);
    return true;
  }
  return false;
}

template <class CK>
void _test_Line_arc(CK ck)
{

  typedef CGAL::Circle_2<CK>                   Circle_2;
  typedef CGAL::Circular_arc_2<CK>             Circular_arc_2;
  typedef CGAL::Point_2<CK>                    Point_2;
  typedef CGAL::Line_2<CK>                     Line_2;
  typedef typename CK::Intersect_2   Intersect_2;
  typedef typename CK::Split_2                 Split_2;
  typedef typename CK::Circular_arc_point_2 Circular_arc_point_2;
  typedef typename CK::Line_arc_2              Line_arc_2;
  typedef typename CK::Do_overlap_2            Do_overlap_2;
  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);
  //int random_max = 127;
  //int random_min = -127;
  Point_2 center_circle1(0, 0);
  int circle1_r = 5;
  Circle_2 circle1(center_circle1, circle1_r * circle1_r);
  Circle_2 circle2(center_circle1, circle1_r * circle1_r * 4);
  Point_2 p2_line_horizontal(1, 0);

  Line_arc_2 line_arc_horizontal = Line_arc_2(Line_2(center_circle1, p2_line_horizontal),
                        circle1,
                        true,
                        circle2,
                        false);
  assert(line_arc_horizontal.source().x() == -5);
  assert(line_arc_horizontal.source().y() == 0);
  assert(line_arc_horizontal.target().x() == 10);
  assert(line_arc_horizontal.target().y() == 0);
  Point_2 p2_line_vertical(0, 1);
  Line_arc_2 line_arc_vertical(Line_2(center_circle1, p2_line_vertical),
                        circle1,
                        true,
                        Circle_2(center_circle1, circle1_r * circle1_r * 4),
                        false);
  assert(line_arc_vertical.source().x() == 0);
  assert(line_arc_vertical.source().y() == -5);
  assert(line_arc_vertical.target().x() == 0);
  assert(line_arc_vertical.target().y() == 10);
  Point_2 p2_line_diagonal(3, 4);
  Line_arc_2 line_arc_diagonal(Line_2(center_circle1, p2_line_diagonal),
                        circle1,
                        true,
                        Circle_2(center_circle1, circle1_r * circle1_r * 4),
                        false);
  Circular_arc_2(circle1, line_arc_diagonal.source(), line_arc_diagonal.source());
  Circular_arc_2(circle2, line_arc_diagonal.target(), line_arc_diagonal.target());

   Line_arc_2 line_arc_2(Line_2(center_circle1, p2_line_horizontal),
                        Line_2(Point_2(0, -1), Point_2(1,0)),
                        Line_2(Point_2(0, -1), Point_2(-1, 0)));

  //////////////Intersection Line_arc Line_arc//////////////////

   Intersect_2 theConstruct_intersect_2 = ck.intersect_2_object();
   typedef std::variant<std::pair<Circular_arc_point_2, unsigned int>, Line_arc_2> Intersection_result;

   std::vector< Intersection_result >
     vector_for_intersection_1;
   theConstruct_intersect_2(line_arc_horizontal,
                            line_arc_vertical,
                            std::back_inserter(vector_for_intersection_1));
   std::pair<Circular_arc_point_2, unsigned int> the_pair;
   assert(assign_variant(the_pair, vector_for_intersection_1[0]));
   Circular_arc_point_2 first = the_pair.first;
   assert(first.x() == 0);
   assert(first.y() == 0);

   Line_arc_2 line_arc_3(Point_2(-1, -2), Point_2(2,1));
   std::vector< Intersection_result >
     vector_for_intersection_2;
   theConstruct_intersect_2(line_arc_horizontal,
                            line_arc_3,
                            std::back_inserter(vector_for_intersection_2));
   assert(assign_variant(the_pair, vector_for_intersection_2[0]));
   first = the_pair.first;
   assert(first.x() == 1);
   assert(first.y() == 0);

   std::vector< Intersection_result >
     vector_for_intersection_3;
   theConstruct_intersect_2(line_arc_vertical,
                            line_arc_3,
                            std::back_inserter(vector_for_intersection_3));
   assert(assign_variant(the_pair, vector_for_intersection_3[0]));
   first = the_pair.first;
   assert(first.x() == 0);
   assert(first.y() == -1);

   Line_arc_2 line_arc_4(Point_2(20, -2), Point_2(20,10));
   std::vector< Intersection_result >
     vector_for_intersection_4;
   theConstruct_intersect_2(line_arc_horizontal,
                            line_arc_4,
                            std::back_inserter(vector_for_intersection_4));
   assert(vector_for_intersection_4.size() == 0);



   ////////////intersection in overlap///////////////
   std::vector< Intersection_result >
     vector_for_intersection_5;
   theConstruct_intersect_2(line_arc_horizontal,
                            line_arc_horizontal,
                            std::back_inserter(vector_for_intersection_5));
   Line_arc_2 line_arc_tmp;
   assert(vector_for_intersection_5.size() == 1);
   assert(assign_variant(line_arc_tmp, vector_for_intersection_5[0]));
   assert(line_arc_tmp == line_arc_horizontal);

   Line_arc_2 line_arc_horizontal2(Line_2(center_circle1, p2_line_horizontal),
                        circle1,
                        true,
                        circle1,
                        false);
   std::vector< Intersection_result >
     vector_for_intersection_6;
   theConstruct_intersect_2(line_arc_horizontal,
                            line_arc_horizontal2,
                            std::back_inserter(vector_for_intersection_6));
   assert(vector_for_intersection_6.size() == 1);
   assert(assign_variant(line_arc_tmp, vector_for_intersection_6[0]));
   assert(line_arc_tmp.source() == line_arc_horizontal.source());
   assert(line_arc_tmp.target() == line_arc_horizontal2.target());



   std::vector< Intersection_result >
     vector_for_intersection_6_reverse;
   theConstruct_intersect_2(line_arc_horizontal2,
                            line_arc_horizontal,
                            std::back_inserter(vector_for_intersection_6_reverse));
   assert(vector_for_intersection_6_reverse.size() == 1);
   assert(assign_variant(line_arc_tmp, vector_for_intersection_6_reverse[0]));
   assert(line_arc_tmp.source() == line_arc_horizontal.source());
   assert(line_arc_tmp.target() == line_arc_horizontal2.target());



   Line_arc_2 line_arc_horizontal3(Line_2(center_circle1, p2_line_horizontal),
                        circle2,
                        true,
                        circle1,
                        false);
   std::vector< Intersection_result >
     vector_for_intersection_7;
   theConstruct_intersect_2(line_arc_horizontal,
                            line_arc_horizontal3,
                            std::back_inserter(vector_for_intersection_7));
   assert(vector_for_intersection_7.size() == 1);
   assert(assign_variant(line_arc_tmp, vector_for_intersection_7[0]));
   assert(line_arc_tmp.source() == line_arc_horizontal.source());
   assert(line_arc_tmp.target() == line_arc_horizontal3.target());

   std::vector< Intersection_result >
     vector_for_intersection_8;
   theConstruct_intersect_2(line_arc_horizontal3,
                            line_arc_horizontal,
                            std::back_inserter(vector_for_intersection_8));
   assert(vector_for_intersection_8.size() == 1);
   assert(assign_variant(line_arc_tmp, vector_for_intersection_8[0]));
   assert(line_arc_tmp.source() == line_arc_horizontal.source());
   assert(line_arc_tmp.target() == line_arc_horizontal3.target());

   std::vector< Intersection_result >
     vector_for_intersection_8_bis_1;
   theConstruct_intersect_2(line_arc_horizontal3,
                            line_arc_horizontal2,
                            std::back_inserter(vector_for_intersection_8_bis_1));
   assert(vector_for_intersection_8_bis_1.size() == 1);
   assert(assign_variant(line_arc_tmp, vector_for_intersection_8_bis_1[0]));
   assert(line_arc_tmp.source() == line_arc_horizontal2.source());
   assert(line_arc_tmp.target() == line_arc_horizontal3.target());

   std::vector< Intersection_result >
     vector_for_intersection_8_bis_2;
   theConstruct_intersect_2(line_arc_horizontal2,
                            line_arc_horizontal3,
                            std::back_inserter(vector_for_intersection_8_bis_2));
   assert(vector_for_intersection_8_bis_2.size() == 1);
   assert(assign_variant(line_arc_tmp, vector_for_intersection_8_bis_2[0]));
   assert(line_arc_tmp.source() == line_arc_horizontal2.source());
   assert(line_arc_tmp.target() == line_arc_horizontal3.target());


   Line_arc_2 line_arc_horizontal4(Line_2(center_circle1, p2_line_horizontal),
                                   circle2,
                                   true,
                                   circle1,
                                   true);
   std::vector< Intersection_result >
     vector_for_intersection_9;
   theConstruct_intersect_2(line_arc_horizontal4,
                            line_arc_horizontal,
                            std::back_inserter(vector_for_intersection_9));
   assert(vector_for_intersection_9.size() == 1);
   assert(assign_variant(the_pair, vector_for_intersection_9[0]));
   assert(the_pair.second == 1);
   assert(the_pair.first == line_arc_horizontal.source());

    std::vector< Intersection_result >
     vector_for_intersection_10;
   theConstruct_intersect_2(line_arc_horizontal,
                            line_arc_horizontal4,
                            std::back_inserter(vector_for_intersection_10));
   assert(vector_for_intersection_10.size() == 1);
   assert(assign_variant(the_pair, vector_for_intersection_10[0]));
   assert(the_pair.second == 1);
   assert(the_pair.first == line_arc_horizontal.source());


   std::vector< Intersection_result >
     vector_for_intersection_11;
   theConstruct_intersect_2(Line_arc_2(Point_2(-1, -1),Point_2(2, 2)),
                            Line_arc_2(Point_2(0, 0), Point_2(1, 1)),
                            std::back_inserter(vector_for_intersection_11));
   assert(vector_for_intersection_11.size() == 1);
   assert(assign_variant(line_arc_tmp, vector_for_intersection_11[0]));
   assert(line_arc_tmp == Line_arc_2(Point_2(0, 0), Point_2(1, 1)));


   //////Split/////////////

   std::vector< Intersection_result >
     vector_for_intersection_split_1;
   theConstruct_intersect_2(line_arc_horizontal,
                            line_arc_3,
                            std::back_inserter(vector_for_intersection_split_1));
   assert(assign_variant(the_pair, vector_for_intersection_split_1[0]));
   Circular_arc_point_2 point_of_split = the_pair.first;
   Split_2 theSplit_2 = ck.split_2_object();
   Line_arc_2 first_part_line_arc_horizontal;
   Line_arc_2 second_part_line_arc_horizontal;
   theSplit_2(line_arc_horizontal, point_of_split, first_part_line_arc_horizontal, second_part_line_arc_horizontal);
   assert(first_part_line_arc_horizontal.source() == line_arc_horizontal.source());
   assert(first_part_line_arc_horizontal.target() == point_of_split);
   assert(second_part_line_arc_horizontal.source() == point_of_split);
   assert(second_part_line_arc_horizontal.target() == line_arc_horizontal.target());

   ////Overlap//////

   Do_overlap_2 theDo_overlap_2 = ck.do_overlap_2_object();
   assert(theDo_overlap_2(line_arc_horizontal, line_arc_horizontal));
   assert(theDo_overlap_2(line_arc_horizontal, line_arc_horizontal2));
   assert(!theDo_overlap_2(line_arc_horizontal, line_arc_vertical));

}


template <class CK>
void _test_intersection_Line_arc_Circle(CK ck)
{

  typedef CGAL::Circle_2<CK>                   Circle_2;

  typedef CGAL::Point_2<CK>                    Point_2;
  typedef CGAL::Line_2<CK>                     Line_2;
  typedef CGAL::Line_arc_2<CK>                 Line_arc_2;
  typedef CGAL::Circular_arc_point_2<CK>       Circular_arc_point_2;
  typedef typename CK::Intersect_2   Intersect_2;
  typedef std::variant<std::pair<Circular_arc_point_2, unsigned int>, Line_arc_2> Intersection_result;

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);

  Intersect_2 theConstruct_intersect_2 = ck.intersect_2_object();
  Point_2 center_circle1(0, 0);
  int circle1_r = 5;
  Circle_2 circle1(center_circle1, circle1_r * circle1_r);
  Circle_2 circle2(center_circle1, circle1_r * circle1_r * 4);
  Point_2 p2_line_horizontal(1, 0);
  Point_2 p2_high_circle2(0, circle1_r * 2);
  Point_2 p2_low_circle2(0, circle1_r * -2);
  Point_2 p2_left_circle2(circle1_r * -2, 0);
  Point_2 p2_right_circle2(circle1_r * 2, 0);
  Point_2 p2_high_right_circle2(circle1_r * 2, circle1_r * 2);
  Point_2 p2_high_left_circle2(circle1_r * -2, circle1_r * 2);
  Point_2 p2_low_right_circle2(circle1_r * 2, circle1_r * -2);
  Point_2 p2_low_left_circle2(circle1_r * -2, circle1_r * -2);
  Line_arc_2 line_arc_horizontal(Line_2(center_circle1,
                                        p2_line_horizontal),
                                 circle1,
                                 true,
                                 circle2,
                                 false);
  Line_arc_2 line_arc_horizontal2(Line_2(center_circle1,
                                         p2_line_horizontal),
                                  circle1,
                                  true,
                                  circle1,
                                  false);
  Line_arc_2 line_arc_horizontal3(Line_2(center_circle1,
                                        p2_line_horizontal),
                                 circle2,
                                 true,
                                 circle1,
                                 false);
  Point_2 p2_line_vertical(0, 1);
  Line_arc_2 line_arc_vertical(Line_2(center_circle1,
                                      p2_line_vertical),
                               circle1,
                               true,
                               circle2,
                               false);
  Line_arc_2 line_arc_vertical2(Line_2(center_circle1,
                                       p2_line_vertical),
                                circle1,
                                true,
                                circle1,
                                false);
  Line_arc_2 line_arc_vertical3(Line_2(center_circle1,
                                       p2_line_vertical),
                                circle2,
                                true,
                                circle1,
                                false);
  Point_2 p2_line_diagonal(2, -2);
  Line_arc_2 line_arc_diagonal(Line_2(center_circle1,
                                      p2_line_diagonal),
                               circle1,
                               true,
                               circle2,
                               false);
   Line_arc_2 line_arc_diagonal2(Line_2(center_circle1,
                                       p2_line_diagonal),
                                circle1,
                                true,
                                circle1,
                                false);
  Line_arc_2 line_arc_diagonal3(Line_2(center_circle1,
                                       p2_line_diagonal),
                                circle2,
                                true,
                                circle1,
                                false);

  ////////test horizontal//////
  std::vector< Intersection_result >
    vector_for_intersection_1;
  theConstruct_intersect_2(line_arc_horizontal,
                           circle1,
                           std::back_inserter(vector_for_intersection_1));
  std::pair<Circular_arc_point_2, unsigned int> the_pair;
  assert(vector_for_intersection_1.size() == 2);
  assert(assign_variant(the_pair, vector_for_intersection_1[0]));
  Circular_arc_point_2 first = the_pair.first;
  assert(first == line_arc_horizontal.source());
  assert(the_pair.second == 1);
  assert(assign_variant(the_pair, vector_for_intersection_1[1]));
  first = the_pair.first;
  assert(first == line_arc_horizontal2.target());
  assert(the_pair.second == 1);

  std::vector< Intersection_result >
    vector_for_intersection_2;
  theConstruct_intersect_2(line_arc_horizontal,
                           circle2,
                           std::back_inserter(vector_for_intersection_2));
  assert(vector_for_intersection_2.size() == 1);
  assert(assign_variant(the_pair, vector_for_intersection_2[0]));
  first = the_pair.first;
  assert(first == line_arc_horizontal.target());
  assert(the_pair.second == 1);

  std::vector< Intersection_result >
    vector_for_intersection_3;
  theConstruct_intersect_2(line_arc_horizontal2,
                           circle2,
                           std::back_inserter(vector_for_intersection_3));
  assert(vector_for_intersection_3.size() == 0);

  std::vector< Intersection_result >
    vector_for_intersection_4;
  theConstruct_intersect_2(line_arc_horizontal3,
                           circle2,
                           std::back_inserter(vector_for_intersection_4));
  assert(vector_for_intersection_4.size() == 1);
  assert(assign_variant(the_pair, vector_for_intersection_4[0]));
  first = the_pair.first;
  assert(first == line_arc_horizontal3.source());
  assert(the_pair.second == 1);

  Line_arc_2 line_arc_aux(p2_high_circle2,p2_high_right_circle2);
  std::vector< Intersection_result >
    vector_for_intersection_5;
  theConstruct_intersect_2(Line_arc_2(p2_high_left_circle2,p2_high_right_circle2),
                           circle2,
                           std::back_inserter(vector_for_intersection_5));
  assert(vector_for_intersection_5.size() == 1);
  theConstruct_intersect_2(line_arc_aux,
                           circle2,
                           std::back_inserter(vector_for_intersection_5));
  assert(vector_for_intersection_5.size() == 2);
  assert(assign_variant(the_pair, vector_for_intersection_5[0]));
  first = the_pair.first;
  assert(first == line_arc_aux.source());
  assert(the_pair.second == 2);
  assert(assign_variant(the_pair, vector_for_intersection_5[1]));
  first = the_pair.first;
  assert(first == line_arc_aux.source());
  assert(the_pair.second == 2);

  line_arc_aux = Line_arc_2(p2_low_circle2,p2_low_right_circle2);
  std::vector< Intersection_result >
    vector_for_intersection_6;
  theConstruct_intersect_2(Line_arc_2(p2_low_left_circle2,p2_low_right_circle2),
                           circle2,
                           std::back_inserter(vector_for_intersection_6));
  assert(vector_for_intersection_6.size() == 1);
  theConstruct_intersect_2(line_arc_aux,
                           circle2,
                           std::back_inserter(vector_for_intersection_6));
  assert(vector_for_intersection_6.size() == 2);
  assert(assign_variant(the_pair, vector_for_intersection_6[0]));
  first = the_pair.first;
  assert(first == line_arc_aux.source());
  assert(the_pair.second == 2);
  assert(assign_variant(the_pair, vector_for_intersection_6[1]));
  first = the_pair.first;
  assert(first == line_arc_aux.source());
  assert(the_pair.second == 2);

  ////////test vertical//////
  std::vector< Intersection_result >
    vector_for_intersection_7;
  theConstruct_intersect_2(line_arc_vertical,
                           circle1,
                           std::back_inserter(vector_for_intersection_7));
  assert(vector_for_intersection_7.size() == 2);
  assert(assign_variant(the_pair, vector_for_intersection_7[0]));
  first = the_pair.first;
  assert(first == line_arc_vertical.source());
  assert(the_pair.second == 1);
  assert(assign_variant(the_pair, vector_for_intersection_7[1]));
  first = the_pair.first;
  assert(first == line_arc_vertical2.target());
  assert(the_pair.second == 1);

  std::vector< Intersection_result >
    vector_for_intersection_8;
  theConstruct_intersect_2(line_arc_vertical,
                           circle2,
                           std::back_inserter(vector_for_intersection_8));
  std::cout << vector_for_intersection_8.size() << std::endl;
  assert(vector_for_intersection_8.size() == 1);
  assert(assign_variant(the_pair, vector_for_intersection_8[0]));
  first = the_pair.first;
  assert(first == line_arc_vertical.target());
  assert(the_pair.second == 1);

  std::vector< Intersection_result >
    vector_for_intersection_9;
  theConstruct_intersect_2(line_arc_vertical2,
                           circle2,
                           std::back_inserter(vector_for_intersection_9));
  assert(vector_for_intersection_9.size() == 0);

  std::vector< Intersection_result >
    vector_for_intersection_10;
  theConstruct_intersect_2(line_arc_vertical3,
                           circle2,
                           std::back_inserter(vector_for_intersection_10));
  assert(vector_for_intersection_10.size() == 1);
  assert(assign_variant(the_pair, vector_for_intersection_10[0]));
  first = the_pair.first;
  assert(first == line_arc_vertical3.source());
  assert(the_pair.second == 1);

  line_arc_aux = Line_arc_2(p2_right_circle2,p2_high_right_circle2);
  std::vector< Intersection_result >
    vector_for_intersection_11;
  theConstruct_intersect_2(Line_arc_2(p2_low_right_circle2,p2_high_right_circle2),
                            circle2,
                            std::back_inserter(vector_for_intersection_11));
  assert(vector_for_intersection_11.size() == 1);
  theConstruct_intersect_2(line_arc_aux,
                            circle2,
                            std::back_inserter(vector_for_intersection_11));
  assert(vector_for_intersection_11.size() == 2);
  assert(assign_variant(the_pair, vector_for_intersection_11[0]));
  first = the_pair.first;
  assert(first == line_arc_aux.source());
  assert(the_pair.second == 2);
  assert(assign_variant(the_pair, vector_for_intersection_11[1]));
  first = the_pair.first;
  assert(first == line_arc_aux.source());
  assert(the_pair.second == 2);

  line_arc_aux = Line_arc_2(p2_left_circle2,p2_high_left_circle2);
  std::vector< Intersection_result >
    vector_for_intersection_12;
  theConstruct_intersect_2(Line_arc_2(p2_low_left_circle2,p2_high_left_circle2),
                           circle2,
                           std::back_inserter(vector_for_intersection_12));
  assert(vector_for_intersection_12.size() == 1);
  theConstruct_intersect_2(line_arc_aux,
                           circle2,
                           std::back_inserter(vector_for_intersection_12));
  assert(vector_for_intersection_12.size() == 2);
  assert(assign_variant(the_pair, vector_for_intersection_12[0]));
  first = the_pair.first;
  assert(first == line_arc_aux.source());
  assert(the_pair.second == 2);
  assert(assign_variant(the_pair, vector_for_intersection_12[1]));
  first = the_pair.first;
  assert(first == line_arc_aux.source());
  assert(the_pair.second == 2);

  //test diagonal
  std::vector< Intersection_result >
    vector_for_intersection_13;
  theConstruct_intersect_2(line_arc_diagonal,
                           circle1,
                           std::back_inserter(vector_for_intersection_13));
  assert(vector_for_intersection_13.size() == 2);
  assert(assign_variant(the_pair, vector_for_intersection_13[0]));
  first = the_pair.first;
  assert(first == line_arc_diagonal.source());
  assert(the_pair.second == 1);
  assert(assign_variant(the_pair, vector_for_intersection_13[1]));
  first = the_pair.first;
  assert(first == line_arc_diagonal2.target());
  assert(the_pair.second == 1);

  std::vector< Intersection_result >
    vector_for_intersection_14;
  theConstruct_intersect_2(line_arc_diagonal,
                           circle2,
                           std::back_inserter(vector_for_intersection_14));
  assert(vector_for_intersection_14.size() == 1);
  assert(assign_variant(the_pair, vector_for_intersection_14[0]));
  first = the_pair.first;
  assert(first == line_arc_diagonal.target());
  assert(the_pair.second == 1);

  std::vector< Intersection_result >
    vector_for_intersection_15;
  theConstruct_intersect_2(line_arc_diagonal2,
                           circle2,
                           std::back_inserter(vector_for_intersection_15));
  assert(vector_for_intersection_15.size() == 0);

  std::vector< Intersection_result >
    vector_for_intersection_16;
  theConstruct_intersect_2(line_arc_diagonal3,
                           circle2,
                           std::back_inserter(vector_for_intersection_16));
  assert(vector_for_intersection_16.size() == 1);
  assert(assign_variant(the_pair, vector_for_intersection_16[0]));
  first = the_pair.first;
  assert(first == line_arc_diagonal3.source());
  assert(the_pair.second == 1);


  //Diagonal tangent

   std::vector< Intersection_result >
     vector_for_intersection_17;
   theConstruct_intersect_2(Line_arc_2(Point_2(-10, -5), Point_2(5, 15)),
                             circle1,
                             std::back_inserter(vector_for_intersection_17));
   assert(vector_for_intersection_17.size() == 1);
   assert(assign_variant(the_pair, vector_for_intersection_17[0]));
   assert(the_pair.second == 2);


   std::vector< Intersection_result >
     vector_for_intersection_18;
   theConstruct_intersect_2(Line_arc_2(Point_2(10, -5), Point_2(-5, 15)),
                             circle1,
                             std::back_inserter(vector_for_intersection_18));
   assert(vector_for_intersection_18.size() == 1);
   assert(assign_variant(the_pair, vector_for_intersection_18[0]));
   assert(the_pair.second == 2);


   std::vector< Intersection_result >
     vector_for_intersection_19;
   theConstruct_intersect_2(Line_arc_2(Point_2(10, 5), Point_2(-5, -15)),
                             circle1,
                             std::back_inserter(vector_for_intersection_19));
   assert(vector_for_intersection_19.size() == 1);
   assert(assign_variant(the_pair, vector_for_intersection_19[0]));
   assert(the_pair.second == 2);

   std::vector< Intersection_result >
     vector_for_intersection_20;
   theConstruct_intersect_2(Line_arc_2(Point_2(-10, 5), Point_2(5, -15)),
                             circle1,
                             std::back_inserter(vector_for_intersection_20));
   assert(vector_for_intersection_20.size() == 1);
   assert(assign_variant(the_pair, vector_for_intersection_20[0]));
   assert(the_pair.second == 2);

}

template <class CK>
void _test_intersection_Line_arc_Circular_arc(CK ck)
{

  typedef CGAL::Circle_2<CK>                   Circle_2;
  typedef CGAL::Circular_arc_2<CK>             Circular_arc_2;
  typedef CGAL::Point_2<CK>                    Point_2;
  typedef CGAL::Line_2<CK>                     Line_2;
  typedef CGAL::Line_arc_2<CK>                 Line_arc_2;
  typedef CGAL::Circular_arc_point_2<CK>       Circular_arc_point_2;
  typedef typename CK::Intersect_2   Intersect_2;
  typedef std::variant<std::pair<Circular_arc_point_2, unsigned int>, Circular_arc_2> Intersection_result;
  typedef typename CK::Make_x_monotone_2           Make_x_monotone_2;
  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);

  Intersect_2 theConstruct_intersect_2 = ck.intersect_2_object();
  Point_2 center_circle1(0, 0);
  int circle1_r = 5;
  Circle_2 circle1(center_circle1, circle1_r * circle1_r);
  Circle_2 circle2(center_circle1, circle1_r * circle1_r * 4);
  Point_2 p2_line_horizontal(1, 0);

  Line_arc_2 line_arc_horizontal(Line_2(center_circle1,
                                        p2_line_horizontal),
                                 circle1,
                                 true,
                                 circle2,
                                 false);
  Line_arc_2 line_arc_horizontal2(Line_2(center_circle1,
                                         p2_line_horizontal),
                                  circle1,
                                  true,
                                  circle1,
                                  false);
  Line_arc_2 line_arc_horizontal3(Line_2(center_circle1,
                                         p2_line_horizontal),
                                  circle2,
                                  true,
                                  circle2,
                                  false);
  Point_2 p2_line_vertical(0, 1);
  Line_arc_2 line_arc_vertical(Line_2(center_circle1,
                                      p2_line_vertical),
                               circle1,
                               true,
                               Circle_2(center_circle1,
                                        circle1_r * circle1_r * 4),
                               false);
  Point_2 p2_line_diagonal(3, 4);
  Line_arc_2 line_arc_diagonal(Line_2(center_circle1,
                                      p2_line_diagonal),
                               circle1,
                               true,
                               Circle_2(center_circle1,
                                        circle1_r * circle1_r * 4),
                               false);
  Circular_arc_2 arc_1(circle2,
                       Line_2(center_circle1,
                              p2_line_horizontal),
                       true,
                       Line_2(center_circle1,
                              p2_line_vertical),
                       false);
  Circular_arc_2 arc_2(circle1,
                       Line_2(center_circle1,
                              p2_line_horizontal),
                       true,
                       Line_2(center_circle1,
                              p2_line_vertical),
                       false);
   std::vector< Intersection_result >
     vector_for_intersection_1;
   theConstruct_intersect_2(line_arc_horizontal3,
                            arc_1,
                            std::back_inserter(vector_for_intersection_1));
   std::pair<Circular_arc_point_2, unsigned int> the_pair;
   assert(vector_for_intersection_1.size() == 2);
   assert(assign_variant(the_pair, vector_for_intersection_1[0]));
   Circular_arc_point_2 first = the_pair.first;
   assert(first == line_arc_horizontal3.source());
   assert(the_pair.second == 1);
   assert(assign_variant(the_pair, vector_for_intersection_1[1]));
   first = the_pair.first;
   assert(first == line_arc_horizontal3.target());
   assert(the_pair.second == 1);

 std::vector< Intersection_result >
     vector_for_intersection_2;
   theConstruct_intersect_2(line_arc_horizontal,
                            arc_1,
                            std::back_inserter(vector_for_intersection_2));
   assert(vector_for_intersection_2.size() == 1);
   assert(assign_variant(the_pair, vector_for_intersection_2[0]));
   first = the_pair.first;
   assert(first == line_arc_horizontal.target());
   assert(the_pair.second == 1);

   std::vector< Intersection_result >
     vector_for_intersection_3;
   theConstruct_intersect_2(line_arc_vertical,
                            arc_1,
                            std::back_inserter(vector_for_intersection_3));
   assert(vector_for_intersection_3.size() == 1);
   assert(assign_variant(the_pair, vector_for_intersection_3[0]));
   first = the_pair.first;
   assert(first == arc_1.target());
   assert(the_pair.second == 1);

   std::vector< Intersection_result >
     vector_for_intersection_4;
   theConstruct_intersect_2(Line_arc_2(Point_2(-10, -5), Point_2(5, 15)),
                            arc_2,
                            std::back_inserter(vector_for_intersection_4));
   assert(vector_for_intersection_4.size() == 0);

   std::vector< Intersection_result >
     vector_for_intersection_5;
   theConstruct_intersect_2(Line_arc_2(Point_2(10, -5), Point_2(-5, 15)),
                            arc_2,
                            std::back_inserter(vector_for_intersection_5));
   assert(vector_for_intersection_5.size() == 1);
   assert(assign_variant(the_pair, vector_for_intersection_5[0]));
   assert(the_pair.second == 2);

   //random
   int random_max = 127;
   int random_min = -127;
   Point_2 center_circle_random1(0,0);
   int circle_random1_r = theRandom.get_int(1, random_max);
   Point_2 p_random1;
   do{
     p_random1 = Point_2(theRandom.get_int(random_min, random_max) * circle_random1_r,
                         theRandom.get_int(random_min, random_max) * circle_random1_r);
   }while(p_random1 == center_circle_random1);
   Point_2 p_random2;
   do{
     p_random2 = Point_2(theRandom.get_int(random_min, random_max) * circle_random1_r,
                         theRandom.get_int(random_min, random_max) * circle_random1_r);
   }while(p_random2 == center_circle_random1);
   Circle_2 circle_random1(center_circle_random1, circle_random1_r * circle_random1_r);
   Circular_arc_2 arc_random_1(circle_random1,
                               Line_2(center_circle_random1,
                                      p_random1),
                               true,
                               Line_2(center_circle_random1,
                                      p_random2),
                               false);

   Circular_arc_point_2 first2;
   std::vector< Intersection_result >
     vector_for_intersection_random1;
   theConstruct_intersect_2(Line_arc_2(Point_2(-p_random1.x(),-p_random1.y()),
                                       p_random1),
                            arc_random_1,
                            std::back_inserter(vector_for_intersection_random1));
   assert(vector_for_intersection_random1.size() > 0);
   assert(assign_variant(the_pair, vector_for_intersection_random1[0]));
   first = the_pair.first;
   assert(the_pair.second == 1);
   assert(first == arc_random_1.source());



   std::vector< Intersection_result >
     vector_for_intersection_random2;
   theConstruct_intersect_2(Line_arc_2(Point_2(-p_random2.x(),-p_random2.y()),p_random2),
                            arc_random_1,
                            std::back_inserter(vector_for_intersection_random2));
   assert(vector_for_intersection_random2.size() > 0);
   assert(assign_variant(the_pair, vector_for_intersection_random2[0]));
   first = the_pair.first;
   assert(the_pair.second == 1);
   if(vector_for_intersection_random2.size() == 2){
     assert(assign_variant(the_pair, vector_for_intersection_random2[1]));
     first2 = the_pair.first;
     assert(the_pair.second == 1);
     assert((first == arc_random_1.target()) || (first2 == arc_random_1.target()));
   }
   else{
     assert(first == arc_random_1.target());
   }


   for (int loop = 0; loop < 200; loop++){
     Make_x_monotone_2 theMake_x_monotone = ck.make_x_monotone_2_object();
     Point_2 p_random3;
     do{
       p_random3 = Point_2(theRandom.get_int(random_min, random_max),
                           theRandom.get_int(random_min, random_max));
     }while(p_random3 == center_circle_random1);
     Point_2 p_random4;
     do{
       p_random4 = Point_2(theRandom.get_int(random_min, random_max),
                           theRandom.get_int(random_min, random_max));
     } while (p_random4 == center_circle_random1 || (p_random3 == p_random4));

     std::vector< Intersection_result >
       vector_for_intersection_random3;
     theConstruct_intersect_2(Line_arc_2(p_random3,p_random4),
                              arc_random_1,
                              std::back_inserter(vector_for_intersection_random3));
     for( std::size_t i = 0; i < vector_for_intersection_random3.size(); i++){
       assert(assign_variant(the_pair, vector_for_intersection_random3[i]));
       first = the_pair.first;
       std::vector<Intersection_result> objects_x_monotone;
       theMake_x_monotone( arc_random_1, std::back_inserter(objects_x_monotone));
       bool is_on_arc = false;
       for(std::size_t j = 0; j < objects_x_monotone.size(); j++){
         Circular_arc_2 aux;
         assign_variant(aux, objects_x_monotone[j]);
         if(CGAL::CircularFunctors::has_on<CK>(aux, first)){
           is_on_arc = true;
           break;
         }
       }
       assert(is_on_arc);
     }
   }
}

template <class CK>
void _test_compare_y_to_right(CK ck)
{

  typedef CGAL::Circle_2<CK>                   Circle_2;
  typedef CGAL::Circular_arc_2<CK>             Circular_arc_2;
  typedef CGAL::Point_2<CK>                    Point_2;
  typedef CGAL::Line_2<CK>                     Line_2;
  typedef CGAL::Line_arc_2<CK>                 Line_arc_2;
  typedef CGAL::Circular_arc_point_2<CK>       Circular_arc_point_2;
  typedef typename CK::Intersect_2   Intersect_2;
  typedef std::variant<std::pair<Circular_arc_point_2, unsigned int>, Line_arc_2> Intersection_result;

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);

  Intersect_2 theConstruct_intersect_2 = ck.intersect_2_object();
  Point_2 center_circle1(0, 0);
  int circle1_r = 5;
  Circle_2 circle1(center_circle1, circle1_r * circle1_r);
  Circle_2 circle2(center_circle1, circle1_r * circle1_r * 4);
  Point_2 p2_line_horizontal(1, 0);
  Point_2 p2_line_diagonal(1, 1);
  Line_arc_2 line_arc_horizontal(Line_2(center_circle1,
                                        p2_line_horizontal),
                                 circle1,
                                 true,
                                 circle2,
                                 false);
  Line_arc_2 line_arc_diagonal(Line_2(center_circle1,
                                      p2_line_diagonal),
                               circle1,true,
                               circle2,false);
  std::vector< Intersection_result >
     vector_for_intersection_1;
   theConstruct_intersect_2(line_arc_horizontal,
                            line_arc_diagonal,
                            std::back_inserter(vector_for_intersection_1));
   std::pair<Circular_arc_point_2, unsigned int> the_pair;
   assert(assign_variant(the_pair, vector_for_intersection_1[0]));
   Circular_arc_point_2 first = the_pair.first;

   assert(CGAL::CircularFunctors::compare_y_to_right<CK>(line_arc_horizontal,
                                               line_arc_diagonal,
                                                         first) == CGAL::SMALLER);

   assert(CGAL::CircularFunctors::compare_y_to_right<CK>(line_arc_diagonal,
                                                         line_arc_horizontal,
                                                         first) == CGAL::LARGER);
   assert(CGAL::CircularFunctors::compare_y_to_right<CK>(line_arc_diagonal,
                                                         line_arc_diagonal,
                                                         first) == CGAL::EQUAL);


  Circular_arc_2 part_high(circle1,
                           line_arc_horizontal.supporting_line(), false,
                           line_arc_horizontal.supporting_line(), true);

  Circular_arc_2 part_low(circle1,
                           line_arc_horizontal.supporting_line(), true,
                           line_arc_horizontal.supporting_line(), false);

  assert(CGAL::CircularFunctors::compare_y_to_right<CK>(line_arc_horizontal,
                                             part_high,
                                             part_high.target()) == CGAL::SMALLER);
  assert(CGAL::CircularFunctors::compare_y_to_right<CK>(line_arc_horizontal,
                                             part_low,
                                             part_low.source()) == CGAL::LARGER);
  Point_2 p2_high(0, circle1_r);
  Point_2 p2_high_right(circle1_r, circle1_r);
  Line_arc_2 line_arc_horizontal_high(p2_high,
                                      p2_high_right);
  assert(CGAL::CircularFunctors::compare_y_to_right<CK>(line_arc_horizontal_high,
                                             part_high,
                                             line_arc_horizontal_high.source()) == CGAL::LARGER);
  Point_2 p2_low(0, -circle1_r);
  Point_2 p2_low_right(circle1_r, -circle1_r);
  Line_arc_2 line_arc_horizontal_low(p2_low,
                                     p2_low_right
                                     );
  assert(CGAL::CircularFunctors::compare_y_to_right<CK>(line_arc_horizontal_low,
                                             part_low,
                                             line_arc_horizontal_low.source()) == CGAL::SMALLER);

   std::vector< Intersection_result >
     vector_for_intersection_2;
   theConstruct_intersect_2(part_high,
                            line_arc_diagonal,
                            std::back_inserter(vector_for_intersection_2));
   assert(assign_variant(the_pair, vector_for_intersection_2[0]));
   first = the_pair.first;
   assert(CGAL::CircularFunctors::compare_y_to_right<CK>(line_arc_diagonal,
                                                         part_high,
                                                         first) == CGAL::LARGER);
   std::vector< Intersection_result >
     vector_for_intersection_3;
   theConstruct_intersect_2(part_low,
                            line_arc_diagonal,
                            std::back_inserter(vector_for_intersection_3));
   assert(assign_variant(the_pair, vector_for_intersection_3[0]));
   first = the_pair.first;
   assert(CGAL::CircularFunctors::compare_y_to_right<CK>(line_arc_diagonal,
                                                         part_low,
                                                         first) == CGAL::LARGER);
   Point_2 p2_line_diagonal2(1, -1);
   Line_arc_2 line_arc_diagonal2(Line_2(center_circle1,
                                       p2_line_diagonal2),
                                circle1,true,
                                circle2,false);
   std::vector< Intersection_result >
     vector_for_intersection_4;
   theConstruct_intersect_2(part_high,
                            line_arc_diagonal2,
                            std::back_inserter(vector_for_intersection_4));
   assert(assign_variant(the_pair, vector_for_intersection_4[0]));
   first = the_pair.first;
   assert(CGAL::CircularFunctors::compare_y_to_right<CK>(line_arc_diagonal2,
                                                         part_high,
                                                         first) == CGAL::SMALLER);
   std::vector< Intersection_result >
     vector_for_intersection_5;
   theConstruct_intersect_2(part_low,
                            line_arc_diagonal2,
                            std::back_inserter(vector_for_intersection_5));
   assert(assign_variant(the_pair, vector_for_intersection_5[0]));
   first = the_pair.first;
   assert(CGAL::CircularFunctors::compare_y_to_right<CK>(line_arc_diagonal2,
                                                         part_low,
                                                         first) == CGAL::SMALLER);
}

template <class CK>
void _test_compare_y_at_x(CK ck)
{



  typedef CGAL::Point_2<CK>                    Point_2;

  typedef CGAL::Line_arc_2<CK>                 Line_arc_2;



  typedef typename CK::Compare_y_at_x_2        Compare_y_at_x_2;

  Line_arc_2 line_arc_horizontal(Point_2(-1, 0), Point_2(1, 0));
  Line_arc_2 line_arc_vertical(Point_2(0, 0), Point_2(0, 1));
  Line_arc_2 line_arc_vertical2(Point_2(0, -1), Point_2(0, 1));
  Line_arc_2 line_arc_diagonal(Point_2(-1, -1), Point_2(1, 1));
  Line_arc_2 line_arc_diagonal2(Point_2(-1, 1), Point_2(1, -1));
  Line_arc_2 line_arc_vertical3(Point_2(0, -1), Point_2(0, 2));
  Compare_y_at_x_2 theCompare_y_at_x_2 = ck.compare_y_at_x_2_object();
  assert(theCompare_y_at_x_2(line_arc_vertical.target(),
                             line_arc_horizontal) == CGAL::LARGER);
  assert(theCompare_y_at_x_2(line_arc_vertical2.source(),
                             line_arc_horizontal) == CGAL::SMALLER);
  assert(theCompare_y_at_x_2(line_arc_vertical.source(),
                             line_arc_horizontal) == CGAL::EQUAL);
  assert(theCompare_y_at_x_2(line_arc_vertical2.source(),
                             line_arc_vertical) == CGAL::SMALLER);
  assert(theCompare_y_at_x_2(line_arc_vertical2.target(),
                             line_arc_vertical) == CGAL::EQUAL);
  assert(theCompare_y_at_x_2(line_arc_vertical3.target(),
                             line_arc_vertical) == CGAL::LARGER);
  assert(theCompare_y_at_x_2(line_arc_vertical2.target(),
                             line_arc_diagonal) == CGAL::LARGER);
  assert(theCompare_y_at_x_2(line_arc_vertical2.source(),
                             line_arc_diagonal) == CGAL::SMALLER);
  assert(theCompare_y_at_x_2(line_arc_vertical2.target(),
                             line_arc_diagonal2) == CGAL::LARGER);
  assert(theCompare_y_at_x_2(line_arc_vertical2.source(),
                             line_arc_diagonal2) == CGAL::SMALLER);
  assert(theCompare_y_at_x_2(line_arc_vertical.source(),
                             line_arc_diagonal) == CGAL::EQUAL);
  assert(theCompare_y_at_x_2(line_arc_diagonal.source(),
                             line_arc_diagonal) == CGAL::EQUAL);
  assert(theCompare_y_at_x_2(line_arc_diagonal.target(),
                             line_arc_diagonal) == CGAL::EQUAL);

}



template <class CK>
void _test_has_on(CK)
{



  typedef CGAL::Point_2<CK>                    Point_2;

  typedef CGAL::Line_arc_2<CK>                 Line_arc_2;





  Line_arc_2 line_arc_horizontal(Point_2(-1, 0), Point_2(1, 0));
  Line_arc_2 line_arc_vertical(Point_2(0, 0), Point_2(0, 1));
  Line_arc_2 line_arc_vertical2(Point_2(0, -1), Point_2(0, 1));
  Line_arc_2 line_arc_diagonal(Point_2(-1, -1), Point_2(1, 1));
  Line_arc_2 line_arc_diagonal2(Point_2(-1, 1), Point_2(1, -1));
  Line_arc_2 line_arc_vertical3(Point_2(0, -1), Point_2(0, 2));
  assert(CGAL::CircularFunctors::has_on<CK>(line_arc_horizontal,
                                            line_arc_vertical.source()));
  assert(!CGAL::CircularFunctors::has_on<CK>(line_arc_horizontal,
                                             line_arc_vertical.target()));
  assert(!CGAL::CircularFunctors::has_on<CK>(line_arc_vertical,
                                             line_arc_vertical2.source()));
  assert(CGAL::CircularFunctors::has_on<CK>(line_arc_vertical,
                                             line_arc_vertical2.target()));
  assert(CGAL::CircularFunctors::has_on<CK>(line_arc_diagonal,
                                            line_arc_vertical.source()));

}

template <class K>
void do_test() {
  _test_Line_arc(K());
  _test_intersection_Line_arc_Circle(K());
  _test_intersection_Line_arc_Circular_arc(K());
  _test_compare_y_to_right(K());
  _test_compare_y_at_x(K());
  _test_has_on(K());
}

int main()
{
  typedef CGAL::Quotient<CGAL::MP_Float>                       NT1;
  typedef CGAL::Cartesian<NT1>                                 Linear_k1;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT1>          Algebraic_k1;
  typedef CGAL::Circular_kernel_2<Linear_k1,Algebraic_k1>      CK1;
        do_test< CK1 >();
        do_test< CGAL::Exact_circular_kernel_2 >();
}
