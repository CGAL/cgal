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

#include <CGAL/Random.h>
#include <cassert>

template <class CK>
void _test_circle_predicat(CK ck)
{
  typedef CGAL::Circle_2<CK>                   Circle_2;
  typedef CGAL::Circular_arc_2<CK>             Circular_arc_2;
  typedef CGAL::Point_2<CK>                    Point_2;
  typedef CGAL::Line_2<CK>                     Line_2;
  typedef CGAL::Line_arc_2<CK>                 Line_arc_2;
  typedef CGAL::Circular_arc_point_2<CK>       Circular_arc_point_2;
  typedef typename CK::Compare_x_2             Compare_x_2;
  typedef typename CK::Compare_y_2             Compare_y_2;
  typedef typename CK::Compare_xy_2            Compare_xy_2;
  typedef typename CK::Compare_y_to_right_2    Compare_y_to_right_2;
  typedef typename CK::Compare_y_at_x_2        Compare_y_at_x_2;
  typedef typename CK::Equal_2                 Equal_2;
  typedef typename CK::In_x_range_2            In_x_range_2;
  typedef typename CK::Do_overlap_2            Do_overlap_2;
  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127;

  //we make the circle1
  int center1_x = theRandom.get_int(random_min, random_max);
  int center1_y = theRandom.get_int(random_min, random_max);
  Point_2 center1(center1_x,center1_y);

  int circ1_r = theRandom.get_int(1, random_max);
  Circle_2 circ1(center1, circ1_r * circ1_r);

  //Comparison between endpoint

  Compare_x_2 theCompare_x_2 = ck.compare_x_2_object();
  Compare_y_2 theCompare_y_2 = ck.compare_y_2_object();
  Compare_xy_2 theCompare_xy_2 = ck.compare_xy_2_object();

  //we create a circle in low right of circle1
  //to have 2 different intersection's points
  //in quarter right inferior of the circle1
  Point_2 center1_low_right(center1_x + circ1_r, center1_y - circ1_r);
  Circle_2 circ1_low_right(center1_low_right, circ1_r * circ1_r);

  //p1 is lefter and lower than que p2
  Circular_arc_point_2 circ1_arc_end_p1 =
           CGAL::circle_intersect<CK>(circ1, circ1_low_right, true);
  Circular_arc_point_2 circ1_arc_end_p2 =
           CGAL::circle_intersect<CK>(circ1, circ1_low_right, false);

  assert(theCompare_x_2(circ1_arc_end_p1,circ1_arc_end_p2 )== CGAL::SMALLER);
  assert(theCompare_y_2(circ1_arc_end_p1,circ1_arc_end_p2) == CGAL::SMALLER);
  assert(theCompare_xy_2(circ1_arc_end_p1,circ1_arc_end_p2) == CGAL::SMALLER);

  assert(theCompare_x_2(circ1_arc_end_p2,circ1_arc_end_p2) == CGAL::EQUAL);
  assert(theCompare_y_2(circ1_arc_end_p2,circ1_arc_end_p2) == CGAL::EQUAL);
  assert(theCompare_xy_2(circ1_arc_end_p2,circ1_arc_end_p2) == CGAL::EQUAL);

  assert(theCompare_x_2(circ1_arc_end_p2,circ1_arc_end_p1) == CGAL::LARGER);
  assert(theCompare_y_2(circ1_arc_end_p2,circ1_arc_end_p1) == CGAL::LARGER);
  assert(theCompare_xy_2(circ1_arc_end_p2,circ1_arc_end_p1) == CGAL::LARGER);

  assert(compare_x(circ1_arc_end_p1,circ1_arc_end_p2 )== CGAL::SMALLER);
  assert(compare_y(circ1_arc_end_p1,circ1_arc_end_p2) == CGAL::SMALLER);
  assert(compare_xy(circ1_arc_end_p1,circ1_arc_end_p2) == CGAL::SMALLER);

  assert(compare_x(circ1_arc_end_p2,circ1_arc_end_p2) == CGAL::EQUAL);
  assert(compare_y(circ1_arc_end_p2,circ1_arc_end_p2) == CGAL::EQUAL);
  assert(compare_xy(circ1_arc_end_p2,circ1_arc_end_p2) == CGAL::EQUAL);

  assert(compare_x(circ1_arc_end_p2,circ1_arc_end_p1) == CGAL::LARGER);
  assert(compare_y(circ1_arc_end_p2,circ1_arc_end_p1) == CGAL::LARGER);
  assert(compare_xy(circ1_arc_end_p2,circ1_arc_end_p1) == CGAL::LARGER);

  //We create a circle in top of the circle1
  Point_2 center1_high(center1_x, center1_y + circ1_r);
  Circle_2 circ1_high(center1_high, circ1_r * circ1_r);

  //p3 is in the quarter superior left
  Circular_arc_point_2 circ1_arc_end_p3 =
           CGAL::circle_intersect<CK>(circ1, circ1_high, true);

  assert(theCompare_x_2(circ1_arc_end_p3, circ1_arc_end_p1 )== CGAL::SMALLER);
  assert(theCompare_y_2(circ1_arc_end_p3, circ1_arc_end_p1) == CGAL::LARGER);
  assert(theCompare_xy_2(circ1_arc_end_p3, circ1_arc_end_p1) == CGAL::SMALLER);

  assert(theCompare_x_2(circ1_arc_end_p1, circ1_arc_end_p3 )== CGAL::LARGER);
  assert(theCompare_y_2(circ1_arc_end_p1, circ1_arc_end_p3) == CGAL::SMALLER);
  assert(theCompare_xy_2(circ1_arc_end_p1, circ1_arc_end_p3) == CGAL::LARGER);


  // Comparison between an endpoint and an arc

  std::cout << "Compare_y_at_x_2_object" << std::endl;
  Compare_y_at_x_2 theCompare_y_at_x_2 = ck.compare_y_at_x_2_object();
  Point_2 point_2_left(center1_x - circ1_r, center1_y);
  Line_2 theLine_2_horizontal(center1, point_2_left);
  //The circ1_arc_high and circ1_arc_low are x_monotone
  Circular_arc_2 circ1_arc_low(circ1,
                               theLine_2_horizontal,true,
                               theLine_2_horizontal, false);
  Circular_arc_2 circ1_arc_high(circ1,
                                theLine_2_horizontal,false,
                                theLine_2_horizontal, true);

  //Comparison between the superior arc and one of its points (p3)
  CGAL::Comparison_result theComparison_result_y_at_x_2 =
    theCompare_y_at_x_2(circ1_arc_end_p3,circ1_arc_high);
  assert(theComparison_result_y_at_x_2 == CGAL::EQUAL);
  assert(compare_y_at_x(circ1_arc_end_p3,circ1_arc_high) == CGAL::EQUAL);
  //Comparison between the inferior arc and a point in top (p3)
  theComparison_result_y_at_x_2 =
    theCompare_y_at_x_2(circ1_arc_end_p3,circ1_arc_low);
  assert(theComparison_result_y_at_x_2 == CGAL::LARGER);

  //Comparison between the superior arc and a point in bottom (p1)
  theComparison_result_y_at_x_2 =
    theCompare_y_at_x_2(circ1_arc_end_p1,circ1_arc_high);
  assert(theComparison_result_y_at_x_2 == CGAL::SMALLER);

  //We create a circle tangent at left of the cercle1
  Point_2 center1_left(center1_x - (2 * circ1_r), center1_y);
  Circle_2 circ1_left(center1_left, circ1_r * circ1_r);
  //P4 is the leftest point of the circle1
  Circular_arc_point_2 circ1_arc_end_p4 =
        CGAL::circle_intersect<CK>(circ1, circ1_left, true);

  //Comparison between the superior arc and the point p4
  theComparison_result_y_at_x_2 =
    theCompare_y_at_x_2(circ1_arc_end_p4,circ1_arc_high);
  assert(theComparison_result_y_at_x_2 == CGAL::EQUAL);

  //Comparison between the inferior arc and the point p4
  theComparison_result_y_at_x_2 =
    theCompare_y_at_x_2(circ1_arc_end_p4,circ1_arc_low);
  assert(theComparison_result_y_at_x_2 == CGAL::EQUAL);

  //Comparison between the superior arc and a point
  //on this arc but create with two others arcs
  Point_2 center1_1(center1_x - circ1_r, center1_y + circ1_r);
  Circle_2 circ1_1(center1_1, circ1_r * circ1_r);
  Point_2 center1_2(center1_x + circ1_r, center1_y + circ1_r);
  Circle_2 circ1_2(center1_2, circ1_r * circ1_r);
  Circular_arc_point_2 circ1_arc_end_1_1_1_2 =
         CGAL::circle_intersect<CK>(circ1_1, circ1_2, true);
  theComparison_result_y_at_x_2 =
    theCompare_y_at_x_2(circ1_arc_end_1_1_1_2,circ1_arc_high);
  assert(theComparison_result_y_at_x_2 == CGAL::EQUAL);

  //Comparison between an inferior arc and a point in top create
  //with 2 others circles (p1)
  Circle_2 circ1_big(center1, circ1_r * circ1_r * 4);
  Circular_arc_2 circ1_big_arc_low(circ1_big,
                               theLine_2_horizontal,true,
                               theLine_2_horizontal, false);
  theComparison_result_y_at_x_2 =
    theCompare_y_at_x_2(circ1_arc_end_p1,circ1_big_arc_low);
  assert(theComparison_result_y_at_x_2 == CGAL::LARGER);



  //Comparison between 2 arc at right of one of their intersection's point

  std::cout << "Compare_y_to_right_2_object" << std::endl;
  Point_2 point_2_left_circ_high(center1_x - circ1_r,
                                 center1_y + circ1_r);
  Line_2 theLine_2_horizontal_circ_high(center1_high,
                                        point_2_left_circ_high);
  //it's the inferior arc of the circle circ_high
  Circular_arc_2 circ_high_arc_low(circ1_high,
                                   theLine_2_horizontal_circ_high, true,
                                   theLine_2_horizontal_circ_high, false);
  Compare_y_to_right_2 theCompare_y_to_right_2 =
    ck.compare_y_to_right_2_object();
  CGAL::Comparison_result theComparison_result_y_to_right_2 =
    theCompare_y_to_right_2(circ1_arc_high,
                            circ_high_arc_low,
                            circ1_arc_end_p3);
        CGAL::Comparison_result theComparison_result_y_to_right_2_l =
          compare_y_to_right(circ1_arc_high,
                                                    circ_high_arc_low,
                                                    circ1_arc_end_p3);
  assert(theComparison_result_y_to_right_2 == CGAL::LARGER);
  assert(theComparison_result_y_to_right_2_l == CGAL::LARGER);

  theComparison_result_y_to_right_2=
    theCompare_y_to_right_2(circ_high_arc_low,
                            circ1_arc_high,
                            circ1_arc_end_p3);
  assert(theComparison_result_y_to_right_2 == CGAL::SMALLER);

  theComparison_result_y_to_right_2 =
    theCompare_y_to_right_2(circ_high_arc_low,
                            circ_high_arc_low,
                            circ1_arc_end_p3);
  assert(theComparison_result_y_to_right_2 == CGAL::EQUAL);

  Circular_arc_point_2 circ1_big_arc_end =
        CGAL::circle_intersect<CK>(circ1_big, circ1_high, true);
  Circular_arc_2 circ1_big_arc_high(circ1_big,
                                   theLine_2_horizontal,false,
                                   theLine_2_horizontal, true);
  Circular_arc_2 circ_high_arc_high(circ1_high,
                                   theLine_2_horizontal_circ_high, false,
                                   theLine_2_horizontal_circ_high, true);
  theComparison_result_y_to_right_2 =
    theCompare_y_to_right_2(circ1_big_arc_high,
                            circ_high_arc_high,
                            circ1_big_arc_end);
  assert(theComparison_result_y_to_right_2 == CGAL::LARGER);

  Point_2 center1_very_high(center1_x, center1_y + 2*circ1_r);
  Circle_2 circ1_very_high(center1_very_high, circ1_r * circ1_r);
  Point_2 point_left_circ1_very_high(center1_x - circ1_r,
                                     center1_y + 2*circ1_r);
  Line_2 theLine_2_horizontal_circ_very_high(center1_very_high,
                                             point_left_circ1_very_high);
  Circular_arc_2 circ_very_high_arc_low (circ1_very_high,
                                   theLine_2_horizontal_circ_very_high, true,
                                   theLine_2_horizontal_circ_very_high, false);
  Circular_arc_point_2 circ1_arc_end_tangent =
          CGAL::circle_intersect<CK>(circ1, circ1_very_high, true);
  theComparison_result_y_to_right_2 =
    theCompare_y_to_right_2(circ1_arc_high,
                            circ_very_high_arc_low,
                            circ1_arc_end_tangent);
  assert(theComparison_result_y_to_right_2 == CGAL::SMALLER);




  //P3 is the leftest point of circle1
  //P3, P4 and P5 are the same
  Circular_arc_point_2 circ1_arc_end_p5 =
         CGAL::circle_intersect<CK>(circ1, circ1_left, true);
  Circular_arc_point_2 circ1_arc_end_p6 =
         CGAL::circle_intersect<CK>(circ1, circ1_left, false);
  Circular_arc_point_2 circ1_arc_end_p7 =
         CGAL::circle_intersect<CK>(circ1_left, circ1, true);

  ///////////////////////EQUAL///////////////////////////////

  std::cout << "Equal" << std::endl;
  Equal_2 theEqual_2 = ck.equal_2_object();
  assert(theEqual_2(circ1_arc_end_p1,circ1_arc_end_p1));
  assert(!theEqual_2(circ1_arc_end_p1,circ1_arc_end_p2));
  assert(theEqual_2(circ1_arc_end_p5,circ1_arc_end_p6));
  assert(theEqual_2(circ1_arc_end_p5,circ1_arc_end_p7));

  assert(theEqual_2(circ1_arc_low, circ1_arc_low));
  assert(!theEqual_2(circ1_arc_low, circ1_arc_high));
  std::cout << std::endl;


  //////////////////////////IN RANGE//////////////////////////

  //We create a circle with center of c1 and a radius = circ1_r / 2.0
  Point_2 center2(center1_x, center1_y);
  double circ2_r = circ1_r / 2.0;
  Circle_2 circ2(center2, circ2_r * circ2_r);
  //we create a same circle than circ2 but offset towards the left
  Point_2 center2_left(center1_x - circ2_r, center1_y);
  Circle_2 circ2_left(center2_left, circ2_r * circ2_r);
  //Point on top
  Circular_arc_point_2 circ2_arc_end_p1 =
          CGAL::circle_intersect<CK>(circ2, circ2_left, false);
  //We create a circle lefter than circ1_left but it cuts it
  Point_2 center1_very_left(center1_x - (3 * circ1_r), center1_y);
  Circle_2 circ1_very_left(center1_very_left, circ1_r * circ1_r);
  //Point on top
  Circular_arc_point_2 circ1_left_arc_end_p1 =
         CGAL::circle_intersect<CK>(circ1_left, circ1_very_left, false);
  std::cout << "In x range" << std::endl;
  In_x_range_2 theIn_x_range_2 = ck.in_x_range_2_object();
  assert(theIn_x_range_2(circ1_arc_low, circ1_arc_end_p7));
  assert(theIn_x_range_2(circ1_arc_high, circ1_arc_end_p7));
  assert(theIn_x_range_2(circ1_arc_high, circ2_arc_end_p1));
  assert(theIn_x_range_2(circ1_arc_low, circ2_arc_end_p1));
  assert(!theIn_x_range_2(circ1_arc_low, circ1_left_arc_end_p1));

  assert(has_in_x_range(circ1_arc_low, circ1_arc_end_p7));
  assert(has_in_x_range(circ1_arc_high, circ1_arc_end_p7));
  assert(has_in_x_range(circ1_arc_high, circ2_arc_end_p1));
  assert(has_in_x_range(circ1_arc_low, circ2_arc_end_p1));
  assert(!has_in_x_range(circ1_arc_low, circ1_left_arc_end_p1));

  std::cout << std::endl;


  //////////////////////////DO OVERLAP/////////////////////////

  std::cout << "Do_overlap_2_object" << std::endl;
  Do_overlap_2 theDo_overlap_2 = ck.do_overlap_2_object();
  Point_2 point_2_high(center1_x, center1_y + circ1_r);
  Line_2 theLine_2_vertical(center1,point_2_high);
  Circular_arc_2 circ1_arc_low_right(circ1,
                                     theLine_2_vertical,true,
                                     theLine_2_horizontal, false);

  //The following commented test is wrong, cause 2 arcs with the
  //same supporting circle can overlaps on one point
  //assert(!theDo_overlap_2(circ1_arc_high, circ1_arc_low_right));
  assert(theDo_overlap_2(circ1_arc_high, circ1_arc_low_right));
  assert(theDo_overlap_2(circ1_arc_low, circ1_arc_low_right));
  assert(theDo_overlap_2(circ1_arc_low_right, circ1_arc_low));
  assert(theDo_overlap_2(circ1_arc_low, circ1_arc_high));
  //The following commented test is wrong
  //1 circle and 1 half-circle overlap
  //assert(!theDo_overlap_2(circ1_arc_low, circ1_arc_high));
  std::cout << std::endl;


  //////////////////////IS_Y_MONOTONE///////////////////////
  std::cout << "is_y_monotone" << std::endl;
  Line_2 theLine_low_right(center1, Point_2(center1_x + circ1_r,
                                            center1_y - circ1_r));
  Line_2 theLine_low_left(center1, Point_2(center1_x - circ1_r,
                                           center1_y - circ1_r));
  Circular_arc_2 arc_y_monotone_1(circ1,
                                  theLine_2_vertical, true,
                                  theLine_2_horizontal, false);
  Circular_arc_2 arc_y_monotone_2(circ1,
                                  theLine_2_vertical, true,
                                  theLine_2_vertical, false);
  Circular_arc_2 arc_y_monotone_3(circ1,
                                  theLine_2_vertical, false,
                                  theLine_2_horizontal, true);
  Circular_arc_2 arc_y_monotone_4(circ1,
                                  theLine_low_right, true,
                                  theLine_low_left, true);
  Circular_arc_2 arc_no_y_monotone_1(circ1,
                                     theLine_2_horizontal, true,
                                     theLine_2_horizontal, false);
  Circular_arc_2 arc_no_y_monotone_2(circ1,
                                     theLine_2_horizontal, false,
                                     theLine_2_horizontal, true);
  Circular_arc_2 arc_no_y_monotone_3(circ1,
                                     theLine_2_horizontal, true,
                                     theLine_2_horizontal, true);
  assert(arc_y_monotone_1.is_y_monotone());
  assert(arc_y_monotone_2.is_y_monotone());
  assert(arc_y_monotone_3.is_y_monotone());
  assert(arc_y_monotone_4.is_y_monotone());
  assert(!arc_no_y_monotone_1.is_y_monotone());
  assert(!arc_no_y_monotone_2.is_y_monotone());
  assert(!arc_no_y_monotone_3.is_y_monotone());

  //////////////////////Bbox of an endpoint////////////////////////////
  for(int i = 0; i < 100; i++){
  Point_2 center_random(theRandom.get_int(random_min, random_max),
                        theRandom.get_int(random_min, random_max));
  Circle_2 circle_random(center_random,
                         theRandom.get_int(1, random_max));
  int x_random1, y_random1;
  int x_random2, y_random2;
  do{
    x_random1 = theRandom.get_int(random_min, random_max);
    y_random1 = theRandom.get_int(random_min, random_max);
  }while(x_random1 == 0 && y_random1 ==0);

  do{
    x_random2 = theRandom.get_int(random_min, random_max);
    y_random2 = theRandom.get_int(random_min, random_max);
  }while(x_random2 == 0 && y_random2 ==0);

  Line_2 line_random_1(center_random,
                       Point_2(center_random.x() +
                               x_random1,
                               center_random.y() +
                               y_random1));
  Line_2 line_random_2(center_random,
                       Point_2(center_random.x() +
                               x_random2,
                               center_random.y() +
                               y_random2));
  Circular_arc_2 arc_random(circle_random,
                            line_random_1, theRandom.get_bool(),
                            line_random_2, theRandom.get_bool());
  CGAL::Bbox_2 bb = arc_random.source().bbox();
  assert(typename CK::FT(bb.xmin()) <= arc_random.source().x());
  assert(typename CK::FT(bb.xmax()) >= arc_random.source().x());
  assert(typename CK::FT(bb.ymin()) <= arc_random.source().y());
  assert(typename CK::FT(bb.ymax()) >= arc_random.source().y());
  bb = arc_random.target().bbox();
  assert(typename CK::FT(bb.xmin()) <= arc_random.target().x());
  assert(typename CK::FT(bb.xmax()) >= arc_random.target().x());
  assert(typename CK::FT(bb.ymin()) <= arc_random.target().y());
  assert(typename CK::FT(bb.ymax()) >= arc_random.target().y());
  }

  // Testing Comparison Operators
  Circular_arc_point_2 p[3];
        p[0] = Circular_arc_point_2(Point_2(1,0));
  p[1] = Circular_arc_point_2(Point_2(1,0));
  p[2] = Circular_arc_point_2(Point_2(0,1));
  std::cout << "Testing lexico_operations(Circular_arc_point, Circular_arc_point)..." << std::endl;
        assert(p[0] > p[2]);
        assert(p[0] >= p[1]);
        assert(p[0] <= p[1]);
        assert(p[2] < p[0]);

        // TEST THE FUNCTOR CALL (VC8 porting mainly reason)
        Circular_arc_2 ccaa =
          typename CK::Construct_circular_arc_2()(Point_2(1, 2), Point_2(2, 2), Point_2(3, 3));
        Line_arc_2 llaa =
          typename CK::Construct_line_arc_2()(Point_2(2, 1), Point_2(2, 2));
        Circular_arc_point_2 ccaapp = typename CK::Construct_circular_arc_point_2()(Point_2(1, 2));

        assert(typename CK::Has_on_2()(ccaa, ccaapp));
        assert(typename CK::Is_vertical_2()(llaa));
        assert(typename CK::Compute_circular_x_2()(ccaapp) == 1);
        assert(typename CK::Compute_circular_y_2()(ccaapp) == 2);
        assert(typename CK::Is_x_monotone_2()(ccaa));
        assert(!(typename CK::Is_y_monotone_2()(ccaa)));
        assert(typename CK::Bounded_side_2()(ccaa.supporting_circle(), ccaapp) != CGAL::ON_BOUNDED_SIDE);
        assert(typename CK::Bounded_side_2()(ccaa.supporting_circle(), ccaapp) != CGAL::ON_UNBOUNDED_SIDE);
        assert(!(typename CK::Has_on_bounded_side_2()(ccaa.supporting_circle(), ccaapp)));
        assert(!(typename CK::Has_on_unbounded_side_2()(ccaa.supporting_circle(), ccaapp)));

}
