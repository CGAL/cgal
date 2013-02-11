// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#include <CGAL/Random.h>
#include <cassert>
#include <CGAL/use.h>

//#define typename 
template <class CK>
void _test_circle_construct(CK ck)
{
  typedef CGAL::Circle_2<CK>                       Circle_2;
  typedef CGAL::Circular_arc_2<CK>                 Circular_arc_2;
  typedef CGAL::Point_2<CK>                        Point_2;
  typedef CGAL::Line_2<CK>                         Line_2;
  typedef CGAL::Line_arc_2<CK>                     Line_arc_2;
  typedef CGAL::Circular_arc_point_2<CK>           Circular_arc_point_2;
  typedef typename CK::RT                          RT;
  typedef typename CK::FT                          FT;
  typedef typename CK::Construct_circle_2          Construct_circle_2;
  typedef typename CK::Intersect_2   Intersect_2;
  typedef typename CK::Make_x_monotone_2           Make_x_monotone_2;
  typedef typename CK::Make_xy_monotone_2           Make_xy_monotone_2;
  typedef typename CK::Split_2                     Split_2;  
  typedef typename CK::Get_equation                Get_equation;
  typedef typename CK::Compare_xy_2                Compare_xy_2;
  typedef typename CK::Do_intersect_2              Do_intersect_2;

  //fix warnings with g++-4.8 [-Wunused-local-typedefs]
  CGAL_USE_TYPE(Construct_circle_2);
  CGAL_USE_TYPE(Get_equation);

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127;
  
  FT sqrt2 = std::sqrt(2.0)/2;

  //test of get_equation_object()
  int x_equation = theRandom.get_int(random_min,random_max);
  int y_equation = theRandom.get_int(random_min,random_max);
  int r_equation = theRandom.get_int(1,random_max);
  Point_2 center_circ_equation(x_equation,y_equation);
  Circle_2 circ_equation(center_circ_equation, r_equation);
  std::cout << "the circle used by the equation :" 
	    << circ_equation << std::endl;
	
	std::cout << "testing {x,y}_extremal_points" << std::endl;
	for(int i=0; i<20; i++) {
	  int x1 = theRandom.get_int(random_min,random_max);
	  int y1 = theRandom.get_int(random_min,random_max);
	  int x2 = theRandom.get_int(random_min,random_max);
	  int y2 = theRandom.get_int(random_min,random_max);
	  int x3 = theRandom.get_int(random_min,random_max);
	  int y3 = theRandom.get_int(random_min,random_max);
		if(x1 == x2 && y1 == y2) continue;
		if(x1 == x3 && y1 == y3) continue;
		if(x2 == x3 && y2 == y3) continue;
		if(CGAL::collinear(Point_2(x1,y1), Point_2(x2,y2), Point_2(x3,y3))) continue;
	  Circular_arc_2 ca(Point_2(x1,y1), Point_2(x2,y2), Point_2(x3,y3));
		Circle_2 c = ca.supporting_circle();
		Circular_arc_point_2 cp_x_min = x_extremal_point(c, true);
		Circular_arc_point_2 cp_x_max = x_extremal_point(c, false);
		Circular_arc_point_2 cp_y_min = y_extremal_point(c, true);
		Circular_arc_point_2 cp_y_max = y_extremal_point(c, false);
		assert(CGAL_NTS square(cp_x_min.x() - c.center().x()) == c.squared_radius());
		assert(CGAL_NTS square(cp_x_max.x() - c.center().x()) == c.squared_radius());
		assert(CGAL_NTS square(cp_y_min.y() - c.center().y()) == c.squared_radius());
		assert(CGAL_NTS square(cp_y_max.y() - c.center().y()) == c.squared_radius());
		assert(cp_x_min.x() < cp_x_max.x());
		assert(cp_y_min.y() < cp_y_max.y());
  }

  //Constuct_intersections_2 with 2 intersection's points
  std::cout << std::endl << "construct_intersection_2" << std::endl;
  Do_intersect_2 theDo_intersect_2 = ck.do_intersect_2_object();
  Intersect_2 theConstruct_intersect_2 
    = ck.intersect_2_object();
  int center_circ_intersection_2_1_x = theRandom.get_int(random_min, random_max);
  int center_circ_intersection_2_1_y = theRandom.get_int(random_min, random_max);
  int circ_intersection_2_1_r = theRandom.get_int(1, random_max);
  Point_2 center_circ_intersections_2_1(center_circ_intersection_2_1_x, 
					center_circ_intersection_2_1_y);
  Circle_2 circ_intersections_2_1(center_circ_intersections_2_1, 
				  circ_intersection_2_1_r 
				  * circ_intersection_2_1_r);
  Point_2 center_circ_intersections_2_2(center_circ_intersection_2_1_x 
					+ circ_intersection_2_1_r,
					center_circ_intersection_2_1_y);
  Circle_2 circ_intersections_2_2(center_circ_intersections_2_2,
				  circ_intersection_2_1_r);
   
  std::vector< CGAL::Object > 
    vector_for_intersection_1, vector_for_intersection_1l;
  
  theConstruct_intersect_2(circ_intersections_2_1, 
			   circ_intersections_2_2,
			   std::back_inserter(vector_for_intersection_1));
	intersection(circ_intersections_2_1, 
					   circ_intersections_2_2,
					   std::back_inserter(vector_for_intersection_1l));
  // there are 2 intersection's points
	assert(theDo_intersect_2(circ_intersections_2_1, circ_intersections_2_1));
	assert(do_intersect(circ_intersections_2_1, circ_intersections_2_1));
  std::pair<Circular_arc_point_2, unsigned > the_pair;
  assert(assign(the_pair, vector_for_intersection_1[0]));
  assert(assign(the_pair, vector_for_intersection_1l[0]));
  Circular_arc_point_2 first = the_pair.first;
  std::cout << first << std::endl;
  assert(assign(the_pair, vector_for_intersection_1[1]));
  assert(assign(the_pair, vector_for_intersection_1l[1]));
  Circular_arc_point_2 second = the_pair.first;
  std::cout << second << std::endl;
  Compare_xy_2 theCompare_xy_2 = ck.compare_xy_2_object();
  assert(theCompare_xy_2(first, second) == CGAL::SMALLER);

  //Constuct_intersections_2 with 1 intersection's point
  Point_2 center_circ_intersections_2_3(center_circ_intersection_2_1_x 
					+ 2 * circ_intersection_2_1_r,
					center_circ_intersection_2_1_y);
  Circle_2 circ_intersections_2_3(center_circ_intersections_2_3,
				  circ_intersection_2_1_r 
				  * circ_intersection_2_1_r);
  Circular_arc_point_2 the_intersection_point_1 =
          CGAL::circle_intersect<CK>(circ_intersections_2_1,
				     circ_intersections_2_3,
				     true);
  Circular_arc_point_2 the_intersection_point_2 =
          CGAL::circle_intersect<CK>(circ_intersections_2_1,
				     circ_intersections_2_3,
				     false); 
  std::vector< CGAL::Object > 
    vector_for_intersection_2;
  theConstruct_intersect_2(circ_intersections_2_1, 
			   circ_intersections_2_3,
			   std::back_inserter(vector_for_intersection_2));
	assert(theDo_intersect_2(circ_intersections_2_1, circ_intersections_2_3));
  assert(vector_for_intersection_2.size() == 1);
  assign(the_pair, vector_for_intersection_2[0]);
  assert(the_pair.first == the_intersection_point_1);
  assert(the_pair.first == the_intersection_point_2);
  
  Point_2 center_circ_intersections_2_3_bis(center_circ_intersection_2_1_x 
					    + circ_intersection_2_1_r,
					    center_circ_intersection_2_1_y);
  Circle_2 circ_intersections_2_3_bis(center_circ_intersections_2_3_bis,
				      circ_intersection_2_1_r 
				      * circ_intersection_2_1_r * 4);
  std::vector< CGAL::Object > 
    vector_for_intersection_2_bis;
  theConstruct_intersect_2(circ_intersections_2_1, 
			   circ_intersections_2_3_bis,
			   std::back_inserter(vector_for_intersection_2_bis));
	assert(theDo_intersect_2(circ_intersections_2_1, circ_intersections_2_3_bis));
  assert(vector_for_intersection_2_bis.size() == 1);
  assign(the_pair, vector_for_intersection_2_bis[0]);
  assert(the_pair.second == 2u);
  

  //With circular arc
  Point_2 center_circ_intersections_2_1_low(center_circ_intersection_2_1_x, 
					center_circ_intersection_2_1_y - 
					    circ_intersection_2_1_r);
  Circle_2 circ_intersections_2_1_low(center_circ_intersections_2_1_low, 
				      circ_intersection_2_1_r 
				      * circ_intersection_2_1_r);
  Line_2 line_horizontal_circ_2_1_low(center_circ_intersections_2_1_low,
				      Point_2(center_circ_intersection_2_1_x - 
					      circ_intersection_2_1_r,
					      center_circ_intersection_2_1_y - 
					      circ_intersection_2_1_r));
  Circular_arc_2 circ_arc_2_1_low_part_high(circ_intersections_2_1_low,
				       line_horizontal_circ_2_1_low, false,
				       line_horizontal_circ_2_1_low, true);

  Line_2 line_horizontal(center_circ_intersections_2_1,
			 Point_2(center_circ_intersection_2_1_x - 
				 circ_intersection_2_1_r,
				 center_circ_intersection_2_1_y));
  Circular_arc_2 circ_arc_2_1_part_low(circ_intersections_2_1,
				       line_horizontal, true,
				       line_horizontal, false);
  //////////////if(circ_arc_2_1_low_part_high.center() == 
  //////////////   circ_arc_2_1_part_low.center()) {
  //////////////  std::cout << "OH NO!" << std::endl;    	
  //////////////} else std::cout << "OK" << std::endl;    	
  
  std::vector< CGAL::Object > 
    vector_for_intersection_3;
  theConstruct_intersect_2(circ_arc_2_1_part_low, 
			   circ_arc_2_1_low_part_high,
			   std::back_inserter(vector_for_intersection_3));
	assert(theDo_intersect_2(circ_arc_2_1_part_low, circ_arc_2_1_low_part_high));
	
  /////////////std::cout << "The size: " << vector_for_intersection_3.size() << std::endl;			   
  assert(vector_for_intersection_3.size() == 2);
  assign(the_pair, vector_for_intersection_3[0]);
  assert(the_pair.second == 1u);
  assert(the_pair.first == CGAL::circle_intersect<CK>(circ_intersections_2_1, 
						      circ_intersections_2_1_low, 
						      true));
  assign(the_pair, vector_for_intersection_3[1]);
  assert(the_pair.second == 1u);
  assert(the_pair.first == CGAL::circle_intersect<CK>(circ_intersections_2_1, 
						      circ_intersections_2_1_low, 
						      false));
  
    


  
  std::cout << "intersection with overlap arc_circle" << std::endl;

  Point_2 point_arc_overlap_left(center_circ_intersection_2_1_x
				 - circ_intersection_2_1_r, 
				 center_circ_intersection_2_1_y);
  Point_2 point_arc_overlap_low_right(center_circ_intersection_2_1_x 
				      + circ_intersection_2_1_r , 
				      center_circ_intersection_2_1_y
				      - circ_intersection_2_1_r);
  Point_2 point_arc_overlap_low_left(center_circ_intersection_2_1_x 
				      - circ_intersection_2_1_r , 
				      center_circ_intersection_2_1_y
				      - circ_intersection_2_1_r);
  Line_2 line_arc_overlap_horizontal(center_circ_intersections_2_1,
				   point_arc_overlap_left);
  Line_2 line_arc_overlap_low_right(center_circ_intersections_2_1,
				   point_arc_overlap_low_right);
  Line_2 line_arc_overlap_low_left(center_circ_intersections_2_1,
				    point_arc_overlap_low_left);
  //circ_arc_overlap_1 and circ_arc_overlap_2 are overlap on a circular_arc
  Circular_arc_2 circ_arc_overlap_1(circ_intersections_2_1,
				    line_arc_overlap_horizontal, true,
				    line_arc_overlap_low_right, false);
  Circular_arc_2 circ_arc_overlap_2(circ_intersections_2_1,
				    line_arc_overlap_low_left, true,
				    line_arc_overlap_horizontal, false);
  //circ_arc_overlap_1 and circ_arc_overlap_3 are overlap in one point
  Circular_arc_2 circ_arc_overlap_3(circ_intersections_2_1,
				    line_arc_overlap_low_right, false,
				    line_arc_overlap_horizontal, false);
  Circular_arc_2 circ_arc_overlap_upper_part(circ_intersections_2_1,
				    line_arc_overlap_horizontal, false,
				    line_arc_overlap_horizontal, true);
  Circular_arc_2 circ_arc_overlap_lower_part(circ_intersections_2_1,
				    line_arc_overlap_horizontal, true,
				    line_arc_overlap_horizontal, false);
  assert(circ_arc_overlap_1.is_x_monotone() 
	 && circ_arc_overlap_2.is_x_monotone()
	 && circ_arc_overlap_3.is_x_monotone()
	 && circ_arc_overlap_upper_part.is_x_monotone()
	 && circ_arc_overlap_lower_part.is_x_monotone()
	 );

  std::cout << "Intersection : same circular arc" << std::endl;
  std::vector< CGAL::Object > 
    vector_for_intersection_the_same_arc;
  theConstruct_intersect_2(circ_arc_overlap_1, 
 			   circ_arc_overlap_1,
 			   std::back_inserter(vector_for_intersection_the_same_arc));
  assert(theDo_intersect_2(circ_arc_overlap_1, circ_arc_overlap_1));
  assert(vector_for_intersection_the_same_arc.size() == 1);
  Circular_arc_2 res_same;
  assert(assign(res_same, vector_for_intersection_the_same_arc[0]));
  assert(res_same.source() == circ_arc_overlap_1.source());
  assert(res_same.target() == circ_arc_overlap_1.target());

  std::cout << "Intersection : overlap on a circular arc" << std::endl;
  Circular_arc_2 circ_arc_in_overlap;
  Circular_arc_2 circ_arc_in_overlap_2;
  std::vector< CGAL::Object > 
    vector_for_intersection_overlap_1_1;
  theConstruct_intersect_2(circ_arc_overlap_2, 
 			   circ_arc_overlap_1,
 			   std::back_inserter(vector_for_intersection_overlap_1_1));
  assert(theDo_intersect_2(circ_arc_overlap_2, circ_arc_overlap_1));
  Circular_arc_2 circ_arc_overlap_result(circ_intersections_2_1,
				    line_arc_overlap_low_left, true,
				    line_arc_overlap_low_right, false);
  assert(vector_for_intersection_overlap_1_1.size() == 1);
  assign(circ_arc_in_overlap, vector_for_intersection_overlap_1_1[0]);
  assert(circ_arc_in_overlap.source() == circ_arc_overlap_result.source());
  assert(circ_arc_in_overlap.target() == circ_arc_overlap_result.target());
  std::vector< CGAL::Object > 
    vector_for_intersection_overlap_1_2;
  theConstruct_intersect_2(circ_arc_overlap_2, 
 			   circ_arc_overlap_1,
 			   std::back_inserter(vector_for_intersection_overlap_1_2));
  assert(theDo_intersect_2(circ_arc_overlap_2, circ_arc_overlap_1));
  assert(vector_for_intersection_overlap_1_2.size() == 1);
  assign(circ_arc_in_overlap, vector_for_intersection_overlap_1_2[0]);
  assert(circ_arc_in_overlap.source() == circ_arc_overlap_result.source());
  assert(circ_arc_in_overlap.target() == circ_arc_overlap_result.target());
  
  
  std::cout << "Intersection : overlap in one point" << std::endl;
  std::vector< CGAL::Object > 
    vector_for_intersection_overlap_2_1;
  theConstruct_intersect_2(circ_arc_overlap_1, 
			   circ_arc_overlap_3,
			   std::back_inserter(vector_for_intersection_overlap_2_1));
  assert(theDo_intersect_2(circ_arc_overlap_1, circ_arc_overlap_3));
  assert(vector_for_intersection_overlap_2_1.size() == 1);
  assign(the_pair, vector_for_intersection_overlap_2_1[0]);
  std::cout << "x = " << the_pair.first.x() << " the result must be = " <<
    center_circ_intersection_2_1_x + circ_intersection_2_1_r * sqrt2 
	    << std::endl;
  
  assert(the_pair.first.x() * 
	 (center_circ_intersection_2_1_x + circ_intersection_2_1_r * sqrt2) >= 0);

  assert(square(the_pair.first.x() - RT(center_circ_intersection_2_1_x))
	 == (circ_intersection_2_1_r * circ_intersection_2_1_r / typename CK::RT(2)));
  std::cout << "y = " << the_pair.first.y() << " the result must be = " <<
    center_circ_intersection_2_1_y - circ_intersection_2_1_r * sqrt2 
	    << std::endl;

  assert(the_pair.first.y() * 
	 (center_circ_intersection_2_1_y - circ_intersection_2_1_r * sqrt2) >= 0);
  assert(square(the_pair.first.y() - RT(center_circ_intersection_2_1_y))
	 == (circ_intersection_2_1_r * circ_intersection_2_1_r / typename CK::RT(2)));

  std::vector< CGAL::Object > 
    vector_for_intersection_overlap_2_2;
  theConstruct_intersect_2(circ_arc_overlap_3, 
			   circ_arc_overlap_1,
			   std::back_inserter(vector_for_intersection_overlap_2_2));
	assert(theDo_intersect_2(circ_arc_overlap_3, circ_arc_overlap_1));
  assert(vector_for_intersection_overlap_2_2.size() == 1);
  assign(the_pair, vector_for_intersection_overlap_2_2[0]);
  std::cout << "x = " << the_pair.first.x() << " the result must be = " <<
    center_circ_intersection_2_1_x + circ_intersection_2_1_r * sqrt2 << std::endl;
  
  assert(the_pair.first.x() * 
	 (center_circ_intersection_2_1_x + circ_intersection_2_1_r * sqrt2) >= 0);
  assert(square(the_pair.first.x() - RT(center_circ_intersection_2_1_x))
	 == (circ_intersection_2_1_r * circ_intersection_2_1_r / typename CK::RT(2)));
  std::cout << "y = " << the_pair.first.y() << " the result must be = " <<
    center_circ_intersection_2_1_y - circ_intersection_2_1_r * sqrt2 << std::endl;
  
  assert(the_pair.first.y() * 
	 (center_circ_intersection_2_1_y - circ_intersection_2_1_r * sqrt2) >= 0);
  assert(square(the_pair.first.y() - RT(center_circ_intersection_2_1_y))
	 == (circ_intersection_2_1_r * circ_intersection_2_1_r / typename CK::RT(2)));

  std::cout << "Intersection : overlap in two points: " <<
    "lower_part_arc , upper_part_arc" << std::endl;
  std::vector< CGAL::Object > 
    vector_for_intersection_overlap_3_1;
  theConstruct_intersect_2(circ_arc_overlap_upper_part, 
			   circ_arc_overlap_lower_part,
			   std::back_inserter(vector_for_intersection_overlap_3_1));
  assert(theDo_intersect_2(circ_arc_overlap_upper_part, circ_arc_overlap_lower_part));
  assert(vector_for_intersection_overlap_3_1.size() == 2);

  assign(the_pair, vector_for_intersection_overlap_3_1[0]); 
  assert(the_pair.first == circ_arc_overlap_lower_part.source());
  //assert(the_pair.first.is_left());
  assign(the_pair, vector_for_intersection_overlap_3_1[1]); 
  assert(the_pair.first == circ_arc_overlap_lower_part.target());
  //assert(!the_pair.first.is_left());
  std::vector< CGAL::Object > 
    vector_for_intersection_overlap_3_2;
  theConstruct_intersect_2(circ_arc_overlap_lower_part, 
			   circ_arc_overlap_upper_part,
			   std::back_inserter(vector_for_intersection_overlap_3_2));
	assert(theDo_intersect_2(circ_arc_overlap_lower_part, circ_arc_overlap_upper_part));
  
  assert(vector_for_intersection_overlap_3_2.size() == 2);
  assign(the_pair, vector_for_intersection_overlap_3_2[0]);
  assert(the_pair.first == circ_arc_overlap_lower_part.source());
  //assert(the_pair.first.is_left());
  assign(the_pair, vector_for_intersection_overlap_3_2[1]);
  assert(the_pair.first == circ_arc_overlap_lower_part.target());
  //assert(!the_pair.first.is_left());

  //Intersection with 2 Circular_arc no x_monotone
  std::cout << "Intersection on two points of 2 Circular_arc no x_monotone" 
	    << std::endl;
  Circular_arc_2 circ_arc_no_x_monotone_1(circ_intersections_2_1,
					  line_arc_overlap_low_right, true,
					  line_arc_overlap_low_right, false);
  Point_2 center_circ_intersections_2_4(center_circ_intersection_2_1_x 
					- circ_intersection_2_1_r,
					center_circ_intersection_2_1_y 
					- circ_intersection_2_1_r);
  Circle_2 circ_intersections_2_4(center_circ_intersections_2_4,
				  circ_intersection_2_1_r 
				  * circ_intersection_2_1_r);
  std::vector< CGAL::Object > 
    vector_for_intersection_no_x_monotone_1_1;
  theConstruct_intersect_2(circ_arc_no_x_monotone_1, 
 			   circ_intersections_2_4,
 			   std::back_inserter(vector_for_intersection_no_x_monotone_1_1));
  assert(theDo_intersect_2(circ_arc_no_x_monotone_1, circ_intersections_2_4));
  assert(vector_for_intersection_no_x_monotone_1_1.size() == 2);
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_1_1[0]));
  assert(the_pair.first == CGAL::circle_intersect<CK>(circ_intersections_2_1,
						      circ_intersections_2_4,
						      true));
  assert(the_pair.second == 1u);
  //assert(the_pair.first.is_left());
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_1_1[1]));
  assert(the_pair.first == CGAL::circle_intersect<CK>(circ_intersections_2_1,
						      circ_intersections_2_4,
						      false));
  assert(the_pair.second == 1u);
  //assert(!the_pair.first.is_left());

  std::cout << "Intersection on one points no tangent of 2 Circular_arc no x_monotone" 
	    << std::endl;
  Point_2 center_circ_intersections_2_5(center_circ_intersection_2_1_x 
					- circ_intersection_2_1_r,
					center_circ_intersection_2_1_y );
  Circle_2 circ_intersections_2_5(center_circ_intersections_2_5,
				  circ_intersection_2_1_r 
				  * circ_intersection_2_1_r);
  std::vector< CGAL::Object > 
    vector_for_intersection_no_x_monotone_1_2;
  theConstruct_intersect_2(circ_arc_no_x_monotone_1, 
 			   circ_intersections_2_5,
 			   std::back_inserter(vector_for_intersection_no_x_monotone_1_2));
  assert(theDo_intersect_2(circ_arc_no_x_monotone_1, circ_intersections_2_5));
  assert(vector_for_intersection_no_x_monotone_1_2.size() == 1);
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_1_2[0]));
  assert(the_pair.first == CGAL::circle_intersect<CK>(circ_intersections_2_1,
						      circ_intersections_2_5,
						      true));
  assert(the_pair.second == 1u);
  //assert(the_pair.first.is_left());

  std::cout << "Intersection on one points tangent of 2 Circular_arc no x_monotone" 
	    << std::endl;
  Point_2 center_circ_intersections_2_6(center_circ_intersection_2_1_x ,
					center_circ_intersection_2_1_y 
					- 2 * circ_intersection_2_1_r );
  Circle_2 circ_intersections_2_6(center_circ_intersections_2_6,
				  circ_intersection_2_1_r 
				  * circ_intersection_2_1_r);
  std::vector< CGAL::Object > 
    vector_for_intersection_no_x_monotone_1_3;
  theConstruct_intersect_2(circ_arc_no_x_monotone_1, 
 			   circ_intersections_2_6,
 			   std::back_inserter(vector_for_intersection_no_x_monotone_1_3));
  assert(theDo_intersect_2(circ_arc_no_x_monotone_1, circ_intersections_2_6));
  assert(vector_for_intersection_no_x_monotone_1_3.size() == 1);
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_1_3[0]));
  assert(the_pair.first == CGAL::circle_intersect<CK>(circ_intersections_2_1,
						      circ_intersections_2_6,
						      true));
  assert(the_pair.second == 2u);
    


  Point_2 center_circ_intersections_2_7(center_circ_intersection_2_1_x
					  - 2 * circ_intersection_2_1_r,
					  center_circ_intersection_2_1_y);
  Circle_2 circ_intersections_2_7(center_circ_intersections_2_7,
				  circ_intersection_2_1_r 
				  * circ_intersection_2_1_r);
  std::vector< CGAL::Object > 
    vector_for_intersection_no_x_monotone_1_4;
  theConstruct_intersect_2(circ_arc_no_x_monotone_1, 
 			   circ_intersections_2_7,
 			   std::back_inserter(vector_for_intersection_no_x_monotone_1_4));
  assert(theDo_intersect_2(circ_arc_no_x_monotone_1, circ_intersections_2_7));
  assert(vector_for_intersection_no_x_monotone_1_4.size() == 1);
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_1_4[0]));
  assert(the_pair.first == CGAL::circle_intersect<CK>(circ_intersections_2_1,
						      circ_intersections_2_7,
						      true));
  assert(the_pair.second == 2u);
  

  Point_2 center_circ_intersections_2_8(center_circ_intersection_2_1_x
					+ circ_intersection_2_1_r,
					center_circ_intersection_2_1_y);
  Circle_2 circ_intersections_2_8(center_circ_intersections_2_8,
				  circ_intersection_2_1_r 
				  * circ_intersection_2_1_r * 4);
  std::vector< CGAL::Object > 
    vector_for_intersection_no_x_monotone_1_5;
  theConstruct_intersect_2(circ_arc_no_x_monotone_1, 
 			   circ_intersections_2_8,
 			   std::back_inserter(vector_for_intersection_no_x_monotone_1_5));
  assert(theDo_intersect_2(circ_arc_no_x_monotone_1, circ_intersections_2_8));
  assert(vector_for_intersection_no_x_monotone_1_5.size() == 1);
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_1_5[0]));
  assert(the_pair.first == CGAL::circle_intersect<CK>(circ_intersections_2_1,
						      circ_intersections_2_8,
						      true));
  assert(the_pair.second == 2u);


  Point_2 center_circ_intersections_2_9(center_circ_intersection_2_1_x,
					center_circ_intersection_2_1_y
					+ circ_intersection_2_1_r);
  Circle_2 circ_intersections_2_9(center_circ_intersections_2_9,
				  circ_intersection_2_1_r 
				  * circ_intersection_2_1_r * 4);
  std::vector< CGAL::Object > 
    vector_for_intersection_no_x_monotone_1_6;
  theConstruct_intersect_2(circ_arc_no_x_monotone_1, 
 			   circ_intersections_2_9,
 			   std::back_inserter(vector_for_intersection_no_x_monotone_1_6));
  assert(theDo_intersect_2(circ_arc_no_x_monotone_1, circ_intersections_2_9));
  assert(vector_for_intersection_no_x_monotone_1_6.size() == 1);
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_1_6[0]));
  assert(the_pair.first == CGAL::circle_intersect<CK>(circ_intersections_2_1,
						      circ_intersections_2_9,
						      true));
  assert(the_pair.second == 2u);
  
  
  std::cout << "Intersection of 2 Circular_arc no x_monotone" 
	    << " : overlap on a circular arc in bottom " 
	    << std::endl;
  Circular_arc_2 circ_arc_no_x_monotone_2(circ_intersections_2_1,
					  line_arc_overlap_low_left, true,
					  line_arc_overlap_low_left, false);
  std::vector< CGAL::Object > 
    vector_for_intersection_no_x_monotone_2_1;
  theConstruct_intersect_2(circ_arc_no_x_monotone_1, 
 			   circ_arc_no_x_monotone_2,
 			   std::back_inserter(vector_for_intersection_no_x_monotone_2_1));
  assert(theDo_intersect_2(circ_arc_no_x_monotone_1, circ_arc_no_x_monotone_2));
  assert(vector_for_intersection_no_x_monotone_2_1.size() == 1);
  assert(assign(circ_arc_in_overlap, vector_for_intersection_no_x_monotone_2_1[0]));
  assert(circ_arc_in_overlap.is_x_monotone());
  assert(circ_arc_in_overlap.source() == circ_arc_overlap_result.source());
  assert(circ_arc_in_overlap.target() == circ_arc_overlap_result.target());
  
  
  std::cout << "Intersection of 2 Circular_arc no x_monotone" 
	    << " : overlap on a circular arc at left " 
	    << std::endl;
  Circular_arc_2 circ_arc_no_x_monotone_3(circ_intersections_2_1,
 					  line_arc_overlap_low_left, false,
 					  line_arc_overlap_low_left, true);
  
  
  std::vector< CGAL::Object > 
    vector_for_intersection_no_x_monotone_2_2;
  theConstruct_intersect_2(circ_arc_no_x_monotone_1, 
 			   circ_arc_no_x_monotone_3,
 			   std::back_inserter(vector_for_intersection_no_x_monotone_2_2));
  assert(theDo_intersect_2(circ_arc_no_x_monotone_1, circ_arc_no_x_monotone_3));
  assert(vector_for_intersection_no_x_monotone_2_2.size() == 1);
  assert(assign(circ_arc_in_overlap, vector_for_intersection_no_x_monotone_2_2[0]));
  assert(circ_arc_in_overlap.source() == circ_arc_no_x_monotone_1.source());
  assert(circ_arc_in_overlap.target() == circ_arc_no_x_monotone_3.target());
  
  
  std::cout << "Intersection of 2 Circular_arc no x_monotone" 
	    << " : overlap on a circular arc in bottom "
	    << "and one endpoint"
	    << std::endl;
  Circular_arc_2 circ_arc_no_x_monotone_4(circ_intersections_2_1,
 					  line_arc_overlap_low_left, true,
 					  line_arc_overlap_low_right, true);
  std::vector< CGAL::Object > 
    vector_for_intersection_no_x_monotone_2_3;
  theConstruct_intersect_2(circ_arc_no_x_monotone_1, 
 			   circ_arc_no_x_monotone_4,
 			   std::back_inserter(vector_for_intersection_no_x_monotone_2_3));
  assert(theDo_intersect_2(circ_arc_no_x_monotone_1, circ_arc_no_x_monotone_4));
  std::cout << vector_for_intersection_no_x_monotone_2_3.size() << std::endl;

  std::cout << vector_for_intersection_no_x_monotone_2_3.size() << std::endl;

  assert(vector_for_intersection_no_x_monotone_2_3.size() == 2);
  assert(assign(circ_arc_in_overlap,
		vector_for_intersection_no_x_monotone_2_3[0]));
  assert(assign(the_pair, 
		vector_for_intersection_no_x_monotone_2_3[1]));
  assert(circ_arc_in_overlap.is_x_monotone());
  assert(circ_arc_in_overlap.source() 
	 == circ_arc_overlap_result.source());
  assert(circ_arc_in_overlap.target() 
	 == circ_arc_overlap_result.target());
  assert(the_pair.first == circ_arc_no_x_monotone_4.target());
  assert(the_pair.second == 1u);
  

  std::cout << "Intersection of 2 Circular_arc no x_monotone"
	     << ": overlap in two points" << std::endl;
  
  Circular_arc_2 circ_arc_no_x_monotone_5(circ_intersections_2_1,
					  line_arc_overlap_low_right, true,
					  line_arc_overlap_low_left, true);
  std::vector< CGAL::Object > 
    vector_for_intersection_no_x_monotone_2_4;
  theConstruct_intersect_2(circ_arc_no_x_monotone_4, 
			   circ_arc_no_x_monotone_5,
			    std::back_inserter(vector_for_intersection_no_x_monotone_2_4));
	assert(theDo_intersect_2(circ_arc_no_x_monotone_4, circ_arc_no_x_monotone_5));
  std::cout << vector_for_intersection_no_x_monotone_2_4.size() << std::endl;
  assert(vector_for_intersection_no_x_monotone_2_4.size() == 2);
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_2_4[0]));
  assert(the_pair.first == circ_arc_no_x_monotone_5.target());
  assert(the_pair.second == 1u);
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_2_4[1]));
  assert(the_pair.first == circ_arc_no_x_monotone_5.source());
  assert(the_pair.second == 1u);


  std::cout << "Intersection of 2 Circular_arc no x_monotone"
	    << ": overlap in one points" << std::endl;
  Circular_arc_2 circ_arc_no_x_monotone_6(circ_intersections_2_1,
					  line_arc_overlap_low_right, false,
					  line_arc_overlap_low_right, true);
  std::vector< CGAL::Object > 
    vector_for_intersection_no_x_monotone_2_5;
  theConstruct_intersect_2(circ_arc_no_x_monotone_6, 
			   circ_arc_no_x_monotone_5,
			    std::back_inserter(vector_for_intersection_no_x_monotone_2_5));
	assert(theDo_intersect_2(circ_arc_no_x_monotone_6, circ_arc_no_x_monotone_5));
  std::cout << vector_for_intersection_no_x_monotone_2_5.size() << std::endl;
  assert(vector_for_intersection_no_x_monotone_2_5.size() == 1);
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_2_5[0]));
  assert(the_pair.first == circ_arc_no_x_monotone_5.source());
  assert(the_pair.second == 1u);

  
  std::cout << "Intersection of 2 Circular_arc no x_monotone" 
	    << " : overlap on 2 circular arcs " << std::endl; 
  Circular_arc_2 circ_arc_no_x_monotone_7(circ_intersections_2_1,
					  line_arc_overlap_low_left, false,
					  line_arc_overlap_low_right, false);
  std::vector< CGAL::Object > 
    vector_for_intersection_no_x_monotone_2_6;
  theConstruct_intersect_2(circ_arc_no_x_monotone_7, 
			   circ_arc_no_x_monotone_4,
			    std::back_inserter(vector_for_intersection_no_x_monotone_2_6));
	assert(theDo_intersect_2(circ_arc_no_x_monotone_7, circ_arc_no_x_monotone_4));
  std::cout << vector_for_intersection_no_x_monotone_2_6.size() << std::endl;
  assert(vector_for_intersection_no_x_monotone_2_6.size() == 2);
  assign(circ_arc_in_overlap,vector_for_intersection_no_x_monotone_2_6[0]);
  assign(circ_arc_in_overlap_2,vector_for_intersection_no_x_monotone_2_6[1]);
  assert((circ_arc_in_overlap.source() == circ_arc_no_x_monotone_7.source() &&
          circ_arc_in_overlap.target() == circ_arc_no_x_monotone_4.target()) ||
         (circ_arc_in_overlap_2.source() == circ_arc_no_x_monotone_7.source() &&
          circ_arc_in_overlap_2.target() == circ_arc_no_x_monotone_4.target()));
  std::cout << "source4 : " << std::endl
	    << circ_arc_no_x_monotone_4.source() << std::endl
	    << "target4 : " << std::endl
	    << circ_arc_no_x_monotone_4.target() << std::endl
	    << "source7 : " << std::endl
	    << circ_arc_no_x_monotone_7.source() << std::endl
	    << "target7 : " << std::endl
	    << circ_arc_no_x_monotone_7.target() << std::endl;
  std::cout << "res source : " << std::endl
	    << circ_arc_in_overlap.source() << std::endl 
	    << "res target : " << std::endl
	    << circ_arc_in_overlap.target() << std::endl;
  assert(circ_arc_in_overlap.is_x_monotone());
  std::cout << "res source : " << std::endl
	    << circ_arc_in_overlap.source() << std::endl 
	    << "res target : " << std::endl
	    << circ_arc_in_overlap.target() << std::endl;
  if(circ_arc_in_overlap.source() == circ_arc_no_x_monotone_7.source() &&
          circ_arc_in_overlap.target() == circ_arc_no_x_monotone_4.target()) {
    assert(circ_arc_in_overlap_2.source() == circ_arc_no_x_monotone_4.source());
    assert(circ_arc_in_overlap_2.target() == circ_arc_no_x_monotone_7.target());
  } else {
    assert(circ_arc_in_overlap.source() == circ_arc_no_x_monotone_4.source());
    assert(circ_arc_in_overlap.target() == circ_arc_no_x_monotone_7.target());
  }
  
  //Make_x_monotone_2 with a full circle
  Make_x_monotone_2 theMake_x_monotone = ck.make_x_monotone_2_object();
  int x = theRandom.get_int(random_min,random_max);
  int y = theRandom.get_int(random_min,random_max);
  int r = theRandom.get_int(1,random_max);
  Point_2 center_circ_monotone(x,y);
  Circle_2 circ_monotone(center_circ_monotone, r*r);
  Circular_arc_2 theCircular_arc_2_full(circ_monotone);
  std::vector< CGAL::Object > outputIterator1, outputIterator1l;
  theMake_x_monotone(theCircular_arc_2_full,
		     std::back_inserter(outputIterator1));
	make_x_monotone(theCircular_arc_2_full,	std::back_inserter(outputIterator1l));
  std::cout << std::endl;
  std::cout << "x_monotone full circle : " 
	    << circ_monotone << std::endl;
  Circular_arc_2 circular_arc_2_full, circular_arc_2_fulll;
  for(std::size_t i = 0; i < outputIterator1.size(); i++){
    assign(circular_arc_2_full,  outputIterator1[i]);
    assign(circular_arc_2_fulll,  outputIterator1l[i]);
    std::cout << "Circular_arc_2_" << i << " : " 
	      << circular_arc_2_full << std::endl;
   std::cout << "Circular_arc_2_" << i << "source : " 
	      << circular_arc_2_full.source() << std::endl;
    std::cout << "Circular_arc_2_" << i << "target : " 
	      << circular_arc_2_full.target() << std::endl;
    assert(circular_arc_2_full.is_x_monotone());
		assert(circular_arc_2_full == circular_arc_2_fulll);
  }

  //Make_xy_monotone_2 with a full circle
  Make_xy_monotone_2 theMake_xy_monotone = ck.make_xy_monotone_2_object();
  outputIterator1.clear();
  outputIterator1l.clear();
  theMake_xy_monotone(theCircular_arc_2_full,
		     std::back_inserter(outputIterator1));
	theMake_xy_monotone(theCircular_arc_2_full,
		     std::back_inserter(outputIterator1l));
	assert(outputIterator1.size() == 4);
  for(std::size_t i = 0; i < outputIterator1.size(); i++){
    assign(circular_arc_2_full,  outputIterator1[i]);
    assign(circular_arc_2_fulll,  outputIterator1l[i]);
    assert(circular_arc_2_full.is_x_monotone());
    assert(circular_arc_2_full.is_y_monotone());
    assert(circular_arc_2_full == circular_arc_2_fulll);
  } 

  //Make_xy_monotone_2 general test
	Point_2 ps[8];
	ps[0] = Point_2(-5, 0);
  ps[1] = Point_2(-3, -4);
  ps[2] = Point_2(0, -5);
  ps[3] = Point_2(3, -4);
  ps[4] = Point_2(5, 0);
  ps[5] = Point_2(3, 4);
  ps[6] = Point_2(0, 5);
  ps[7] = Point_2(-3, 4);
	Circle_2 tc = Circle_2(Point_2(0,0),25);
	
	unsigned isize[2][8];
	isize[0][1] = 1;
	isize[0][2] = 1;
	isize[0][3] = 2;
	isize[0][4] = 2;
	isize[0][5] = 3;
	isize[0][6] = 3;
	isize[0][7] = 4;
	isize[1][2] = 1;
	isize[1][3] = 2;
	isize[1][4] = 2;
	isize[1][5] = 3;
	isize[1][6] = 3;
	isize[1][7] = 4;
	isize[1][0] = 4; 
	
	for(int i=0; i<2; i++) {
		for(int j=i+1; j!=i; j = (j+1)%8) {
	    Circular_arc_2 ca;
	    ca = Circular_arc_2(tc, ps[i], ps[j]); 
      outputIterator1.clear();
      theMake_xy_monotone(ca, std::back_inserter(outputIterator1));
			std::cout << "T: " << i << " " << j << std::endl;
	    assert(outputIterator1.size() == isize[i][j]);
      for(std::size_t k = 0; k < outputIterator1.size(); k++) {
        assign(circular_arc_2_full,  outputIterator1[k]);
        assert(circular_arc_2_full.is_x_monotone());
        assert(circular_arc_2_full.is_y_monotone());
      }
    }
  }

  //Make_x_monotone_2 with a three quarter of last circle
  Point_2 pointLine_2_1(x,y+r);
  Line_2  theLine_2_1(center_circ_monotone, pointLine_2_1);
  Point_2 pointLine_2_2(x+r,y);
  Line_2  theLine_2_2(center_circ_monotone, pointLine_2_2);
  Circular_arc_2 theCircular_arc_2_quarter(circ_monotone,
					   theLine_2_1, true,
					   theLine_2_2, true);
  std::vector< CGAL::Object > vector_of_object_1;
  theMake_x_monotone(theCircular_arc_2_quarter,
		     std::back_inserter(vector_of_object_1));
  std::cout << std::endl;
  std::cout << "x_monotone a three quarter of last circle: "
	    << circ_monotone << std::endl;
  std::cout << vector_of_object_1.size() << std::endl;
  Circular_arc_2 circular_arc_2_quarter;
  for(std::size_t i = 0; i < vector_of_object_1.size(); i++){
    assign(circular_arc_2_quarter,  vector_of_object_1[i]);
    std::cout << "Circular_arc_2_" << i << " : " 
	      << circular_arc_2_quarter << std::endl;
    assert(circular_arc_2_quarter.is_x_monotone());
  }

  //Make_x_monotone_2 with half circle
  Circular_arc_2 theCircular_arc_2_half(circ_monotone,
					theLine_2_2, false,
					theLine_2_2, true);
   std::vector< CGAL::Object > vector_of_object_1_half;
  theMake_x_monotone(theCircular_arc_2_half,
		     std::back_inserter(vector_of_object_1_half));
  std::cout << std::endl;
  std::cout << "x_monotone a half circle" << std::endl;
  assert(vector_of_object_1_half.size() == 1);
  assign(circular_arc_2_quarter,  vector_of_object_1_half[0]);
  assert(circular_arc_2_quarter.is_x_monotone());
  assert(circular_arc_2_quarter.source() == theCircular_arc_2_half.source());
  assert(circular_arc_2_quarter.target() == theCircular_arc_2_half.target());					    


  //Make_x_monotone_2 with a random circular arc
  int pointLine_2_3_x = theRandom.get_int(random_min,random_max);
  int pointLine_2_3_y = theRandom.get_int(random_min,random_max);
  while((pointLine_2_3_x == x) && (pointLine_2_3_y == y)){
    if(pointLine_2_3_x == x) 
      pointLine_2_3_x = theRandom.get_int(random_min,random_max);
    else 
      pointLine_2_3_y = theRandom.get_int(random_min,random_max);
  }
  Point_2 pointLine_2_3(pointLine_2_3_x, pointLine_2_3_y);
  Line_2  theLine_2_3(center_circ_monotone, pointLine_2_3);
  int pointLine_2_4_x = theRandom.get_int(random_min,random_max);
  int pointLine_2_4_y = theRandom.get_int(random_min,random_max);
  while((pointLine_2_4_x == x) && (pointLine_2_4_y == y)){
    if(pointLine_2_4_x == x) 
      pointLine_2_4_x = theRandom.get_int(random_min,random_max);
    else 
      pointLine_2_4_y = theRandom.get_int(random_min,random_max);
  }
  Point_2 pointLine_2_4(pointLine_2_4_x, pointLine_2_4_y);
  Line_2  theLine_2_4(center_circ_monotone, pointLine_2_4);
  Circular_arc_2 theCircular_arc_2_random(circ_monotone,
					  theLine_2_3, true,
					  theLine_2_4, true);
  std::vector< CGAL::Object > vector_of_object_2;
  theMake_x_monotone(theCircular_arc_2_random,
		     std::back_inserter(vector_of_object_2));
  std::cout << std::endl;
  std::cout << "x_monotone random circular arc: " 
	    << circ_monotone << std::endl;
  Circular_arc_2 circular_arc_2_random;
  for(std::size_t i = 0; i < vector_of_object_2.size(); i++){
    assign(circular_arc_2_random,  vector_of_object_2[i]);
    std::cout << "Circular_arc_2_" << i << " : "
	      << circular_arc_2_random << std::endl;
    assert(circular_arc_2_random.is_x_monotone());
  }

	std::cout << "Split_2_object " << std::endl;
	//we make the circle1
	int center1_x = theRandom.get_int(random_min, random_max);
	int center1_y = theRandom.get_int(random_min, random_max);
	Point_2 center1(center1_x,center1_y); 
	int circ1_r = theRandom.get_int(1, random_max);
	Circle_2 circ1(center1, circ1_r * circ1_r);
	Point_2 center1_low_right(center1_x + circ1_r, center1_y - circ1_r);
	Circle_2 circ1_low_right(center1_low_right, circ1_r * circ1_r);
	Point_2 center1_low_left(center1_x - circ1_r, center1_y - circ1_r);
	Circle_2 circ1_low_left(center1_low_left, circ1_r * circ1_r);
	Point_2 point_2_left(center1_x - circ1_r, center1_y);
	Line_2 theLine_2_horizontal(center1, point_2_left);
	//The circ1_arc_high and circ1_arc_low are x_monotone
	Circular_arc_2 circ1_arc_low(circ1,
	                              theLine_2_horizontal,true,
	                              theLine_2_horizontal, false);
	//p1 is lefter and lower than p2
	Circular_arc_point_2 circ1_arc_end_p1 =
	           CGAL::circle_intersect<CK>(circ1, circ1_low_right, true);
	Split_2 theSplit_2 = ck.split_2_object();
	Circular_arc_2 circ_arc_split_1;
	Circular_arc_2 circ_arc_split_2;
	theSplit_2(circ1_arc_low, circ1_arc_end_p1,
	           circ_arc_split_1, circ_arc_split_2);
	assert(circ_arc_split_1.target() == circ1_arc_end_p1);
	assert(circ1_arc_low.source() == circ_arc_split_1.source());
	assert(circ_arc_split_1.target() == circ_arc_split_2.source());
	assert(circ1_arc_low.target() == circ_arc_split_2.target());
	
	//We used a point created without the support circle
	Circular_arc_point_2 circ1_arc_end_p2 =
	            CGAL::circle_intersect<CK>(circ1_low_left, circ1_low_right, true);
	theSplit_2(circ1_arc_low, circ1_arc_end_p2,
	            circ_arc_split_1, circ_arc_split_2);
	assert(circ_arc_split_1.target() == circ1_arc_end_p2);
	assert(circ1_arc_low.source() == circ_arc_split_1.source());
	assert(circ_arc_split_1.target() == circ_arc_split_2.source());
	assert(circ1_arc_low.target() == circ_arc_split_2.target());
  
  //The commented code in bottom must create an error
  ////We used a point which is not on the arc
  //Circular_arc_2 arc_aux(circ1_low_right,
  //			 Line_2(center1, center1_low_right),true,
  //			 Line_2(center1, center1_low_right),false);
  //theSplit_2(circ1_arc_low, arc_aux.source(),
  //	     circ_arc_split_1, circ_arc_split_2);

  // testing intersect_2(Line_2, Circular_arc_2)
  Line_2 lo1 = Line_2(Point_2(0,0), Point_2(0,10));
  Circular_arc_2 cao1 = Circular_arc_2(Point_2(0,0), Point_2(-5,5), Point_2(0,10)); // = two intersection p/arc
  Circular_arc_2 cao2 = Circular_arc_2(Point_2(0,0), Point_2(5,5), Point_2(0,10));
  Circular_arc_2 cao3 = Circular_arc_2(Point_2(-5,5), Point_2(0,0), Point_2(5,5)); // = one intersection p/ arc
  Circular_arc_2 cao4 = Circular_arc_2(Point_2(-5,5), Point_2(0,10), Point_2(5,5));
  Circular_arc_2 cao5 = Circular_arc_2(Point_2(1,1), Point_2(-1,2), Point_2(1,4)); // = zero-two intersections
  Circular_arc_2 cao6 = Circular_arc_2(Point_2(1,4), Point_2(-1,2), Point_2(1,1)); 
  Circular_arc_2 cao7 = Circular_arc_2(Point_2(10,10), Point_2(0,0), Point_2(10,-10)); // = tangency
  Circular_arc_2 cao8 = Circular_arc_2(Point_2(11,10), Point_2(1,0), Point_2(11,-10)); // = no intersection

  std::cout << "Testing intersect with lines" << std::endl;
  std::vector< CGAL::Object > v_ll1, v_ll2, v_ll3, v_ll4, v_ll5, v_ll6, v_ll7, v_ll8;
  theConstruct_intersect_2(lo1, cao1, std::back_inserter(v_ll1));
  theConstruct_intersect_2(lo1, cao2, std::back_inserter(v_ll2));
  theConstruct_intersect_2(lo1, cao3, std::back_inserter(v_ll3));
  theConstruct_intersect_2(lo1, cao4, std::back_inserter(v_ll4));
  theConstruct_intersect_2(lo1, cao5, std::back_inserter(v_ll5));
  theConstruct_intersect_2(lo1, cao6, std::back_inserter(v_ll6));
  theConstruct_intersect_2(lo1, cao7, std::back_inserter(v_ll7));
  theConstruct_intersect_2(lo1, cao8, std::back_inserter(v_ll8));
  assert(v_ll1.size() == 2);
  assert(theDo_intersect_2(lo1, cao1));
  assert(v_ll2.size() == 2);
  assert(theDo_intersect_2(lo1, cao2));
  assert(v_ll3.size() == 1);
  assert(theDo_intersect_2(lo1, cao3));
  assert(v_ll4.size() == 1);
  assert(theDo_intersect_2(lo1, cao4));
  assert(v_ll5.size() == 0);
  assert(!theDo_intersect_2(lo1, cao5));
  assert(v_ll6.size() == 2);
  assert(theDo_intersect_2(lo1, cao6));
  assert(v_ll7.size() == 1);
  assert(theDo_intersect_2(lo1, cao7));
  assert(v_ll8.size() == 0);
  assert(!theDo_intersect_2(lo1, cao8));

	Line_arc_2 llu1 = Line_arc_2(Point_2(-1,-1), Point_2(1,1));
	Line_arc_2 llu2 = Line_arc_2(Point_2(-1,-1), Point_2(-1,1));
	Line_arc_2 llu3 = Line_arc_2(Point_2(-2,-1), Point_2(-2,1));
	Circle_2 ccu = Circle_2(Point_2(0,-1), Point_2(-1,0), Point_2(0,1));
	
	std::vector< CGAL::Object > v_llc1, v_llc2, v_llc3;
  theConstruct_intersect_2(llu1, ccu, std::back_inserter(v_llc1));
  theConstruct_intersect_2(llu2, ccu, std::back_inserter(v_llc2));
  theConstruct_intersect_2(llu3, ccu, std::back_inserter(v_llc3));

  assert(v_llc1.size() == 2);
  assert(theDo_intersect_2(llu1, ccu));
  assert(CGAL::do_intersect(llu1, ccu));
  assert(v_llc2.size() == 1);
  assert(theDo_intersect_2(llu2, ccu));
  assert(CGAL::do_intersect(llu2, ccu));
  assert(v_llc3.size() == 0);
  assert(!theDo_intersect_2(llu3, ccu));
  assert(!CGAL::do_intersect(llu3, ccu));

	std::vector< CGAL::Object > v_rllc1, v_rllc2, v_rllc3, v_rllc1l, v_rllc2l, v_rllc3l;
  theConstruct_intersect_2(llu1.supporting_line(), ccu, std::back_inserter(v_rllc1));
  theConstruct_intersect_2(llu2.supporting_line(), ccu, std::back_inserter(v_rllc2));
  theConstruct_intersect_2(llu3.supporting_line(), ccu, std::back_inserter(v_rllc3));
  theConstruct_intersect_2(ccu, llu1.supporting_line(), std::back_inserter(v_rllc1l));
  theConstruct_intersect_2(ccu, llu2.supporting_line(), std::back_inserter(v_rllc2l));
  theConstruct_intersect_2(ccu, llu3.supporting_line(), std::back_inserter(v_rllc3l));

  assert(v_rllc1.size() == 2);
  assert(theDo_intersect_2(llu1.supporting_line(), ccu));
  assert(CGAL::do_intersect(llu1.supporting_line(), ccu));

  assert(v_rllc1l.size() == 2);
  assert(theDo_intersect_2(ccu, llu1.supporting_line()));
  assert(CGAL::do_intersect(ccu, llu1.supporting_line()));

  assert(v_rllc2.size() == 1);
  assert(theDo_intersect_2(llu2.supporting_line(), ccu));
  assert(CGAL::do_intersect(llu2.supporting_line(), ccu));

  assert(v_rllc2l.size() == 1);
  assert(theDo_intersect_2(ccu, llu2.supporting_line()));
  assert(CGAL::do_intersect(ccu, llu2.supporting_line()));

  assert(v_rllc3.size() == 0);
  assert(!theDo_intersect_2(llu3.supporting_line(), ccu));
  assert(!CGAL::do_intersect(llu3.supporting_line(), ccu));

  assert(v_rllc3l.size() == 0);
  assert(!theDo_intersect_2(ccu, llu3.supporting_line()));
  assert(!CGAL::do_intersect(ccu, llu3.supporting_line()));

  // TEST THE FUNCTOR CALL (VC8 porting mainly reason)
	Circular_arc_2 ccaa = 
	  typename CK::Construct_circular_arc_2()(Point_2(1, 2), Point_2(2, 2), Point_2(3, 3));
	Line_arc_2 llaa = 
	  typename CK::Construct_line_arc_2()(Point_2(1, 2), Point_2(2, 2));
	Circular_arc_point_2 ccaapp = typename CK::Construct_circular_arc_point_2()(Point_2(1, 2));
	typename CK::Construct_circular_min_vertex_2()(llaa);
	typename CK::Construct_circular_max_vertex_2()(llaa);
	typename CK::Construct_circular_source_vertex_2()(llaa);
	typename CK::Construct_circular_target_vertex_2()(llaa);
	
#ifndef CGAL_NO_DEPRECATED_CODE
	// testing the deprecate stuffs
	typename CK::Construct_supporting_circle_2()(ccaa);
	typename CK::Construct_supporting_line_2()(llaa);
#endif
}
