#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Circular_kernel.h>
#include <CGAL/Circular_arc_traits.h>
#include <CGAL/Algebraic_kernel_2_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Curved_kernel/function_objects_polynomial_circular.h>
#include <CGAL/Curved_kernel/Circular_arc_2.h>

#include <CGAL/NT_extensions_Root_of/CGAL_Quotient.h>
#include <CGAL/NT_extensions_Root_of/CGAL_Gmpq.h>


#include <CGAL/Random.h>

template <class CK>
void _test_circle_construct(CK ck)
{
  typedef typename CK::Circle_2                    Circle_2;
  typedef typename CK::Circular_arc_2              Circular_arc_2;
  typedef typename CK::Point_2                     Point_2;
  typedef typename CK::Line_2                      Line_2;
  typedef typename CK::Circular_arc_endpoint_2     Circular_arc_endpoint_2;
  typedef typename CK::Construct_circle_2          Construct_circle_2;
  typedef typename CK::Construct_intersections_2   Construct_intersections_2;
  typedef typename CK::Make_x_monotone_2           Make_x_monotone_2;
  typedef typename CK::Split_2                     Split_2;  
  typedef typename CK::Get_equation                Get_equation;
  typedef typename CK::Polynomial_for_circles_2_2  Polynomial_for_circles_2_2;
  typedef typename CK::Compare_xy_2                Compare_xy_2;

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127;
  
  //test of get_equation_object()
  int x_equation = theRandom.get_int(random_min,random_max);
  int y_equation = theRandom.get_int(random_min,random_max);
  int r_equation = theRandom.get_int(1,random_max);
  Point_2 center_circ_equation(x_equation,y_equation);
  Circle_2 circ_equation(center_circ_equation, r_equation);
  std::cout << "the circle used by the equation :" 
	    << circ_equation << std::endl;
  Get_equation theEquation = ck.get_equation_object();
  Polynomial_for_circles_2_2 theResult_equation = theEquation(circ_equation);
  std::cout << "a= " << theResult_equation.a() << ", b= " <<
    theResult_equation.b() << ", r_sq= " <<
    theResult_equation.r_sq() <<std::endl;



  //Construct_circle_2
  Construct_circle_2 theConstruct_circle = ck.construct_circle_2_object();
  //We use the Polynomial_for_circles_2_2 fund before
  Circle_2 theConstruct_circle_2 = theConstruct_circle(theResult_equation);
  std::cout << "the circle made with the equation :" <<
    theConstruct_circle_2 << std::endl;
  assert(circ_equation == theConstruct_circle_2);


  //not declare 29/06/2005
  //ck.construct_circular_arc_2_object();
  //ck.construct_circular_arc_endpoint_2_object();

  //Constuct_intersections_2 with 2 intersection's points
  std::cout << std::endl << "construct_intersection_2" << std::endl;
  Construct_intersections_2 theConstruct_intersect_2 
    = ck.construct_intersections_2_object();
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
    vector_for_intersection_1;
  
  theConstruct_intersect_2(circ_intersections_2_1, 
			   circ_intersections_2_2,
			   std::back_inserter(vector_for_intersection_1));
  // there are 2 intersection's points
  std::pair<Circular_arc_endpoint_2, std::size_t> the_pair;
  assign(the_pair, vector_for_intersection_1[0]);
  Circular_arc_endpoint_2 first = the_pair.first;
  std::cout << first << std::endl;
  assign(the_pair, vector_for_intersection_1[1]);
  Circular_arc_endpoint_2 second = the_pair.first;
  std::cout << second << std::endl;
  assert((first.is_left() && !second.is_left()) 
	 || (!first.is_left() && second.is_left()));
  Compare_xy_2 theCompare_xy_2 = ck.compare_xy_2_object();
  assert(theCompare_xy_2(first, second) == CGAL::SMALLER);

  //Constuct_intersections_2 with 1 intersection's point
  Point_2 center_circ_intersections_2_3(center_circ_intersection_2_1_x 
					+ 2 * circ_intersection_2_1_r,
					center_circ_intersection_2_1_y);
  Circle_2 circ_intersections_2_3(center_circ_intersections_2_3,
				  circ_intersection_2_1_r 
				  * circ_intersection_2_1_r);
  Circular_arc_endpoint_2 the_intersection_point_1(circ_intersections_2_1,
						   circ_intersections_2_3,
						   true);
  Circular_arc_endpoint_2 the_intersection_point_2(circ_intersections_2_1,
						   circ_intersections_2_3,
						   false); 
  std::vector< CGAL::Object > 
    vector_for_intersection_2;
  theConstruct_intersect_2(circ_intersections_2_1, 
			   circ_intersections_2_3,
			   std::back_inserter(vector_for_intersection_2));
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
  std::vector< CGAL::Object > 
    vector_for_intersection_3;
  theConstruct_intersect_2(circ_arc_2_1_part_low, 
			   circ_arc_2_1_low_part_high,
			   std::back_inserter(vector_for_intersection_3));
  assert(vector_for_intersection_3.size() == 2);
  assign(the_pair, vector_for_intersection_3[0]);
  assert(the_pair.second == 1u);
  assert(the_pair.first == Circular_arc_endpoint_2(circ_intersections_2_1, 
						   circ_intersections_2_1_low, 
						   true));
  assign(the_pair, vector_for_intersection_3[1]);
  assert(the_pair.second == 1u);
  assert(the_pair.first == Circular_arc_endpoint_2(circ_intersections_2_1, 
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
  std::cout << "Intersection : overlap on a circular arc" << std::endl;
  Circular_arc_2 circ_arc_in_overlap;
  std::vector< CGAL::Object > 
    vector_for_intersection_overlap_1_1;
  theConstruct_intersect_2(circ_arc_overlap_2, 
 			   circ_arc_overlap_1,
 			   std::back_inserter(vector_for_intersection_overlap_1_1));
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
  assert(vector_for_intersection_overlap_2_1.size() == 1);
  assign(the_pair, vector_for_intersection_overlap_2_1[0]);
  std::cout << "x = " << the_pair.first.x() << " the result must be = " <<
    center_circ_intersection_2_1_x + circ_intersection_2_1_r * sqrt(2)/2 
	    << std::endl;
  
  assert(the_pair.first.x() * 
	 (center_circ_intersection_2_1_x + circ_intersection_2_1_r * sqrt(2)/2) >= 0);
  assert(square(the_pair.first.x() - center_circ_intersection_2_1_x)
	 == (circ_intersection_2_1_r * circ_intersection_2_1_r / typename CK::RT(2)));
  std::cout << "y = " << the_pair.first.y() << " the result must be = " <<
    center_circ_intersection_2_1_y - circ_intersection_2_1_r * sqrt(2)/2 
	    << std::endl;

  assert(the_pair.first.y() * 
	 (center_circ_intersection_2_1_y - circ_intersection_2_1_r * sqrt(2)/2) >= 0);
  assert(square(the_pair.first.y() - center_circ_intersection_2_1_y)
	 == (circ_intersection_2_1_r * circ_intersection_2_1_r / typename CK::RT(2)));

  std::vector< CGAL::Object > 
    vector_for_intersection_overlap_2_2;
  theConstruct_intersect_2(circ_arc_overlap_3, 
			   circ_arc_overlap_1,
			   std::back_inserter(vector_for_intersection_overlap_2_2));
  assert(vector_for_intersection_overlap_2_2.size() == 1);
  assign(the_pair, vector_for_intersection_overlap_2_2[0]);
  std::cout << "x = " << the_pair.first.x() << " the result must be = " <<
    center_circ_intersection_2_1_x + circ_intersection_2_1_r * sqrt(2)/2 << std::endl;
  
  assert(the_pair.first.x() * 
	 (center_circ_intersection_2_1_x + circ_intersection_2_1_r * sqrt(2)/2) >= 0);
  assert(square(the_pair.first.x() - center_circ_intersection_2_1_x)
	 == (circ_intersection_2_1_r * circ_intersection_2_1_r / typename CK::RT(2)));
  std::cout << "y = " << the_pair.first.y() << " the result must be = " <<
    center_circ_intersection_2_1_y - circ_intersection_2_1_r * sqrt(2)/2 << std::endl;
  
  assert(the_pair.first.y() * 
	 (center_circ_intersection_2_1_y - circ_intersection_2_1_r * sqrt(2)/2) >= 0);
  assert(square(the_pair.first.y() - center_circ_intersection_2_1_y)
	 == (circ_intersection_2_1_r * circ_intersection_2_1_r / typename CK::RT(2)));


  

  std::cout << "Intersection : overlap in two points: " <<
    "lower_part_arc , upper_part_arc" << std::endl;
  std::vector< CGAL::Object > 
    vector_for_intersection_overlap_3_1;
  theConstruct_intersect_2(circ_arc_overlap_upper_part, 
			   circ_arc_overlap_lower_part,
			   std::back_inserter(vector_for_intersection_overlap_3_1));
  
  assert(vector_for_intersection_overlap_3_1.size() == 2);
  assign(the_pair, vector_for_intersection_overlap_3_1[0]);
  assert(the_pair.first == circ_arc_overlap_lower_part.source());
  assert(the_pair.first.is_left());
  assign(the_pair, vector_for_intersection_overlap_3_1[1]);
  assert(the_pair.first == circ_arc_overlap_lower_part.target());
  assert(!the_pair.first.is_left());
  
  std::vector< CGAL::Object > 
    vector_for_intersection_overlap_3_2;
  theConstruct_intersect_2(circ_arc_overlap_lower_part, 
			   circ_arc_overlap_upper_part,
			   std::back_inserter(vector_for_intersection_overlap_3_2));
  
  assert(vector_for_intersection_overlap_3_2.size() == 2);
  assign(the_pair, vector_for_intersection_overlap_3_2[0]);
  assert(the_pair.first == circ_arc_overlap_lower_part.source());
  assert(the_pair.first.is_left());
  assign(the_pair, vector_for_intersection_overlap_3_2[1]);
  assert(the_pair.first == circ_arc_overlap_lower_part.target());
  assert(!the_pair.first.is_left());

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
  assert(vector_for_intersection_no_x_monotone_1_1.size() == 2);
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_1_1[0]));
  assert(the_pair.first == Circular_arc_endpoint_2(circ_intersections_2_1,
						   circ_intersections_2_4,
						   true));
  assert(the_pair.second == 1u);
  assert(the_pair.first.is_left());
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_1_1[1]));
  assert(the_pair.first == Circular_arc_endpoint_2(circ_intersections_2_1,
						   circ_intersections_2_4,
						   false));
  assert(the_pair.second == 1u);
  assert(!the_pair.first.is_left());

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
  assert(vector_for_intersection_no_x_monotone_1_2.size() == 1);
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_1_2[0]));
  assert(the_pair.first == Circular_arc_endpoint_2(circ_intersections_2_1,
						   circ_intersections_2_5,
						   true));
  assert(the_pair.second == 1u);
  assert(the_pair.first.is_left());


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
  assert(vector_for_intersection_no_x_monotone_1_3.size() == 1);
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_1_3[0]));
  assert(the_pair.first == Circular_arc_endpoint_2(circ_intersections_2_1,
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
  assert(vector_for_intersection_no_x_monotone_1_4.size() == 1);
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_1_4[0]));
  assert(the_pair.first == Circular_arc_endpoint_2(circ_intersections_2_1,
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
  assert(vector_for_intersection_no_x_monotone_1_5.size() == 1);
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_1_5[0]));
  assert(the_pair.first == Circular_arc_endpoint_2(circ_intersections_2_1,
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
  assert(vector_for_intersection_no_x_monotone_1_6.size() == 1);
  assert(assign(the_pair, vector_for_intersection_no_x_monotone_1_6[0]));
  assert(the_pair.first == Circular_arc_endpoint_2(circ_intersections_2_1,
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
  std::cout << vector_for_intersection_no_x_monotone_2_6.size() << std::endl;
  assert(vector_for_intersection_no_x_monotone_2_6.size() == 2);
  assert(assign(circ_arc_in_overlap,
		vector_for_intersection_no_x_monotone_2_6[0]));
  assert(circ_arc_in_overlap.source() == circ_arc_no_x_monotone_7.source());
  std::cout << "source4 : " << std::endl
	    << circ_arc_no_x_monotone_4.source() << std::endl
	    << "target4 : " << std::endl
	    << circ_arc_no_x_monotone_4.target() << std::endl
	    << "source7 : " << std::endl
	    << circ_arc_no_x_monotone_7.source() << std::endl
	    << "target7 : " << std::endl
	    << circ_arc_no_x_monotone_7.target() << std::endl;
  std::cout << "res source : " << std::endl
	    <<circ_arc_in_overlap.source() << std::endl 
	    << "res target : " << std::endl
	    << circ_arc_in_overlap.target() << std::endl;
  assert(circ_arc_in_overlap.target() == circ_arc_no_x_monotone_4.target());
  assert(circ_arc_in_overlap.is_x_monotone());
  assert(assign(circ_arc_in_overlap,
		vector_for_intersection_no_x_monotone_2_6[1]));
  std::cout << "res source : " << std::endl
	    <<circ_arc_in_overlap.source() << std::endl 
	    << "res target : " << std::endl
	    << circ_arc_in_overlap.target() << std::endl;
  assert(circ_arc_in_overlap.source() == circ_arc_no_x_monotone_4.source());
  assert(circ_arc_in_overlap.target() == circ_arc_no_x_monotone_7.target());
  

  
  //Make_x_monotone_2 with a full circle
  Make_x_monotone_2 theMake_x_monotone = ck.make_x_monotone_2_object();
  int x = theRandom.get_int(random_min,random_max);
  int y = theRandom.get_int(random_min,random_max);
  int r = theRandom.get_int(1,random_max);
  Point_2 center_circ_monotone(x,y);
  Circle_2 circ_monotone(center_circ_monotone, r*r);
  Circular_arc_2 theCircular_arc_2_full(circ_monotone);
  std::vector< CGAL::Object > outputIterator1;
  theMake_x_monotone(theCircular_arc_2_full,
		     std::back_inserter(outputIterator1));
  std::cout << std::endl;
  std::cout << "x_monotone full circle : " 
	    << circ_monotone << std::endl;
  Circular_arc_2 circular_arc_2_full;
  for(std::size_t i = 0; i < outputIterator1.size(); i++){
    assign(circular_arc_2_full,  outputIterator1[i]);
    std::cout << "Circular_arc_2_" << i << " : " 
	      << circular_arc_2_full << std::endl;
    std::cout << "Circular_arc_2_" << i << "source : " 
	      << circular_arc_2_full.source() << std::endl;
    std::cout << "Circular_arc_2_" << i << "target : " 
	      << circular_arc_2_full.target() << std::endl;
    assert(circular_arc_2_full.is_x_monotone());
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
  Circular_arc_endpoint_2 circ1_arc_end_p1(circ1, circ1_low_right, true);
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
  Circular_arc_endpoint_2 circ1_arc_end_p2(circ1_low_left, circ1_low_right, true);
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



}
