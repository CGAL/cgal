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
void _test_circle_predicat(CK ck)
{
  typedef typename CK::FT                      FT;
  typedef typename CK::Circle_2                Circle_2;
  typedef typename CK::Circular_arc_2          Circular_arc_2;
  typedef typename CK::Point_2                 Point_2;
  typedef typename CK::Line_2                  Line_2;
  typedef typename CK::Compare_x_2             Compare_x_2;
  typedef typename CK::Compare_y_2             Compare_y_2;
  typedef typename CK::Compare_xy_2            Compare_xy_2;
  typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;
  typedef typename CK::Compare_y_to_right_2    Compare_y_to_right_2;
  typedef typename CK::Compare_y_at_x_2        Compare_y_at_x_2;
  typedef typename CK::Equal_2                 Equal_2;
  typedef typename CK::In_range_2              In_range_2;
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
  Circular_arc_endpoint_2 circ1_arc_end_p1(circ1, circ1_low_right, true);
  Circular_arc_endpoint_2 circ1_arc_end_p2(circ1, circ1_low_right, false);

  assert(theCompare_x_2(circ1_arc_end_p1,circ1_arc_end_p2 )== CGAL::SMALLER);
  assert(theCompare_y_2(circ1_arc_end_p1,circ1_arc_end_p2) == CGAL::SMALLER);
  assert(theCompare_xy_2(circ1_arc_end_p1,circ1_arc_end_p2) == CGAL::SMALLER);

  assert(theCompare_x_2(circ1_arc_end_p2,circ1_arc_end_p2) == CGAL::EQUAL);
  assert(theCompare_y_2(circ1_arc_end_p2,circ1_arc_end_p2) == CGAL::EQUAL);
  assert(theCompare_xy_2(circ1_arc_end_p2,circ1_arc_end_p2) == CGAL::EQUAL);

  assert(theCompare_x_2(circ1_arc_end_p2,circ1_arc_end_p1) == CGAL::LARGER);
  assert(theCompare_y_2(circ1_arc_end_p2,circ1_arc_end_p1) == CGAL::LARGER);
  assert(theCompare_xy_2(circ1_arc_end_p2,circ1_arc_end_p1) == CGAL::LARGER);

  //We create a circle in top of the circle1
  Point_2 center1_high(center1_x, center1_y + circ1_r);
  Circle_2 circ1_high(center1_high, circ1_r * circ1_r);

  //p3 is in the quarter superior left
  Circular_arc_endpoint_2 circ1_arc_end_p3(circ1, circ1_high, true);
  
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
  Circular_arc_endpoint_2 circ1_arc_end_p4(circ1, circ1_left, true);

  //Comparison between the superior arc and the point p4 
  theComparison_result_y_at_x_2 = 
    theCompare_y_at_x_2(circ1_arc_end_p4,circ1_arc_high);
  assert(theComparison_result_y_at_x_2 == CGAL::EQUAL);

  //Comparison between the inferior arc and the point p4 
  theComparison_result_y_at_x_2 = 
    theCompare_y_at_x_2(circ1_arc_end_p4,circ1_arc_low);
  assert(theComparison_result_y_at_x_2 == CGAL::EQUAL);
  

  
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
  assert(theComparison_result_y_to_right_2 == CGAL::LARGER);
  
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


    
  //P3 is the leftest point of circle1
  //P3, P4 and P5 are the same
  Circular_arc_endpoint_2 circ1_arc_end_p5(circ1, circ1_left, true);
  Circular_arc_endpoint_2 circ1_arc_end_p6(circ1, circ1_left, false);
  Circular_arc_endpoint_2 circ1_arc_end_p7(circ1_left, circ1, true);
  
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
  Circular_arc_endpoint_2 circ2_arc_end_p1(circ2, 
					   circ2_left, 
					   false);
  //We create a circle lefter than circ1_left but it cuts it
  Point_2 center1_very_left(center1_x - (3 * circ1_r), center1_y);
  Circle_2 circ1_very_left(center1_very_left, circ1_r * circ1_r);
  //Point on top
  Circular_arc_endpoint_2 circ1_left_arc_end_p1(circ1_left,
						circ1_very_left,
						false);
  std::cout << "In range" << std::endl;
  In_range_2 theIn_range_2 = ck.in_range_2_object();
  assert(theIn_range_2(circ1_arc_low, circ1_arc_end_p7));
  assert(theIn_range_2(circ1_arc_high, circ1_arc_end_p7));
  assert(theIn_range_2(circ1_arc_high, circ2_arc_end_p1));
  assert(theIn_range_2(circ1_arc_low, circ2_arc_end_p1));
  assert(!theIn_range_2(circ1_arc_low, circ1_left_arc_end_p1));
  std::cout << std::endl;

  
  //////////////////////////DO OVERLAP/////////////////////////

  std::cout << "Do_overlap_2_object" << std::endl;
  Do_overlap_2 theDo_overlap_2 = ck.do_overlap_2_object();
  Point_2 point_2_high(center1_x, center1_y + circ1_r);
  Line_2 theLine_2_vertical(center1,point_2_high);
  Circular_arc_2 circ1_arc_low_right(circ1,
				     theLine_2_vertical,true,
				     theLine_2_horizontal, false);
  assert(!theDo_overlap_2(circ1_arc_high, circ1_arc_low_right));
  assert(theDo_overlap_2(circ1_arc_low, circ1_arc_low_right));
  assert(theDo_overlap_2(circ1_arc_low_right, circ1_arc_low));
  assert(!theDo_overlap_2(circ1_arc_low, circ1_arc_high));
  std::cout << std::endl;


  
}
