#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <fstream>

// Making sure test doesn't fail if LEDA is not installed
#if ! defined(CGAL_USE_LEDA) && \
      (CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || \
       CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS)

int main(int argc, char* argv[])
{
  std::cout << "A try to run test with LEDA traits but LEDA is not installed.";
  std::cout << std::endl;
  std::cout << "Test is not performed.";
  std::cout << std::endl;

  return 0;
}
#else

#include <CGAL/Circle_2.h>
#include <CGAL/Arr_circles_real_traits.h>
#include <CGAL/leda_real.h> 

#include "include/Circles_traits_test_base.h"

typedef leda_real                                     NT;
typedef CGAL::Arr_circles_real_traits<NT>             Traits;
 
typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

class Arr_circles_real_traits_test : public Arr_traits_test<Traits> {

  typedef Arr_traits_test<Traits>::Traits Traits;

  void build_curve_list(std::list<Curve_with_info>& curve_list)
  { 
    // Require:
    CGAL_precondition_msg(curve_list.empty(), \
			  "list is not empty.");
       
    const Point CENTER(NT(0), NT(0));
    const NT    RADIUS = NT(100); 
    const NT    SQUARED_RADIUS = RADIUS * RADIUS;
    const NT    sin45 = CGAL::sqrt( NT(2) ) / NT(2);
    const NT    cos45 = sin45;
     
    // Defining four points p1, p2, p3, p4 in quadrant I, II, III, IV
    // resepectively.
    Point p1(  RADIUS *  cos45, RADIUS *  sin45);
    Point p2(  RADIUS * -cos45, RADIUS *  sin45);
    Point p3(  RADIUS * -cos45, RADIUS * -sin45);
    Point p4(  RADIUS *  cos45, RADIUS * -sin45);
    
    CGAL::Circle_2<Curve::R> my_circle(CENTER, SQUARED_RADIUS, 
				       CGAL::COUNTERCLOCKWISE);

    // Check that points are indeed on boundry of circle
    CGAL_assertion( my_circle.has_on_boundary(p1) );
    CGAL_assertion( my_circle.has_on_boundary(p2) );
    CGAL_assertion( my_circle.has_on_boundary(p3) );
    CGAL_assertion( my_circle.has_on_boundary(p4) );

    Curve_with_info cv;

    // Defining curves and inserting into container
    // Curves are oriented counterclockwise.
   
    // X-monotone curves
    cv = Curve_with_info(Curve(my_circle, p1, p2),
			 true,
			 1, 
			 "arc of quadrants I - II");
    curve_list.push_back(cv);

    cv = Curve_with_info(Curve(my_circle, p3, p4),
			 true,
			 1,
			 "arc of quadrants III - IV");
    curve_list.push_back(cv);

  // Non x-monotone circular curves with 2 x-monotone parts
    cv = Curve_with_info(Curve(my_circle),
			 false,
			 2,
			 "a whole circle");
    curve_list.push_back(cv);

    cv = Curve_with_info(Curve(my_circle, p4, p1),
			 false,
			 2,
			 "arc of quadrants IV - I");
    curve_list.push_back(cv);

    cv = Curve_with_info(Curve(my_circle, p3, p1),
			 false,
			 2,
			 "arc of quadrants III - IV - I");
    curve_list.push_back(cv);


    cv = Curve_with_info(Curve(my_circle, p4, p2),
			 false,
			 2,
			 "arc of quadrants IV - I - II");
    curve_list.push_back(cv);

    cv = Curve_with_info(Curve(my_circle, p3, p2),
			 false,  
			 2,
			 "arc of quadrants III - IV - I - II");
    curve_list.push_back(cv);
 
    cv = Curve_with_info(Curve(my_circle, p4, p2),
			 false,
			 2,
			 "arc of quadrants IV - I - II");
    curve_list.push_back(cv);
 
    cv = Curve_with_info(Curve(my_circle, p2, p3),
 			 false,
			 2,
			 "arc of quadrants II - III");
    curve_list.push_back(cv);

    cv = Curve_with_info(Curve(my_circle, p1, p3),
			 false,
			 2,
			 "arc of quadrants I - II - III");
    curve_list.push_back(cv);


    cv = Curve_with_info(Curve(my_circle, p2, p4),
			 false,
			 2,
			 "arc of quadrants II - III - IV");
    curve_list.push_back(cv);

    cv = Curve_with_info(Curve(my_circle, p1, p4),
			 false,
			 2,
			 "arc of quadrants I - II - III - IV");
    curve_list.push_back(cv);

  // Non x-monotone circular curves with 3 x-monotone parts
    cv = Curve_with_info(Curve(my_circle, p2, p1),
			 false,
			 3,
			 "arc of quadrants II - III - IV - I");
    curve_list.push_back(cv);

    cv = Curve_with_info(Curve(my_circle, p4, p3),
			 false,
			 3,
			 "arc of quadrants II - III");
    curve_list.push_back(cv);
  }

}; // Arr_traits_test

int main()  
{
  Arr_circles_real_traits_test test;
 
  if ( test.start() )    
    return 0; 
  else
    return 1;
}
 
#endif // ! defined(CGAL_USE_LEDA) ...
