// file: examples/Sweep_line_2/example3.C

#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Sweep_line_2.h> 

#include <iostream>
#include <vector>

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;

typedef Traits::Point_2                                 Point_2;
typedef Traits::Curve_2                                 Curve_2;

typedef std::vector<Curve_2>                            CurveList;
typedef CurveList::iterator                             CurveListIter;
typedef CGAL::Sweep_line_2<CurveListIter, Traits>       Sweep_line;

int main()
{
  CurveList segments;
  
  Curve_2 c1(Point_2(10,1), Point_2(20,1));
  Curve_2 c2(Point_2(10, -4), Point_2(20,6));

  segments.push_back(c1);
  segments.push_back(c2);

  std::vector<Point_2> points;
  Sweep_line sl;
  sl.get_intersection_points(segments.begin(),
			     segments.end(), 
			     std::back_inserter(points), false);

  // Write results...

  std::cout << " Demonstrating Sweep_line_2::get_intersection_points " 
	    << std::endl
	    << " Curves: " << std::endl
	    << " " << c1 << std::endl
	    << " " << c2 << std::endl << std::endl;
  
  std::cout << " Intersection point (not including end points): " << std::endl;
  
  for (std::vector<Point_2>::iterator p_iter = points.begin();
       p_iter != points.end(); ++p_iter)
    std::cout<< " " << *p_iter << std::endl;

  std::cout << std::endl << std::endl;

  return 0;
}
