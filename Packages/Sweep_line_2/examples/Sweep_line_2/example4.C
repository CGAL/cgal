// examples/Sweep_line/example6.C
// ------------------------------

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
  std::vector<Curve_2>  segments;
  Sweep_line sl;

  std::cout << "Demonstrating Sweep_line_2::do_curves_intersect " 
	    << std::endl;

  // case 1 
  Curve_2 c1(Point_2(10,1), Point_2(20,1));
  Curve_2 c2(Point_2(10, -4), Point_2(20,6));
  segments.push_back(c1);
  segments.push_back(c2);

  bool b = sl.do_curves_intersect(segments.begin(),segments.end());

  std::cout << "Curves: " << std::endl
	    << segments[0] << std::endl
	    << segments[1] << std::endl;

  if (b)
    std::cout << "Curves intersect"<<std::endl<<std::endl;
  else
    std::cout << "Curves do NOT intersect"<<std::endl<<std::endl;

  // case 2  
  segments.clear();
  Curve_2 c3(Point_2(10,1), Point_2(20,1));
  Curve_2 c4(Point_2(16, 2), Point_2(20,6));
  segments.push_back(c3);
  segments.push_back(c4);
  b = sl.do_curves_intersect(segments.begin(), segments.end());

  std::cout << "Curves: " << std::endl
	    << segments[0] << std::endl
	    << segments[1] << std::endl;

  if (b)
    std::cout << "Curves intersect"<<std::endl<<std::endl;
  else
    std::cout << "Curves do NOT intersect"<<std::endl<<std::endl;
  
  return 0;
}





