// examples/Sweep_line/example5.C
// ------------------------------

#if defined(CGAL_CFG_NO_LONGNAME_PROBLEM) || defined(_MSC_VER)
#define Quotient                        Qt
#define Cartesian                       Cn
#define Arr_segment_exact_traits        AST
#define allocator                       All
#define Sweep_curves_to_subcurves_utils SCSU
#define Sweep_curves_base_2             SCB
#define Intersection_point_node         IPN
#endif

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Sweep_line_2.h> 

#include <iostream>
#include <vector>

typedef CGAL::Quotient<CGAL::MP_Float>          NT;
typedef CGAL::Cartesian<NT>                     Kernel;
typedef CGAL::Arr_segment_exact_traits<Kernel>  Traits;

typedef Traits::Point_2                         Point_2;
typedef Traits::Curve_2                         Curve_2;

typedef std::vector<Curve_2> CurveList;
typedef CurveList::iterator CurveListIter;
typedef CGAL::Sweep_line_2<CurveListIter, Traits> Sweep_line;

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
