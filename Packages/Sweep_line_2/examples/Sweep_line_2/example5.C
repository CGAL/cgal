// examples/Sweep_line/example5.C
// ------------------------------

#if defined(CGAL_CFG_NO_LONGNAME_PROBLEM) || defined(_MSC_VER)
#define Quotient                        Qt
#define Cartesian                       Cn
#define Arr_segment_exact_traits        AST
#define allocator                       All
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
typedef Traits::X_curve_2                       X_curve_2;

typedef std::pair<Point_2, std::list<X_curve_2> > Pair_point_curves;

typedef std::vector<Curve_2>                    CurveList;
typedef CurveList::iterator                     CurveListIter;
typedef CGAL::Sweep_line_2<CurveListIter, Traits> Sweep_line;

int main()
{
  std::vector<Curve_2>  segments;
  std::vector<Pair_point_curves> intersecting_curves;
  Traits traits;
  Sweep_line sl(&traits);
  
  Curve_2 c1(Point_2(10,1), Point_2(20,1));
  Curve_2 c2(Point_2(10, -4), Point_2(20,6));
  segments.push_back(c1);
  segments.push_back(c2);

  sl.get_intersecting_curves(segments.begin(), segments.end(), 
			     std::back_inserter(intersecting_curves), false);

  std::cout << "Demonstrating Sweep_line_2::get_intersecting curves " 
	    << std::endl
	    << "Input curves: " << std::endl
	    << c1 << std::endl
	    << c2 << std::endl << std::endl;

  std::cout << "Output: " << std::endl;

  for (std::vector<Pair_point_curves>::iterator 
         p_iter = intersecting_curves.begin();
       p_iter != intersecting_curves.end(); ++p_iter)
  {
    std::cout<<"The curves"<< std::endl;
    for (std::list<X_curve_2>::iterator cv_iter = (p_iter->second).begin();
	 cv_iter != (p_iter->second).end(); ++cv_iter)
      std::cout<< *cv_iter << std::endl;
    std::cout<<"intersect at "<< p_iter->first << std::endl << std::endl;
  }
  
  return 0;
}
