// file: examples/Sweep_line_2/example8.C

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_cached_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Sweep_line_2.h>
#include <iostream>
#include <vector>

typedef int Data;

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_cached_traits_2<Kernel>       Seg_traits;
typedef CGAL::Arr_curve_data_traits_2<Seg_traits,Data>  Traits;

typedef Traits::Point_2                                 Point_2;
typedef Traits::Curve_2                                 Curve_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;

typedef std::list<Curve_2>                              CurveList;
typedef CurveList::iterator                             CurveListIter;
typedef CGAL::Sweep_line_2<CurveListIter, Traits>       Sweep_line;

typedef Seg_traits::Curve_2                             Segment_2;

int main()
{
  CurveList  segments;
  Segment_2 s1(Point_2(10,1), Point_2(20,1));
  Segment_2 s2(Point_2(10,-4), Point_2(20,6));
  Segment_2 s3(Point_2(12,-4), Point_2(12,3));
  Segment_2 s4(Point_2(20,6), Point_2(20,1));
  Curve_2 c1(s1, 1);
  Curve_2 c2(s2, 2);
  Curve_2 c3(s3, 3);
  Curve_2 c4(s4, 4);

  segments.push_back(c1);
  segments.push_back(c2);
  segments.push_back(c3);
  segments.push_back(c4);

  // Use a sweep to create the sub curves  
  Traits traits;
  std::list<X_monotone_curve_2> subcurves;
  Sweep_line sl(&traits);
  sl.get_subcurves(segments.begin(), 
		   segments.end(), 
		   std::back_inserter(subcurves), true);
  
  // Write output
  std::cout << std::endl << "Demonstrating Sweep_line_2::get_subcurves " 
	    << std::endl << std::endl << "Curves: " << std::endl
	    << c1 << std::endl << c2 << std::endl
	    << c3 << std::endl << c4 << std::endl << std::endl;
  std::cout <<"Number of sub segments: " << subcurves.size()
            << std::endl<< std::endl;
  for (std::list<X_monotone_curve_2>::iterator scv_iter = subcurves.begin(); 
       scv_iter != subcurves.end(); scv_iter++)
    std::cout << *scv_iter << " Data: " << (*scv_iter).get_data() << std::endl;
  return 0;
}
