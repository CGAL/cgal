// examples/Pm_with_intersections/example4
// ---------------------------------------

#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>
#include <algorithm>
#include <iostream>

typedef CGAL::Quotient<int>                             NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
typedef Traits::Point_2                                 Point;
typedef Traits::X_monotone_curve_2                      Curve;
typedef CGAL::Pm_default_dcel<Traits>                   Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                 PM;
typedef CGAL::Planar_map_with_intersections_2<PM>       Pmwx;
typedef std::vector<Curve>                              CurveContainer;

int main()
{
  // Create an instance of a Planar_map_with_intersections:
  Pmwx pm;
  CurveContainer cv;

  cv.push_back(Curve(Point(0, 1), Point(1, 0)));
  cv.push_back(Curve(Point(0, 0), Point(1, 1)));
  cv.push_back(Curve(Point(0, 1), Point(1, 1)));
  
  std::cout << "The curves of the map :" << std::endl;
  std::copy(cv.begin(), cv.end(), 
	    std::ostream_iterator<Curve>(std::cout, "\n"));
  std::cout << std::endl;

  // Insert the curves into the Planar_map:
  std::cout << "Inserting the curves to the map ... \n";
  pm.insert(cv.begin(), cv.end());
  
  // insert another segment...
  cv.clear();
  cv.push_back(Curve(Point(0, 0), Point(1, 2)));
  pm.insert(cv.begin(), cv.end());
  
  std::cout << ((pm.is_valid()) ? "map valid!" : "map invalid!") << std::endl
            << std::endl;

  std::cout << "Edges of the planar map:" << std::endl;
  Pmwx::Halfedge_iterator eit;
  for (eit = pm.halfedges_begin(); eit != pm.halfedges_end(); ++eit, ++eit) {
    std::cout << eit->source()->point() << " --- " << eit->target()->point()
              << std::endl;
  }

  return 0;
}
