// examples/Pm_with_intersections/example1
// ---------------------------------------

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>

typedef CGAL::Quotient<CGAL::MP_Float>                      NT;
typedef CGAL::Cartesian<NT>                                 Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                  Traits;
typedef Traits::Point_2                                     Point_2;
typedef Traits::X_monotone_curve_2                          X_monotone_curve_2;
typedef CGAL::Pm_default_dcel<Traits>                       Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                     Planar_map_2;
typedef CGAL::Planar_map_with_intersections_2<Planar_map_2> Pmwx;

int main() {
  
  Pmwx pm;
  X_monotone_curve_2 cv1(Point_2(0, 0), Point_2(1, 1));
  X_monotone_curve_2 cv2(Point_2(0, 1), Point_2(1, 0)); 

  //insertion of the curves
  std::cout << "Inserting the segments:" << std::endl << cv1 << std::endl;
  pm.insert(cv1);
  std::cout << cv2 << std::endl << std::endl;
  pm.insert(cv2);
  
  //traversal of the curves
  std::cout << "Edges of the planar map:" << std::endl;

  Pmwx::Halfedge_iterator eit;
  for (eit = pm.halfedges_begin(); eit != pm.halfedges_end(); ++eit, ++eit) {
    std::cout << eit->source()->point()
              << " --- " << eit->target()->point() << std::endl;
  }

  return 0;
}
