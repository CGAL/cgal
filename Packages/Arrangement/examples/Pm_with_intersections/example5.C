// file: examples/Pm_with_intersections/example1.C

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Arr_segment_cached_traits_2.h>
#include <CGAL/Arr_curve_origin_traits_2.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_cached_traits_2<Kernel>       Base_traits;
typedef CGAL::Arr_curve_origin_traits_2<Base_traits>    Traits;
typedef Traits::Point_2                                 Point_2;
typedef Traits::Curve_2                                 Curve_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef CGAL::Pm_default_dcel<Traits>                   Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                 Pm;
typedef CGAL::Planar_map_with_intersections_2<Pm>       Pmwx;

int main() {
  
  Pmwx pm;
  Curve_2 cv1(Point_2(0, 0), Point_2(1, 1));
  Curve_2 cv2(Point_2(0, 1), Point_2(1, 0)); 

  //insertion of the curves
  std::cout << "Inserting the segments:" << std::endl << cv1 << std::endl;
  pm.insert(cv1);
  std::cout << cv2 << std::endl << std::endl;
  pm.insert(cv2);
  
  //traversal of the curves
  std::cout << "Edges of the planar map:" << std::endl;

  Pmwx::Halfedge_iterator eit;
  for (eit = pm.halfedges_begin(); eit != pm.halfedges_end(); ++eit, ++eit) {
    const X_monotone_curve_2 & curve = eit->curve();
    std::cout << "Curve: " << curve << std::endl;
    const X_monotone_curve_2 * xcv = &curve;
    bool origin_x_monotone = curve.is_origin_x_monotone();
    while (origin_x_monotone) {
      xcv = xcv->get_origin_x_monotone_curve();
      std::cout << "  X monotone curve: " << *xcv << std::endl;
      origin_x_monotone = xcv->is_origin_x_monotone();
    }
    
    const Curve_2 * cv = xcv->get_origin_curve();
    std::cout << "  Curve: " << *cv << std::endl;
  }

  return 0;
}
