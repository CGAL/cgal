// file: examples/Pm_with_intersections/example2.C

#include "short_names.h"

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
typedef Pmwx::Pmwx_change_notification Pmwx_change_notification;

class My_notification : public Pmwx_change_notification 
{
public:

  My_notification()
  {i = 0;}

  void add_edge(const  Traits::X_monotone_curve_2 &,
                Planar_map::Halfedge_handle, 
                bool /* left_to_right */, bool overlap = false)
  {
    (void) overlap;
    std::cout << "add_edge" << std::endl;
    i++;
  }

  void split_edge(Planar_map::Halfedge_handle /* orig_edge */, 
                  Planar_map::Halfedge_handle /* new_edge */,
                  const Traits::X_monotone_curve_2 &,
                  const Traits::X_monotone_curve_2 &)
  {
    std::cout << "split_edge" << std::endl;
    i++;
  }

  void split_face(Planar_map::Face_handle /* orig_face */, 
                  Planar_map::Face_handle /* new_face */)
  {
    std::cout << "split_face" << std::endl;
  }

  void add_hole(Planar_map::Face_handle /* in_face */, 
                Planar_map::Halfedge_handle /* new_hole */)
  {
    std::cout << "add_hole" << std::endl;
  }

  int i;
};

int main() {
  
  Pmwx pm;
  My_notification notif;

  //insertion of the curves
  X_monotone_curve_2 c1(Point_2(0, 1), Point_2(1, 0));
  X_monotone_curve_2 c2(Point_2(0, 0), Point_2(1, 1));
  X_monotone_curve_2 c3(Point_2(0, 1), Point_2(1, 1));

  std::cout << "inserting " << c1 << std::endl;
  pm.insert(c1, &notif);
  std::cout << "inserting " << c2 << std::endl;
  pm.insert(c2, &notif);
  std::cout << "inserting " << c3 << std::endl;
  pm.insert(c3, &notif);

  std::cout << "Total number of edges " << notif.i << std::endl;

  return 0;
}
