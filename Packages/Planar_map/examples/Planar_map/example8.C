// file: examples/Planar_map/example8.C

#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <iostream>

typedef CGAL::Quotient<CGAL::MP_Float>    Number_type;
typedef CGAL::Cartesian<Number_type>      Kernel;
typedef CGAL::Pm_segment_traits_2<Kernel> Traits;
typedef Traits::Point_2                   Point_2;
typedef Traits::X_monotone_curve_2        X_monotone_curve_2;
typedef CGAL::Pm_default_dcel<Traits>     Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>   Planar_map;
typedef Planar_map::Halfedge_handle       Halfedge_handle;

int main()
{
  Planar_map pm;
  X_monotone_curve_2 cv1(Point_2(1.0, 0.0), Point_2(3.0, 2.0));
  X_monotone_curve_2 cv2(Point_2(4.0, -1.0), Point_2(3.0, -2.0));
  X_monotone_curve_2 cv3(Point_2(4.0, -1.0), Point_2(1.0, 0.0));
  X_monotone_curve_2 cv4(Point_2(1.0, 0.0), Point_2(4.0, 1.0));
  X_monotone_curve_2 cv5(Point_2(3.0, 2.0), Point_2(4.0, 1.0));
  X_monotone_curve_2 cv6(Point_2(6.0, 0.0), Point_2(4.0, -1.0));
  X_monotone_curve_2 cv7(Point_2(4.0, 1.0), Point_2(6.0, 0.0));

  Halfedge_handle h1 = pm.insert_in_face_interior(cv1, pm.unbounded_face());
  Halfedge_handle h2 = pm.insert_in_face_interior(cv2, pm.unbounded_face());
  Halfedge_handle h3 = pm.insert_at_vertices(cv3, h2->twin(), h1->twin());
  Halfedge_handle h4 = pm.insert_from_vertex(cv4, h1->twin());
  Halfedge_handle h5 = pm.insert_at_vertices(cv5, h1, h4);
  Halfedge_handle h6 = pm.insert_from_vertex(cv6, h3->twin());
  Halfedge_handle h7 = pm.insert_at_vertices(cv7, h5, h6);
  (void) h7;
  
  std::cout << "Faces = " << pm.number_of_faces() << "\n";
  std::cout << "Vertices = " << pm.number_of_vertices() << "\n";
  std::cout << "Halfedges = " << pm.number_of_halfedges() << "\n";
}
