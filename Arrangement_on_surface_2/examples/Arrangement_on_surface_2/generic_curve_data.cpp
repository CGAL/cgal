//! \file examples/Arrangement_on_surface_2/generic_curve_data.cpp
// Associating a name attribute with segments using the generic curve-data
// traits.

#include <CGAL/basic.h>
#include <CGAL/Arr_curve_data_traits_2.h>

#include "arr_polylines.h"

typedef std::string Name;               // The name-field type.

struct Merge_names {
  Name operator() (const Name& s1, const Name& s2) const
  { return (s1 + " " + s2); }
};

typedef CGAL::Arr_curve_data_traits_2<Traits, Name, Merge_names>
                                                    Ex_traits;
typedef Ex_traits::Curve_2                          Ex_polyline;
typedef Ex_traits::X_monotone_curve_2               Ex_x_monotone_polyline;
typedef CGAL::Arrangement_2<Ex_traits>              Ex_arrangement;

int main() {
  // Construct an arrangement of four polylines named A--D.
  Ex_traits traits;
  Ex_arrangement arr(&traits);
  auto ctr_curve = traits.construct_curve_2_object();

  Point pts1[5] =
    {Point(0, 0), Point(2, 4), Point(3, 3), Point(4, 4), Point(6, 0)};
  insert(arr, Ex_polyline(ctr_curve(pts1, pts1 + 5), "A"));

  Point pts2[3] = {Point(1, 5), Point(3, 3), Point(5, 5)};
  insert(arr, Ex_polyline(ctr_curve(pts2, pts2 + 3), "B"));

  Point pts3[4] = {Point(1, 0), Point(2, 2), Point(4, 2), Point(5, 0)};
  insert(arr, Ex_polyline(ctr_curve(pts3, pts3 + 4), "C"));

  Point pts4[2] = {Point(0, 2), Point(6, 2)};
  insert(arr, Ex_polyline(ctr_curve(pts4, pts4 + 2), "D"));

  // Print all edges that correspond to an overlapping polyline.
  std::cout << "The overlapping subcurves:\n";
  for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    if (eit->curve().data().length() > 1) {
      std::cout << "  [" << eit->curve() << "]  "
                << "named: " << eit->curve().data() << std::endl;

      // Modify the curve associated with the edge.
      arr.modify_edge(eit, Ex_x_monotone_polyline(eit->curve(), "overlap"));
    }
  }
  return 0;
}
