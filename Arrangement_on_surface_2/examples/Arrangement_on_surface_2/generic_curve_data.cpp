//! \file examples/Arrangement_on_surface_2/generic_curve_data.cpp
// Associating a name attribute with segments using the generic curve-data
// traits.

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <string>

// Define a functor for concatenating name fields.
typedef std::string   Name;

struct Merge_names
{
  Name operator() (const Name& s1, const Name& s2) const
  {
    return (s1 + " " + s2);
  }
};

typedef CGAL::Cartesian<CGAL::Exact_rational>           Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>   Polyline_traits_2;
typedef Polyline_traits_2::Curve_2                      Polyline_2;
typedef CGAL::Arr_curve_data_traits_2<Polyline_traits_2, Name, Merge_names>
                                                        Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Curve_2                               Curve_2;
typedef Traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;

int main ()
{
  Polyline_traits_2 traits;
  Polyline_traits_2::Construct_curve_2 poly_const =
    traits.construct_curve_2_object();

  // Construct an arrangement of four polylines named A--D.
  Arrangement_2    arr;

  Point_2          points1[5] = {Point_2(0,0), Point_2(2,4), Point_2(3,3),
                                 Point_2(4,4), Point_2(6,0)};
  insert (arr, Curve_2 (poly_const (points1, points1 + 5), "A"));

  Point_2          points2[3] = {Point_2(1,5), Point_2(3,3), Point_2(5,5)};
  insert (arr, Curve_2 (poly_const (points2, points2 + 3), "B"));

  Point_2          points3[4] = {Point_2(1,0), Point_2(2,2),
                                 Point_2(4,2), Point_2(5,0)};
  insert (arr, Curve_2 (poly_const (points3, points3 + 4), "C"));

  Point_2          points4[2] = {Point_2(0,2), Point_2(6,2)};
  insert (arr, Curve_2 (poly_const (points4, points4 + 2), "D"));

  // Print all edges that correspond to an overlapping polyline.
  Arrangement_2::Edge_iterator    eit;

  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    if (eit->curve().data().length() > 1) {
      std::cout << "[" << eit->curve() << "]  "
                << "named: " << eit->curve().data() << std::endl;

      // Rename the curve associated with the edge.
      arr.modify_edge (eit, X_monotone_curve_2 (eit->curve(), "overlap"));
    }
  }
  return 0;
}
