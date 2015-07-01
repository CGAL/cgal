//! \file examples/Arrangement_on_surface_2/polylines.cpp
// Constructing an arrangement of polylines.

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_directional_non_caching_segment_basic_traits_2.h>
#include <CGAL/Arr_polycurve_basic_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <vector>
#include <list>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Arr_directional_non_caching_segment_basic_traits_2<Kernel>
                                                          Subcurve_traits_2;
typedef CGAL::Arr_polycurve_basic_traits_2<Subcurve_traits_2>
                                                          Geom_traits_2;
typedef Geom_traits_2::Point_2                            Point_2;
typedef Subcurve_traits_2::X_monotone_curve_2             X_monotone_subcurve_2;
typedef Geom_traits_2::X_monotone_curve_2                 X_monotone_curve_2;
typedef CGAL::Arrangement_2<Geom_traits_2>                Arrangement_2;

int main()
{
  Geom_traits_2 traits;
  Arrangement_2 arr(&traits);

  Geom_traits_2::Construct_x_monotone_curve_2 ctr =
    traits.construct_x_monotone_curve_2_object();

  std::vector<X_monotone_subcurve_2> segs1;
  segs1.push_back(X_monotone_subcurve_2(Point_2(0, 0), Point_2(1, 1)));
  segs1.push_back(X_monotone_subcurve_2(Point_2(1, 1), Point_2(2, 2)));
  segs1.push_back(X_monotone_subcurve_2(Point_2(2, 2), Point_2(3, 1)));
  segs1.push_back(X_monotone_subcurve_2(Point_2(3, 1), Point_2(4, 0)));
  X_monotone_curve_2 pc1 = ctr(segs1.begin(), segs1.end());

  std::vector<X_monotone_subcurve_2> segs2;
  segs2.push_back(X_monotone_subcurve_2(Point_2(0, 0), Point_2(1, 1)));
  segs2.push_back(X_monotone_subcurve_2(Point_2(1, 1), Point_2(2, 2)));
  segs2.push_back(X_monotone_subcurve_2(Point_2(2, 2), Point_2(3, 1)));
  segs2.push_back(X_monotone_subcurve_2(Point_2(3, 1), Point_2(4, 0)));
  X_monotone_curve_2 pc2 = ctr(segs2.begin(), segs2.end());

  insert_non_intersecting_curve(arr, pc1);
  insert_non_intersecting_curve(arr, pc2);

  std::cout << "# vertices: " << arr.number_of_vertices() << std::endl;;
  std::cout << "# halfedges: " << arr.number_of_halfedges() << std::endl;;
  std::cout << "# faces: " << arr.number_of_faces() << std::endl;;
  return 0;
}
