// Regression test for the Container_ template parameter of Partition_traits_2.
// Verifies that a user-specified container type (std::vector) can be used in
// place of the default std::list for the output Polygon_2 vertex storage.
// Fixes: https://github.com/CGAL/cgal/issues/7545

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <list>
#include <vector>
#include <cassert>

typedef CGAL::Simple_cartesian<double>  K;
typedef K::Point_2                      Point_2;

// Traits that stores polygon vertices in a std::vector instead of std::list.
typedef CGAL::Partition_traits_2<K,
          CGAL::Identity_property_map<Point_2>,
          std::vector<Point_2>>          VecTraits;
typedef VecTraits::Polygon_2            VecPolygon_2;
typedef std::list<VecPolygon_2>         VecPolygon_list;

int main()
{
  // Build a simple non-convex CCW input polygon (same as the standard
  // Partition_2 test suite polygon — proven to exercise all partitioners).
  VecPolygon_2 polygon;
  polygon.push_back(Point_2(227,423));
  polygon.push_back(Point_2(123,364));
  polygon.push_back(Point_2(129,254));
  polygon.push_back(Point_2(230,285));
  polygon.push_back(Point_2(231,128));
  polygon.push_back(Point_2(387,205));
  polygon.push_back(Point_2(417,331));
  polygon.push_back(Point_2(319,225));
  polygon.push_back(Point_2(268,293));
  polygon.push_back(Point_2(367,399));
  polygon.push_back(Point_2(298,418));
  polygon.push_back(Point_2(196,326));

  // Partition using y_monotone — output polygons use std::vector vertices.
  VecPolygon_list partition_polys;
  CGAL::y_monotone_partition_2(polygon.vertices_begin(),
                               polygon.vertices_end(),
                               std::back_inserter(partition_polys),
                               VecTraits());

  assert(!partition_polys.empty());

  // Verify the vertex container is indeed std::vector<Point_2>.
  static_assert(
    std::is_same<VecTraits::Container, std::vector<Point_2>>::value,
    "Container typedef must alias the Container_ template parameter");

  // Partition using approx_convex.
  partition_polys.clear();
  CGAL::approx_convex_partition_2(polygon.vertices_begin(),
                                  polygon.vertices_end(),
                                  std::back_inserter(partition_polys),
                                  VecTraits());
  assert(!partition_polys.empty());

  return 0;
}
