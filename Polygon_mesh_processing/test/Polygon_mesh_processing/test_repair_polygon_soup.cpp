#define CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <algorithm>
#include <deque>
#include <iostream>
#include <fstream>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef K::Point_3                                              Point_3;

typedef CGAL::Surface_mesh<Point_3>                             Mesh;
typedef std::vector<std::size_t>                                CGAL_polygon;

void test_polygon_canonicalization(const bool verbose = false)
{
  std::cout << "test polygon canonicalization... " << std::endl;

  std::vector<Point_3> points;
  points.push_back(Point_3(0,0,0)); // #0
  points.push_back(Point_3(1,0,0)); // #1
  points.push_back(Point_3(0,1,0)); // #2
  points.push_back(Point_3(1,1,0)); // #3
  points.push_back(Point_3(1,1,2)); // #4
  points.push_back(Point_3(1,1,-2)); // #5

  // empty
  CGAL_polygon polygon;
  CGAL_polygon canonical_polygon = PMP::internal::construct_canonical_polygon(points, polygon, K());
  assert(canonical_polygon.empty());

  // 1 point
  polygon.push_back(5);

  canonical_polygon = PMP::internal::construct_canonical_polygon(points, polygon, K());
  assert(canonical_polygon == polygon);

  // 2 points
  polygon.clear();
  polygon.push_back(4); polygon.push_back(3);

  canonical_polygon = PMP::internal::construct_canonical_polygon(points, polygon, K());
  assert(canonical_polygon[0] == 3 && canonical_polygon[1] == 4);
  std::swap(polygon[0], polygon[1]);
  canonical_polygon = PMP::internal::construct_canonical_polygon(points, polygon, K());
  assert(canonical_polygon[0] == 3 && canonical_polygon[1] == 4);

  // 3 points
  polygon.clear();
  polygon.push_back(4); polygon.push_back(1); polygon.push_back(3);

  canonical_polygon = PMP::internal::construct_canonical_polygon(points, polygon, K());
  assert(canonical_polygon[0] == 1 && canonical_polygon[1] == 3 && canonical_polygon[2] == 4);
  std::swap(polygon[0], polygon[2]);
  canonical_polygon = PMP::internal::construct_canonical_polygon(points, polygon, K());
  assert(canonical_polygon[0] == 1 && canonical_polygon[1] == 3 && canonical_polygon[2] == 4);

  // Generic case
  polygon.clear();
  polygon.push_back(0); polygon.push_back(2); polygon.push_back(5); polygon.push_back(4); polygon.push_back(1);

  // Whether reversed or with cyclic permutations, it should always yield the same canonical polygon
  canonical_polygon = PMP::internal::construct_canonical_polygon(points, polygon, K());
  assert(canonical_polygon.size() == 5);
  assert(canonical_polygon[0] == 0);

  // all cyclic permutations
  for(std::size_t i=0, end=polygon.size(); i<end; ++i)
  {
    CGAL_polygon cpol = PMP::internal::construct_canonical_polygon(points, polygon, K());
    if(verbose)
    {
      std::cout << "Input polygon:";
      PMP::internal::print_polygon(std::cout, polygon);
      std::cout << "Canonical polygon:";
      PMP::internal::print_polygon(std::cout, canonical_polygon);
    }
    assert(cpol == canonical_polygon);

    std::rotate(polygon.begin(), polygon.begin()+1, polygon.end());
  }

  // reverse and all cyclic permutations
  std::reverse(polygon.begin(), polygon.end());
  for(std::size_t i=0, end=polygon.size(); i<end; ++i)
  {
    CGAL_polygon cpol = PMP::internal::construct_canonical_polygon(points, polygon, K());
    if(verbose)
    {
      std::cout << "Input polygon:";
      PMP::internal::print_polygon(std::cout, polygon);
      std::cout << "Canonical polygon:";
      PMP::internal::print_polygon(std::cout, canonical_polygon);
    }
    assert(cpol == canonical_polygon);

    std::rotate(polygon.begin(), polygon.begin()+1, polygon.end());
  }
}

void test_merge_duplicate_points(const bool /*verbose*/ = false)
{
  std::cout << "test merge duplicate points... " << std::endl;

  std::vector<Point_3> points;
  std::vector<CGAL_polygon> polygons;

  // empty
  std::size_t res = PMP::merge_duplicate_points_in_polygon_soup(points, polygons);
  assert(res == 0 && points.empty() && polygons.empty());

  points.push_back(Point_3(0,0,0)); // #0
  points.push_back(Point_3(1,2,0)); // #1
  points.push_back(Point_3(1,1,0)); // #2
  points.push_back(Point_3(1,1,0)); // #3 // identical to #2
  points.push_back(Point_3(0,1,0)); // #4
  points.push_back(Point_3(1,1,0)); // #5 // identical to #2
  points.push_back(Point_3(0,0,0)); // #6 // idental to #0

  CGAL_polygon polygon;
  polygon.push_back(0); polygon.push_back(1); polygon.push_back(2);
  polygons.push_back(polygon);

  polygon.clear();
  polygon.push_back(6); polygon.push_back(3); polygon.push_back(1); polygon.push_back(0);
  polygons.push_back(polygon);

  polygon.clear();
  polygon.push_back(5);
  polygons.push_back(polygon);

  res = PMP::merge_duplicate_points_in_polygon_soup(points, polygons, params::geom_traits(K()));
  assert(res == 3 && points.size() == 4 && polygons.size() == 3);

  assert(polygons[0][0] == 0 && polygons[0][1] == 1 && polygons[0][2] == 2);
  assert(polygons[1][0] == 0 && polygons[1][1] == 2 && polygons[1][2] == 1 && polygons[1][3] == 0);

  for(std::size_t i=0, psn=polygons.size(); i<psn; ++i)
  {
    const CGAL_polygon& polygon = polygons[i];
    for(std::size_t j=0, pn=polygon.size(); j<pn; ++j)
    {
      assert(polygon[j] < points.size());
    }
  }
}

void test_merge_duplicate_polygons(const bool /*verbose*/ = false)
{
  std::cout << "test duplicate polygons merging..." << std::endl;

  std::vector<Point_3> points;
  std::vector<CGAL_polygon> polygons;

  points.push_back(Point_3(0,0,0)); // #0
  points.push_back(Point_3(1,0,0)); // #1
  points.push_back(Point_3(0,1,0)); // #2
  points.push_back(Point_3(1,1,0)); // #3
  points.push_back(Point_3(1,1,2)); // #4

  // -------------------------------------------------------
  // empty
  std::vector<std::vector<std::size_t> > all_duplicate_polygons;
  PMP::internal::collect_duplicate_polygons(points, polygons,
                                            std::back_inserter(all_duplicate_polygons),
                                            K(), true /*only equal if same orientation*/);
  assert(polygons.empty() && all_duplicate_polygons.empty());

  std::size_t res = PMP::merge_duplicate_polygons_in_polygon_soup(points, polygons, params::geom_traits(K()));
  assert(res == 0 && polygons.empty());

  // -------------------------------------------------------
  // 1 polygon
  CGAL_polygon polygon;
  polygon.push_back(0); polygon.push_back(1); polygon.push_back(2);
  polygons.push_back(polygon);

  PMP::internal::collect_duplicate_polygons(points, polygons,
                                            std::back_inserter(all_duplicate_polygons),
                                            K(), false /*equal regardless of orientation*/);
  assert(all_duplicate_polygons.empty());

  PMP::internal::collect_duplicate_polygons(points, polygons,
                                            std::back_inserter(all_duplicate_polygons),
                                            K(), true /*only equal if same orientation*/);
  assert(all_duplicate_polygons.empty());

  res = PMP::merge_duplicate_polygons_in_polygon_soup(points, polygons, params::geom_traits(K()));
  assert(res == 0 && polygons.size() == 1);

  // -------------------------------------------------------
  // 2 different polygons
  polygon.clear();
  polygon.push_back(0); polygon.push_back(1); polygon.push_back(3); polygon.push_back(4);
  polygons.push_back(polygon);

  PMP::internal::collect_duplicate_polygons(points, polygons,
                                            std::back_inserter(all_duplicate_polygons),
                                            K(), false /*equal regardless of orientation*/);
  assert(all_duplicate_polygons.empty());

  PMP::internal::collect_duplicate_polygons(points, polygons,
                                            std::back_inserter(all_duplicate_polygons),
                                            K(), true /*only equal if same orientation*/);
  assert(all_duplicate_polygons.empty());

  res = PMP::merge_duplicate_polygons_in_polygon_soup(points, polygons);
  assert(res == 0 && polygons.size() == 2);

  // -------------------------------------------------------
  // Multiple duplicates
    // duplicate, same orientations
  polygon.clear();
  polygon.push_back(2); polygon.push_back(0); polygon.push_back(1);
  polygons.push_back(polygon);

    // duplicate, different orientations
  polygon.clear();
  polygon.push_back(2); polygon.push_back(1); polygon.push_back(0);
  polygons.push_back(polygon);

    // duplicate, different orientations
  polygon.clear();
  polygon.push_back(4); polygon.push_back(3); polygon.push_back(1); polygon.push_back(0);
  polygons.push_back(polygon);

    // same vertices, not a duplicate
  polygon.clear();
  polygon.push_back(3); polygon.push_back(4); polygon.push_back(1); polygon.push_back(0);
  polygons.push_back(polygon);

  PMP::internal::collect_duplicate_polygons(points, polygons,
                                            std::back_inserter(all_duplicate_polygons),
                                            K(), true /*only equal if same orientation*/);
  assert(all_duplicate_polygons.size() == 1); // one duplication
  assert(all_duplicate_polygons[0].size() == 2); // two polygons are equal
  all_duplicate_polygons.clear();

  PMP::internal::collect_duplicate_polygons(points, polygons,
                                            std::back_inserter(all_duplicate_polygons),
                                            K(), false /*equal regardless of orientation*/);
  assert(all_duplicate_polygons.size() == 2);

  // not sure which duplicate is output first
  if(all_duplicate_polygons[0][0] == 0)
    assert(all_duplicate_polygons[0].size() == 3 && all_duplicate_polygons[1].size() == 2);
  else
    assert(all_duplicate_polygons[0].size() == 2 && all_duplicate_polygons[1].size() == 3);

  // Keep one for each duplicate
  std::vector<CGAL_polygon> polygons_copy(polygons);
  res = PMP::merge_duplicate_polygons_in_polygon_soup(points, polygons_copy,
                                                      params::default_values());
  assert(res == 3 && polygons_copy.size() == 3);

  // Remove all duplicates
  polygons_copy = polygons;
  res = PMP::merge_duplicate_polygons_in_polygon_soup(points, polygons_copy,
                                                      params::erase_all_duplicates(true)
                                                             .require_same_orientation(false));
  assert(res == 5 && polygons_copy.size() == 1);

  // Remove all duplicates but different orientations are different polygons
  res = PMP::merge_duplicate_polygons_in_polygon_soup(points, polygons,
                                                      params::erase_all_duplicates(true)
                                                             .require_same_orientation(true));
  assert(res == 2 && polygons.size() == 4);
}

void test_simplify_polygons(const bool /*verbose*/ = false)
{
  std::cout << "test simplify_polygons... " << std::endl;

  std::vector<Point_3> points;
  std::vector<CGAL_polygon> polygons;

  points.push_back(Point_3(0,0,0)); // #0
  points.push_back(Point_3(1,2,0)); // #1
  points.push_back(Point_3(1,0,0)); // #2
  points.push_back(Point_3(1,3,0)); // #3
  points.push_back(Point_3(0,1,0)); // #4
  points.push_back(Point_3(1,1,0)); // #5
  points.push_back(Point_3(0,0,0)); // #6 == #0

  // ------
  CGAL_polygon polygon;
  polygon.push_back(0); polygon.push_back(2); polygon.push_back(4); polygon.push_back(0); polygon.push_back(0);
  polygons.push_back(polygon);

  std::size_t res = PMP::internal::simplify_polygons_in_polygon_soup<K>(points, polygons);
  std::cout << "res: " << res << " / size: " << polygons.back().size() << std::endl;
  assert(res == 1 && polygons.back().size() == 3);

  // ------
  polygon.clear();
  polygon.push_back(0); polygon.push_back(2); polygon.push_back(4);
  polygons.push_back(polygon);

  res = PMP::internal::simplify_polygons_in_polygon_soup<K>(points, polygons);
  std::cout << "res: " << res << " / size: " << polygons.back().size() << std::endl;
  assert(res == 0 && polygons.back().size() == 3);

  // ------
  polygon.clear();
  polygon.push_back(0); polygon.push_back(0);
  polygons.push_back(polygon);

  res = PMP::internal::simplify_polygons_in_polygon_soup(points, polygons, K());
  std::cout << "res: " << res << " / size: " << polygons.back().size() << std::endl;
  assert(res == 1 && polygons.back().size() == 1);

  // ------
  polygon.clear();
  polygon.push_back(0); polygon.push_back(0); polygon.push_back(0); polygon.push_back(0);
  polygons.push_back(polygon);

  res = PMP::internal::simplify_polygons_in_polygon_soup(points, polygons, K());
  std::cout << "res: " << res << " / size: " << polygons.back().size() << std::endl;
  assert(res == 1 && polygons.back().size() == 1);

  // ------
  polygon.clear();
  polygon.push_back(0); polygon.push_back(3); polygon.push_back(3); polygon.push_back(3); polygon.push_back(0);
  polygons.push_back(polygon);

  res = PMP::internal::simplify_polygons_in_polygon_soup(points, polygons, K());
  std::cout << "res: " << res << " / size: " << polygons.back().size() << std::endl;
  assert(res == 1 && polygons.back().size() == 2);

  // ------
  polygon.clear();
  polygon.push_back(0); polygon.push_back(2); polygon.push_back(0); polygon.push_back(4);
  polygons.push_back(polygon);

  res = PMP::internal::simplify_polygons_in_polygon_soup(points, polygons, K());
  std::cout << "res: " << res << " / size: " << polygons.back().size() << std::endl;
  assert(res == 0 && polygons.back().size() == 4);

  // ------
  // Now with the same geometric positions, but different combinatorial information
  polygon.clear();
  polygon.push_back(0); polygon.push_back(2); polygon.push_back(1); polygon.push_back(6);
  polygons.push_back(polygon);

  res = PMP::internal::simplify_polygons_in_polygon_soup(points, polygons, K());
  std::cout << "res: " << res << " / size: " << polygons.back().size() << std::endl;
  assert(res == 1 && polygons.back().size() == 3);

  // ------
  polygon.clear();
  polygon.push_back(0); polygon.push_back(6); polygon.push_back(0); polygon.push_back(5);
  polygons.push_back(polygon);

  res = PMP::internal::simplify_polygons_in_polygon_soup(points, polygons, K());
  std::cout << "res: " << res << " / size: " << polygons.back().size() << std::endl;
  assert(res == 1 && polygons.back().size() == 2);

  // ------
  polygon.clear();
  polygon.push_back(0); polygon.push_back(6); polygon.push_back(5); polygon.push_back(3); polygon.push_back(3);
  polygons.push_back(polygon);

  res = PMP::internal::simplify_polygons_in_polygon_soup(points, polygons, K());
  std::cout << "res: " << res << " / size: " << polygons.back().size() << std::endl;
  assert(res == 1 && polygons.back().size() == 3);
}

void test_remove_invalid_polygons(const bool /*verbose*/ = false)
{
  std::cout << "test remove_invalid_polygons... " << std::endl;

  // points are not actually needed since only the size of the polygons is considered
  std::vector<Point_3> points;
  std::vector<CGAL_polygon> polygons;

  std::size_t res = PMP::internal::remove_invalid_polygons_in_polygon_soup(points, polygons);
  assert(res == 0 && polygons.size() == 0);

  // non-trivial polygon
  CGAL_polygon polygon;
  polygon.push_back(0); polygon.push_back(2); polygon.push_back(4);
  polygons.push_back(polygon);

  res = PMP::internal::remove_invalid_polygons_in_polygon_soup(points, polygons);
  assert(res == 0 && polygons.size() == 1);

  // another non-trivial polygon
  polygon.clear();
  polygon.push_back(1); polygon.push_back(3); polygon.push_back(3); polygon.push_back(7);
  polygons.push_back(polygon);

  // empty polygon
  polygon.clear();
  polygons.push_back(polygon);

  // 1-vertex polygon
  polygon.push_back(0);
  polygons.push_back(polygon);

  // 2-vertex polygon
  polygon.push_back(1);
  polygons.push_back(polygon);

  // another non-trivial polygon
  polygon.clear();
  polygon.push_back(8); polygon.push_back(9); polygon.push_back(7); polygon.push_back(6);
  polygons.push_back(polygon);

  res = PMP::internal::remove_invalid_polygons_in_polygon_soup(points, polygons);
  assert(res == 3 && polygons.size() == 3);
}

template <typename PointRange, typename PolygonRange>
std::size_t test_remove_isolated_points_data_set(PointRange& points,
                                                 PolygonRange& polygons,
                                                 const bool verbose = false)
{
  if(verbose)
  {
    std::cout << "Input:" << std::endl;
    std::cout << points.size() << " points" << std::endl;
    std::cout << polygons.size() << " polygons" << std::endl;
  }

  std::size_t rm_nv = PMP::remove_isolated_points_in_polygon_soup(points, polygons);

  if(verbose)
  {
    std::cout << "Removed " << rm_nv << " points" << std::endl;
    std::cout << points.size() << " points" << std::endl;
    std::cout << polygons.size() << " polygons" << std::endl;

    std::cout << "Polygons:" << std::endl;
    for(std::size_t i=0; i<polygons.size(); ++i)
    {
      for(std::size_t j=0; j<polygons[i].size(); ++j)
        std::cout << polygons[i][j] << " ";
      std::cout << std::endl;
    }
  }

  return rm_nv;
}

void test_remove_isolated_points(const bool verbose = false)
{
  std::cout << "test remove_isolated_points... " << std::endl;

  std::vector<Point_3> points;
  std::vector<CGAL_polygon> polygons;

  // everything empty
  std::size_t res = test_remove_isolated_points_data_set(points, polygons, verbose);
  assert(res == 0 && points.empty() && polygons.empty());

  points.push_back(Point_3(0,0,0));
  points.push_back(Point_3(1,2,0));
  points.push_back(Point_3(1,0,0));
  points.push_back(Point_3(1,3,0));
  points.push_back(Point_3(0,1,0));
  points.push_back(Point_3(1,1,0));

  // no polygons (all points are unused)
  std::deque<Point_3> points_copy(points.begin(), points.end());
  res = test_remove_isolated_points_data_set(points_copy, polygons, verbose);
  assert(res == 6 && points_copy.empty() && polygons.empty());

  CGAL_polygon polygon;
  polygon.push_back(0); polygon.push_back(2); polygon.push_back(4);
  polygons.push_back(polygon);

  polygon.clear();
  polygon.push_back(4); polygon.push_back(0); polygon.push_back(2); polygon.push_back(5);
  polygons.push_back(polygon);

  // generic test
  res = test_remove_isolated_points_data_set(points, polygons, verbose);
  assert(res == 2 && points.size() == 4 && polygons.size() == 2);
}

void test_slit_pinched_polygons(const bool /*verbose*/ = false)
{
  std::cout << "test split_pinched_polygons... " << std::endl;

  std::vector<Point_3> points;
  std::vector<CGAL_polygon> polygons;

  // everything empty
  std::size_t res = PMP::internal::split_pinched_polygons_in_polygon_soup<K>(points, polygons);
  assert(res == 0 && points.empty() && polygons.empty());

  points.push_back(Point_3(0,0,0)); // #0
  points.push_back(Point_3(1,0,0)); // #1
  points.push_back(Point_3(0,1,0)); // #2
  points.push_back(Point_3(1,1,0)); // #3
  points.push_back(Point_3(1,3,0)); // #4
  points.push_back(Point_3(1,1,0)); // #5

  // no pinch
  CGAL_polygon polygon;
  polygon.push_back(0); polygon.push_back(2); polygon.push_back(4);
  polygons.push_back(polygon);

  res = PMP::internal::split_pinched_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 0 && polygons.size() == 1);

  // pinch via same ID (2)
  polygon.clear();
  polygons.clear();
  polygon.push_back(0); polygon.push_back(2); polygon.push_back(4); polygon.push_back(3); polygon.push_back(2);
  polygons.push_back(polygon);

  res = PMP::internal::split_pinched_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 1 && polygons.size() == 2);
  assert(polygons[0].size() == 3 && polygons[1].size() == 2);

  // pinch via same point (5 & 3)
  polygon.clear();
  polygons.clear();
  polygon.push_back(5); polygon.push_back(1); polygon.push_back(0);
  polygon.push_back(3); polygon.push_back(4); polygon.push_back(2);
  polygons.push_back(polygon);

  res = PMP::internal::split_pinched_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 1 && polygons.size() == 2);
  assert(polygons[0].size() == 3 && polygons[1].size() == 3);

  // pinch on last point
  polygon.clear();
  polygons.clear();
  polygon.push_back(1); polygon.push_back(5); polygon.push_back(0);
  polygon.push_back(2); polygon.push_back(4); polygon.push_back(3);
  polygons.push_back(polygon);

  res = PMP::internal::split_pinched_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 1 && polygons.size() == 2);
  assert(polygons[0].size() == 4 && polygons[1].size() == 2);

  // multiple pinches
  polygon.clear();
  polygons.clear();
  polygon.push_back(5); polygon.push_back(1); polygon.push_back(3); // pinch 5&3, pinch 3...
  polygon.push_back(2); polygon.push_back(4); polygon.push_back(0);
  polygon.push_back(1); polygon.push_back(3); polygon.push_back(1); // ... and 3, pinch 1...
  polygon.push_back(2); polygon.push_back(4); polygon.push_back(1); // ... and 1
  polygons.push_back(polygon);

  res = PMP::internal::split_pinched_polygons_in_polygon_soup(points, polygons, K());
  assert(res == 3 && polygons.size() == 4);
  assert(polygons[0].size() == 2); // 5 1
  assert(polygons[1].size() == 5); // 3 2 4 5 1
  assert(polygons[2].size() == 3); // 1 2 4
  assert(polygons[3].size() == 2); // 1 3
}

int main()
{
  // test compilation with different polygon soup types
  std::vector<Point_3> vpoints;
  std::vector<std::vector<std::size_t> > vpolygons;
  PMP::repair_polygon_soup(vpoints, vpolygons);

  std::vector<std::deque<std::size_t> > dpolygons;
  PMP::repair_polygon_soup(vpoints, dpolygons);

  std::deque<std::vector<std::size_t> > dvpolygons;
  PMP::repair_polygon_soup(vpoints, dvpolygons);

  std::deque<std::array<std::size_t, 3> > apolygons;
  PMP::repair_polygon_soup(vpoints, apolygons);

  // test functions
  test_polygon_canonicalization(true);
  test_merge_duplicate_points(false);
  test_merge_duplicate_polygons(false);
  test_simplify_polygons(false);
  test_remove_invalid_polygons(false);
  test_remove_isolated_points(false);
  test_slit_pinched_polygons(false);

  return EXIT_SUCCESS;
}
