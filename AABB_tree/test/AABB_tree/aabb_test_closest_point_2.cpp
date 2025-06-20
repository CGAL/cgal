#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_rational.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_triangle_primitive_2.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <array>
#include <iostream>
#include <string>
#include <vector>

template<typename Kernel>
void test(const std::vector<CGAL::Simple_cartesian<double>::Point_2> &points, const std::vector<std::array<std::size_t, 3> > &faces) {
  using Point_2 = typename Kernel::Point_2;
  using Triangle_2 = typename Kernel::Triangle_2;
  using Iterator = typename std::vector<Triangle_2>::const_iterator;
  using Primitive = CGAL::AABB_triangle_primitive_2<Kernel, Iterator>;
  using Tree_traits = CGAL::AABB_traits_2<Kernel, Primitive>;
  using Tree = CGAL::AABB_tree<Tree_traits>;

  std::vector<Triangle_2> triangles(faces.size());
  for (std::size_t i = 0; i < faces.size(); ++i) {
    const auto& f = faces[i];
    triangles[i] = Triangle_2(Point_2(points[f[0]].x(), points[f[0]].y()), Point_2(points[f[1]].x(), points[f[1]].y()), Point_2(points[f[2]].x(), points[f[2]].y()));
  }

  Tree tree(triangles.begin(), triangles.end());

  // Without hint
  Point_2 query(-0.092372499264859229, -0.5067061545706153);
  Point_2 closest_point = tree.closest_point(query);
  std::cout << "Closest point to " << query << " is " << closest_point << std::endl;

  // With hint
  Point_2 hint(-0.077185400000000001, -0.42269299999999999);
  Point_2 closest_point_hint = tree.closest_point(query, hint);
  std::cout << "Closest point to " << query << " with hint " << hint << " is " << closest_point_hint << std::endl << std::endl;

  assert(closest_point == closest_point_hint);
}

int main(int argc, char** argv)
{
  std::cout.precision(17);

  // Read the input
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/camel.off");
  std::cout << "Reading " << filename << "..." << std::endl;

  std::vector<CGAL::Simple_cartesian<double>::Point_3> points;
  std::vector<std::array<std::size_t, 3> > faces;
  if (!CGAL::IO::read_polygon_soup(filename, points, faces) || faces.empty())
  {
    std::cerr << "Invalid input:" << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input: " << points.size() << " points, " << faces.size() << " faces" << std::endl;

  // Project onto the XY plane
  std::vector<CGAL::Simple_cartesian<double>::Point_2> points_2(points.size());
  for (std::size_t i = 0; i < points.size(); ++i)
    points_2[i] = CGAL::Simple_cartesian<double>::Point_2(points[i].x(), points[i].y());

  std::cout << "Testing closest point with Simple_cartesian<double>:" << std::endl;
  test<CGAL::Simple_cartesian<double> >(points_2, faces);
  std::cout << "Testing closest point with Epick:" << std::endl;
  test<CGAL::Exact_predicates_inexact_constructions_kernel>(points_2, faces);
  std::cout << "Testing closest point with Epeck:" << std::endl;
  test<CGAL::Exact_predicates_exact_constructions_kernel>(points_2, faces);
  std::cout << "Testing closest point with Simple_cartesian<Exact_rational>:" << std::endl;
  test<CGAL::Simple_cartesian<CGAL::Exact_rational>>(points_2, faces);

  std::cout << "Done." << std::endl;

  return EXIT_SUCCESS;
}
