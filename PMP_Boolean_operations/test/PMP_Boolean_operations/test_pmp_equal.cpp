#include <CGAL/Polygon_mesh_processing/equal.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

#include <cassert>
#include <iostream>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
typedef CGAL::Polyhedron_3<K> Polyhedron;

template <class Mesh>
void test()
{
  Mesh a, b;
  CGAL::make_hexahedron(K::Point_3(0, 0, 0), K::Point_3(1, 0, 0),
                        K::Point_3(1, 1, 0), K::Point_3(0, 1, 0),
                        K::Point_3(0, 0, 1), K::Point_3(1, 0, 1),
                        K::Point_3(1, 1, 1), K::Point_3(0, 1, 1), a);
  CGAL::make_hexahedron(K::Point_3(0, 0, 0), K::Point_3(1, 0, 0),
                        K::Point_3(1, 1, 0), K::Point_3(0, 1, 0),
                        K::Point_3(0, 0, 1), K::Point_3(1, 0, 1),
                        K::Point_3(1, 1, 1), K::Point_3(0, 1, 1), b);

  // Identical meshes.
  assert(PMP::are_equal_meshes(a, b));
  assert(PMP::are_equal_meshes(b, a));

  // Different geometry.
  PMP::transform(K::Aff_transformation_3(CGAL::TRANSLATION,
                                         K::Vector_3(2, 0, 0)), b);
  assert(!PMP::are_equal_meshes(a, b));

  PMP::transform(K::Aff_transformation_3(CGAL::TRANSLATION,
                                         K::Vector_3(-2, 0, 0)), b);
  assert(PMP::are_equal_meshes(a, b));

  // Opposite orientation.
  PMP::reverse_face_orientations(b);
  assert(!PMP::are_equal_meshes(a, b));
  PMP::reverse_face_orientations(b);
  assert(PMP::are_equal_meshes(a, b));

  // Different vertex/face indexing.
  Mesh c;
  std::vector<K::Point_3> points;
  std::vector<std::vector<std::size_t> > polygons;
  PMP::polygon_mesh_to_polygon_soup(a, points, polygons);

  std::reverse(points.begin(), points.end());
  // Reindex polygons according to the reversed point order.
  for (auto& p : polygons)
    for (std::size_t& i : p)
      i = points.size() - 1 - i;

  PMP::polygon_soup_to_polygon_mesh(points, polygons, c);
  assert(PMP::are_equal_meshes(a, c));
  assert(PMP::are_equal_meshes(c, a));
}

void test_polygon_soup()
{
  std::vector<K::Point_3> points_a = {K::Point_3(0, 0, 0),
                                      K::Point_3(1, 0, 0),
                                      K::Point_3(0, 1, 0)};
  std::vector<K::Point_3> points_b = points_a;
  std::vector<std::vector<std::size_t> > polygons_a{{0, 1, 2}};
  std::vector<std::vector<std::size_t> > polygons_b{{0, 1, 2}};

  assert(PMP::are_equal_polygon_soups(points_a, polygons_a,
                                      points_b, polygons_b));

  // Different vertex indexing.
  std::swap(points_b[0], points_b[1]);
  for (auto& p : polygons_b)
    for (std::size_t& i : p)
      i = (i == 0 ? 1 : i == 1 ? 0 : i);

  assert(PMP::are_equal_polygon_soups(points_a, polygons_a,
                                      points_b, polygons_b));

  // Different cyclic ordering of the vertices of a polygon.
  points_b = points_a;
  polygons_b = polygons_a;
  std::rotate(polygons_b[0].begin(),
              polygons_b[0].begin() + 1,
              polygons_b[0].end());

  assert(PMP::are_equal_polygon_soups(points_a, polygons_a,
                                      points_b, polygons_b));

  // Different orientation.
  points_b = points_a;
  polygons_b = polygons_a;
  std::reverse(polygons_b[0].begin(), polygons_b[0].end());

  assert(!PMP::are_equal_polygon_soups(points_a, polygons_a,
                                       points_b, polygons_b));

  // Different geometry.
  points_b = points_a;
  polygons_b = polygons_a;
  points_b[0] = K::Point_3(2, 0, 0);

  assert(!PMP::are_equal_polygon_soups(points_a, polygons_a,
                                       points_b, polygons_b));

  // Different number of polygons.
  points_b = points_a;
  polygons_b = polygons_a;
  polygons_b.push_back(polygons_b.front());

  assert(!PMP::are_equal_polygon_soups(points_a, polygons_a,
                                       points_b, polygons_b));
}

int main(int argc, char* argv[])
{
  if (argc == 3)
  {
    Surface_mesh a, b;
    if (!CGAL::IO::read_polygon_mesh(argv[1], a) ||
        !CGAL::IO::read_polygon_mesh(argv[2], b))
    {
      std::cerr << "Could not read polygon meshes\n";
      return EXIT_FAILURE;
    }

    const bool equal = PMP::are_equal_meshes(a, b);
    std::cout << (equal ? "Meshes are equal\n" : "Meshes are not equal\n");
    return equal ? EXIT_SUCCESS : EXIT_FAILURE;
  }

  test_polygon_soup();
  test<Surface_mesh>();
  test<Polyhedron>();
  return EXIT_SUCCESS;
}
