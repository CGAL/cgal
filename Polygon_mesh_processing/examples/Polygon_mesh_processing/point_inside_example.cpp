#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <iostream>
#include <limits>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
typedef K::Point_3                                             Point;
typedef CGAL::Polyhedron_3<K>                                  Polyhedron;

namespace PMP = CGAL::Polygon_mesh_processing;

double max_coordinate(const Polyhedron& poly)
{
  double max_coord = -std::numeric_limits<double>::infinity();
  for(Polyhedron::Vertex_handle v : vertices(poly))
  {
    Point p = v->point();
    max_coord = (std::max)(max_coord, p.x());
    max_coord = (std::max)(max_coord, p.y());
    max_coord = (std::max)(max_coord, p.z());
  }
  return max_coord;
}

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/eight.off");

  Polyhedron poly;
  if(!PMP::IO::read_polygon_mesh(filename, poly) || CGAL::is_empty(poly) || !CGAL::is_triangle_mesh(poly))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  CGAL::Side_of_triangle_mesh<Polyhedron, K> inside(poly);

  double size = max_coordinate(poly);

  unsigned int nb_points = 100;
  std::vector<Point> points;
  points.reserve(nb_points);
  CGAL::Random_points_in_cube_3<Point> gen(size);
  for (unsigned int i = 0; i < nb_points; ++i)
    points.push_back(*gen++);

  std::cout << "Test " << nb_points << " random points in cube "
    << "[-" << size << "; " << size <<"]" << std::endl;

  int nb_inside = 0;
  int nb_boundary = 0;
  for (std::size_t i = 0; i < nb_points; ++i)
  {
    CGAL::Bounded_side res = inside(points[i]);

    if (res == CGAL::ON_BOUNDED_SIDE) { ++nb_inside; }
    if (res == CGAL::ON_BOUNDARY) { ++nb_boundary; }
  }

  std::cerr << "Total query size: " << points.size() << std::endl;
  std::cerr << "  " << nb_inside << " points inside " << std::endl;
  std::cerr << "  " << nb_boundary << " points on boundary " << std::endl;
  std::cerr << "  " << points.size() - nb_inside - nb_boundary << " points outside " << std::endl;

  return 0;
}
