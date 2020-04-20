#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/optimal_bounding_box.h>

#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>

typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
typedef K::Point_3                                             Point;
typedef K::Vector_3                                            Vector;

typedef CGAL::Surface_mesh<Point>                              Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor   vertex_descriptor;

namespace CP = CGAL::parameters;

int main(int argc, char** argv)
{
  std::ifstream input((argc > 1) ? argv[1] : "data/pig.off");
  Surface_mesh sm;
  if (!input || !(input >> sm) || sm.is_empty())
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  // a typical call
  std::array<Point, 8> obb_points;
  CGAL::oriented_bounding_box(sm, obb_points);

  // one can associate positions to the vertices of the mesh without changing the mesh
  std::unordered_map<vertex_descriptor, Point> translated_positions;
  for(const vertex_descriptor v : vertices(sm))
    translated_positions[v] = sm.point(v) + Vector(1, 2, 3);

  CGAL::oriented_bounding_box(sm, obb_points,
                              CP::vertex_point_map(boost::make_assoc_property_map(translated_positions)));

  // using a range of points
  std::vector<Point> points;
  for(const vertex_descriptor v : vertices(sm))
    points.push_back(sm.point(v));
  CGAL::oriented_bounding_box(points, obb_points);

  // one can associate positions to the range without changing the range
  std::map<Point, Point> scaled_positions;
  for(const Point& p : points)
    scaled_positions[p] = p + (p - CGAL::ORIGIN);

  CGAL::oriented_bounding_box(points, obb_points,
                              CP::point_map(boost::make_assoc_property_map(scaled_positions)));

  return EXIT_SUCCESS;
}
