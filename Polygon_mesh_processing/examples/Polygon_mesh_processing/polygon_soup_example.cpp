#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/OFF_reader.h>

#include <CGAL/orient_polygon_soup.h>
#include <CGAL/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/is_oriented.h>

#include <vector>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;

int main()
{
  std::ifstream input("data/elephant-shuffled.off");
  if (!input)
  {
    std::cerr << "Cannot open file " << std::endl;
    return 1;
  }

  std::vector<K::Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;
  if (!CGAL::read_OFF(input, points, polygons))
  {
    std::cerr << "Error parsing the OFF file " << std::endl;
    return 1;
  }

  std::size_t initial_nb_points = points.size();

  CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);

  Polyhedron mesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);

  if (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh))
    mesh.inside_out();

  std::ofstream out("elephant-oriented.off");
  out << mesh;
  out.close();

  return 0;
}
