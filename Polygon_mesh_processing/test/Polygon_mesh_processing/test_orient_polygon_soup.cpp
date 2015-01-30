#include <CGAL/Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/orient_polygon_soup.h>
#include <CGAL/polygon_soup_to_polygon_mesh.h>
#include <CGAL/IO/OFF_reader.h>

#include <string>
#include <fstream>
#include <iostream>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Polyhedron_3<K> Polyhedron;

void test(std::string fname, std::size_t expected_duplicated_vertices)
{
  std::vector<K::Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;
  std::ifstream input(fname.c_str());
  if (!input)
  {
    std::cerr << "Cannot open file " << fname << "\n";
    exit(EXIT_FAILURE);
  }
  
  if (!CGAL::read_OFF(input, points, polygons))
  {
    std::cerr << "Error parsing the OFF file " << fname << "\n";
    exit(EXIT_FAILURE);
  }
  
  std::size_t initial_nb_points = points.size();
  CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);

  assert(expected_duplicated_vertices == points.size()-initial_nb_points);

  Polyhedron P;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, P);
  assert(P.is_valid());
  std::cout << fname << " OK\n";
}


int main()
{
  test("data_polygon_soup/bad_cube.off", 4);
  test("data_polygon_soup/isolated_singular_vertex_one_cc.off", 1);
  test("data_polygon_soup/isolated_vertices.off", 1);
  test("data_polygon_soup/nm_vertex_and_edge.off", 6);
  test("data_polygon_soup/one_duplicated_edge.off", 6);
  test("data_polygon_soup/one_duplicated_edge_sharing_vertex.off", 8);
  test("data_polygon_soup/partial_overlap.off", 4);
  test("data_polygon_soup/incompatible_orientation.off", 2);
}