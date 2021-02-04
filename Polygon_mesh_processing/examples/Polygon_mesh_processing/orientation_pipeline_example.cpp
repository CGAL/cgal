#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup_extension.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/algorithm.h>
#include <CGAL/Timer.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef K::Point_3                                            Point_3;

typedef CGAL::Surface_mesh<Point_3>                           Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  const std::string input_filename = (argc < 2) ? CGAL::data_file_path("meshes/blobby-shuffled.off") : argv[1];
  const std::string reference_filename = (argc < 2) ? CGAL::data_file_path("meshes/blobby.off") : argv[2];

  std::vector<Point_3> points;
  std::vector<std::vector<std::size_t> > polygons;
  if(!CGAL::IO::read_polygon_soup(input_filename, points, polygons) ||
     points.size() == 0 || polygons.size() == 0)
  {
    std::cerr << "Error: can not read input file.\n";
    return 1;
  }

  Mesh ref1;
  if(!PMP::IO::read_polygon_mesh(reference_filename, ref1))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  std::cout << "Is the soup a polygon mesh ? : " << PMP::is_polygon_soup_a_polygon_mesh(polygons)  << std::endl;

  PMP::orient_triangle_soup_with_reference_triangle_mesh<CGAL::Sequential_tag>(ref1, points, polygons);

  std::cout << "And now ? : " << PMP::is_polygon_soup_a_polygon_mesh(polygons) << std::endl;

  PMP::duplicate_non_manifold_edges_in_polygon_soup(points, polygons);

  std::cout << "And now ? : " << PMP::is_polygon_soup_a_polygon_mesh(polygons) << std::endl;

  Mesh poly;
  PMP::polygon_soup_to_polygon_mesh(points, polygons, poly);

  typedef boost::property_map<Mesh, CGAL::dynamic_face_property_t<std::size_t> >::type Fccmap;
  Fccmap fccmap = get(CGAL::dynamic_face_property_t<std::size_t>(), poly);

  std::cout << PMP::connected_components(poly, fccmap) << " CCs before merge." << std::endl;
  PMP::merge_reversible_connected_components(poly);

  std::cout<<PMP::connected_components(poly, fccmap) << " remaining CCs." << std::endl;

  return 0;
}
