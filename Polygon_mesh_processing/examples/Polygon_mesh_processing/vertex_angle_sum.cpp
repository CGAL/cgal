#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/vertex_angle_sum.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_3                                          Point;
typedef K::Vector_3                                         Vector;

typedef CGAL::Surface_mesh<Point>                           Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor        vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor          face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/eight.off");

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  auto v_angle_sums = mesh.add_property_map<vertex_descriptor, double>("v:angle_sum", 0).first;

  PMP::vertex_angle_sums(mesh, v_angle_sums);

  std::cout << "Vertex normals :" << std::endl;
  for(vertex_descriptor vd: vertices(mesh))
    std::cout << v_angle_sums[vd] << std::endl;

  return 0;
}
