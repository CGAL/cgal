#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <cstring>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_3                                       Point;

typedef CGAL::Surface_mesh<Point>                             SMesh;
typedef boost::graph_traits<SMesh>::vertex_descriptor         vertex_descriptor;
typedef boost::graph_traits<SMesh>::face_descriptor           face_descriptor;

int main()
{
  std::ifstream in("colored_tetra.ply");
  SMesh mesh;
  CGAL::IO::read_PLY(in, mesh);

  std::cerr << "Read mesh with " << mesh.number_of_vertices() << " vertices and "
            << mesh.number_of_faces() << " faces" << std::endl;

  std::cerr << "Properties associated with vertices:" << std::endl;
  std::vector<std::string> properties = mesh.properties<SMesh::Vertex_index>();
  for(std::size_t i = 0; i < properties.size(); ++ i)
    std::cerr << " * " << properties[i] << std::endl;

  std::cerr << "Properties associated with faces:" << std::endl;
  properties = mesh.properties<SMesh::Face_index>();
  for(std::size_t i = 0; i < properties.size(); ++ i)
    std::cerr << " * " << properties[i] << std::endl;

  mesh.add_property_map<SMesh::Edge_index, short>("id", 42);
  mesh.add_property_map<SMesh::Halfedge_index, float>("u", 13.f);
  mesh.add_property_map<SMesh::Halfedge_index, float>("v", 37.f);

  // Append second mesh
  std::ifstream in2("tetra.ply");
  CGAL::IO::read_PLY(in2, mesh);

  std::ofstream out("out.ply");
//  CGAL::IO::set_binary_mode(out);
  CGAL::IO::write_PLY(out, mesh);

  return 0;
}
