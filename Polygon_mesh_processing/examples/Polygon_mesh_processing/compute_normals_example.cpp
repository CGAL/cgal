#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;

typedef CGAL::Surface_mesh<Point> Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;

int main()
{
  Surface_mesh mesh;
  std::ifstream input("data/eight.off");
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  Surface_mesh::Property_map<face_descriptor, Vector>   fnormals;
  Surface_mesh::Property_map<vertex_descriptor, Vector> vnormals;
  bool created;
  boost::tie(fnormals, created)
    = mesh.add_property_map<face_descriptor, Vector>("f:normals", Vector(0, 0, 0));
  boost::tie(vnormals, created)
    = mesh.add_property_map<vertex_descriptor, Vector>("v:normals", Vector(0, 0, 0));

  CGAL::Polygon_mesh_processing::compute_normals(mesh,
                                                 vnormals,
                                                 fnormals,
                                                 mesh.points(),//optional
                                                 K());         //optional

  std::cout << "Face normals :" << std::endl;
  BOOST_FOREACH(face_descriptor fd, faces(mesh)){
    std::cout << fnormals[fd] << std::endl;
  }
  std::cout << "Normals at vertices :" << std::endl;
  BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
    std::cout << vnormals[vd] << std::endl;
  }
  return 0;
}
