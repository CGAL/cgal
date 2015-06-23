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
typedef boost::graph_traits<Surface_mesh>::face_descriptor   face_descriptor;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/eight.off";
  std::ifstream input(filename);

  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  Surface_mesh::Property_map<face_descriptor, Vector> fnormals =
    mesh.add_property_map<face_descriptor, Vector>
      ("f:normals", CGAL::NULL_VECTOR).first;
  Surface_mesh::Property_map<vertex_descriptor, Vector> vnormals =
    mesh.add_property_map<vertex_descriptor, Vector>
      ("v:normals", CGAL::NULL_VECTOR).first;

  CGAL::Polygon_mesh_processing::compute_normals(mesh,
        vnormals,
        fnormals,
        CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh.points()).
        geom_traits(K()));

  std::cout << "Face normals :" << std::endl;
  BOOST_FOREACH(face_descriptor fd, faces(mesh)){
    std::cout << fnormals[fd] << std::endl;
  }
  std::cout << "Vertex normals :" << std::endl;
  BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
    std::cout << vnormals[vd] << std::endl;
  }
  return 0;
}
