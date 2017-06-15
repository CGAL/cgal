#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <boost/property_map/property_map.hpp>

#include <map>
// #include <CGAL/Unique_hash_map.h>
// #include <boost/unordered_map.hpp>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;

typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::face_descriptor   face_descriptor;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/eight.off";
  std::ifstream input(filename);

  Polyhedron mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  std::map<face_descriptor,Vector> fnormals;
  std::map<vertex_descriptor,Vector> vnormals;
  // Instead of std::map you may use std::unordered_map, boost::unordered_map
  // or CGAL::Unique_hash_map
  // CGAL::Unique_hash_map<face_descriptor,Vector> fnormals;
  // boost::unordered_map<vertex_descriptor,Vector> vnormals;

  CGAL::Polygon_mesh_processing::compute_normals(mesh,
                                                 boost::make_assoc_property_map(vnormals),
                                                 boost::make_assoc_property_map(fnormals));

  std::cout << "Face normals :" << std::endl;
  for(face_descriptor fd: faces(mesh)){
    std::cout << fnormals[fd] << std::endl;
  }
  std::cout << "Vertex normals :" << std::endl;
  for(vertex_descriptor vd: vertices(mesh)){
    std::cout << vnormals[vd] << std::endl;
  }
  return 0;
}
