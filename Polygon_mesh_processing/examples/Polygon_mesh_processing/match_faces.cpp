#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;

typedef K::Point_3                                                Point;

typedef CGAL::Surface_mesh<Point>                                 Surface_mesh;
typedef CGAL::Polyhedron_3<K>                                     Polyhedron;
typedef boost::graph_traits<Surface_mesh>::face_descriptor        face_descriptor_1;
typedef boost::graph_traits<Polyhedron>::face_descriptor          face_descriptor_2;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename1 = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/P.off");

  Surface_mesh mesh1;
  Polyhedron mesh2;
  if(!PMP::IO::read_polygon_mesh(filename1, mesh1))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }
  CGAL::copy_face_graph(mesh1, mesh2);
  CGAL::Euler::add_center_vertex(*halfedges(mesh2).begin(),mesh2);
  std::vector<std::pair<face_descriptor_1, face_descriptor_2> > common;
  std::vector<face_descriptor_1> m1_only;
  std::vector<face_descriptor_2> m2_only;

  PMP::match_faces(mesh1, mesh2, std::back_inserter(common), std::back_inserter(m1_only), std::back_inserter(m2_only));

  std::cout <<"Faces only in m1 :"<< std::endl;
  for(const face_descriptor_1& f : m1_only)
    std::cout << " " << f;

  std::cout <<"\n\nFaces only in m2:" << std::endl;
  for(const face_descriptor_2& f : m2_only)
    std::cout << " " << &(*f);

  std::cout << "\n\nFaces in both:" << std::endl;
  for(const auto& f_pair : common)
    std::cout << " (" << f_pair.first << ", " << &(*f_pair.second);
  std::cout << std::endl;

  return 0;
}
