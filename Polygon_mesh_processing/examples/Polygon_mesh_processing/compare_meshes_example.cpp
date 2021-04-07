#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

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
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor      vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor        face_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const char* filename1 = (argc > 1) ? argv[1] : "data/P.off";

  Surface_mesh mesh1, mesh2;
  if(!PMP::read_polygon_mesh(filename1, mesh1))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }
  mesh2 = mesh1;
  CGAL::Euler::add_center_vertex(*halfedges(mesh2).begin(),mesh2);
  std::vector<std::pair<face_descriptor, face_descriptor> > common;
  std::vector<face_descriptor> m1_only, m2_only;
  PMP::match_faces(mesh1, mesh2, std::back_inserter(common), std::back_inserter(m1_only), std::back_inserter(m2_only));
  std::cout<<"Faces only in m1 : "<<std::endl;
  for(const auto& f : m1_only)
  {
    std::cout<<f<<", ";
  }
  std::cout<<"\n Faces only in m2: "<<std::endl;
  for(const auto& f : m2_only)
  {
    std::cout<<f<<", ";
  }
  std::cout<<"\n Faces in both: "<<std::endl;
  for(const auto& f_pair : common)
  {
    std::cout<<f_pair.first<<", "<<f_pair.second<<";;"<<std::endl;
  }
  return 0;
}
