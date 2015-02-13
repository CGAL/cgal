#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/Connected_components.h>
#include <iostream>
#include <fstream>


typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                           Mesh;


typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

int main(int, char* argv[]) 
{
  Mesh sm;
  std::ifstream in(argv[1]);
  in >> sm;

  std::vector<face_descriptor> cc;
  face_descriptor fd = *faces(sm).first;
  CGAL::Polygon_mesh_processing::discover_connected_component(fd,
                                                              sm,
                                                              std::back_inserter(cc));

  std::cerr << cc.size() << " faces in the CC of " << fd << std::endl;
  return 0;
}
