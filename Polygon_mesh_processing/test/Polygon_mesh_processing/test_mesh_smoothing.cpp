#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/mesh_smoothing.h>
#include <boost/graph/graph_traits.hpp>


int main(){

  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

  const char* filename = "data/simple_polygon.off";
  std::ifstream input(filename);

  Mesh mesh;
  input >> mesh;
  input.close();

  boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  CGAL::Polygon_mesh_processing::smooth_angles(mesh);

  /*
  for (vertex_descriptor v : vertices(mesh))
  {

  }
  */




  //0.720343 0.5 0


  std::ofstream out("data/output.off");
  out << mesh;
  out.close();





}
