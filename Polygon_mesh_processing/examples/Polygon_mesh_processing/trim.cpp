#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

#include <CGAL/Polygon_mesh_processing/trim.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
typedef K::Point_3                                                      Point_3;

typedef CGAL::Surface_mesh<Point_3>                                     Mesh;
typedef typename boost::graph_traits<Mesh>::vertex_descriptor           vertex_descriptor;
typedef typename boost::graph_traits<Mesh>::face_descriptor             face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;


int main(int argc, char* argv[])
{
  const std::string mesh_filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/tooth.off");
  const std::string polyline_filename = (argc > 2) ? argv[2] : CGAL::data_file_path("polylines/tooth.polylines.txt");
  Mesh tm, other;
  CGAL::IO::read_polygon_mesh(mesh_filename, tm);

  // Read the polyline and create a mesh for the extruded polyline
  // The polyline file must start with the number of points
  // and the first and last point must be the same

  std::ifstream ifs(polyline_filename);
  int n;
  ifs >> n;

  Point_3 p;
  std::vector<Point_3> points;
  points.reserve(n);
  while(ifs >> p){
    points.push_back(p);
  }
  bool b = PMP::trim(tm, other, points);

  if(b){
    CGAL::IO::write_polygon_mesh("out0.off", tm);
    CGAL::IO::write_polygon_mesh("out1.off", other);
  }else{
    std::cout<< "failed" << std::endl;
  }
  return 0;
}
