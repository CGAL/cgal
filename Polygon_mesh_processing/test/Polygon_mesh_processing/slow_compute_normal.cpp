#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Timer.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
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
  const std::string filename = "./data/slow_compute_normal.off";

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  CGAL::Timer t;
  t.start();
  std::cout << "compute_vertex_normal" << std::endl;
  PMP::compute_vertex_normal(vertex_descriptor(523), mesh);
  std::cout << t.time() << " sec." << std::endl;

  return 0;
}
