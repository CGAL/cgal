#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/distance.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Point_set_3.h>
#include <iostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef K::Point_3                                            Point;

typedef CGAL::Point_set_3<Point>                              Point_set;

typedef CGAL::Surface_mesh<Point>                             Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor            face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/eight.off");

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  const int points_per_face = (argc > 2) ? std::stoi(argv[2]) : 10;

  std::vector<Point> points;
  PMP::sample_triangle_mesh(mesh,
                            std::back_inserter(points),
                            CGAL::parameters::number_of_points_per_face(points_per_face));


  Point_set point_set;
  PMP::sample_triangle_mesh(mesh,
                            point_set.point_back_inserter());

  std::cout << point_set.number_of_points() << std::endl;
  return 0;
}
