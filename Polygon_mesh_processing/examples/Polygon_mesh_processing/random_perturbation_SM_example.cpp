#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/random_perturbation.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef K::Point_3                                            Point;

typedef CGAL::Surface_mesh<Point>                             Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor  vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor    face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/eight.off";

  Surface_mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  const double max_size = (argc > 2) ? atof(argv[2]) : 0.02;

  PMP::random_perturbation(mesh, max_size, CGAL::parameters::vertex_point_map(mesh.points()).geom_traits(K()));

  CGAL::IO::write_polygon_mesh("data/eight_perturbed.off", mesh, CGAL::parameters::stream_precision(17));

  return 0;
}
