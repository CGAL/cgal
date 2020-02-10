#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/random_perturbation.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;

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

  const double max_size = (argc > 2) ? atof(argv[2]) : 0.02;

  namespace PMP = CGAL::Polygon_mesh_processing;
  PMP::random_perturbation(
    mesh,
    max_size,
    PMP::parameters::vertex_point_map(mesh.points()).geom_traits(K()));

  std::ofstream out("data/eight_perturbed.off");
  out.precision(17);
  out << mesh;
  out.close();

  return 0;
}
