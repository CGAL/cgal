#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const char* filename = argc > 1 ? argv[1] : "data/pig.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid .off file." << std::endl;
    return EXIT_FAILURE;
  }

  const double time = 0.1;
  const unsigned int nb_iterations = 2;

  PMP::smooth_along_curvature_flow(mesh, time, PMP::parameters::number_of_iterations(nb_iterations));

  std::ofstream output("mesh_shape_smoothed.off");
  output << mesh;

  return EXIT_SUCCESS;
}
