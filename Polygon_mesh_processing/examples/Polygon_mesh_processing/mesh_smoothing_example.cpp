#define CGAL_PMP_SMOOTHING_VERBOSE
#define CGAL_PMP_SMOOTHING_DEBUG

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>

#include "glog/logging.h"

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef CGAL::Surface_mesh<K::Point_3>                          Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  google::InitGoogleLogging(argv[0]);

  const char* filename = argc > 1 ? argv[1] : "data/grid.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid .off file." << std::endl;
    return EXIT_FAILURE;
  }

  const unsigned int repeat = 1;
  const unsigned int nb_iterations = 50;

  for(unsigned int t=0 ; t<repeat; ++t)
  {
#if 0
//    std::cout << "Smooth areas..." << std::endl;
//    PMP::smooth_areas(mesh, PMP::parameters::number_of_iterations(nb_iterations)
//                                            .use_safety_constraints(false));

    std::cout << "Smooth angles..." << std::endl;
    PMP::smooth_angles(mesh, PMP::parameters::number_of_iterations(nb_iterations)
                                             .use_safety_constraints(false));
#else
    PMP::smooth(mesh, PMP::parameters::number_of_iterations(nb_iterations)
                                      .use_safety_constraints(false));
#endif
  }

  std::ofstream output("mesh_smoothed.off");
  output << mesh;

  return EXIT_SUCCESS;
}
