#define CGAL_OPTIMAL_BOUNDING_BOX_DEBUG

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Real_timer.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
typedef K::Point_3                                             Point;
typedef CGAL::Surface_mesh<Point>                              Surface_mesh;

int main(int argc, char** argv)
{
  std::ifstream input(argv[1]);
  Surface_mesh sm;
  if (!input || !(input >> sm) || sm.is_empty())
  {
    std::cerr << argv[1] << " is not a valid off file.\n";
    return EXIT_FAILURE;
  }

  CGAL::Real_timer timer;
  timer.start();

  std::array<Point, 8> obb_points;
  CGAL::oriented_bounding_box(sm, obb_points);

  std::cout << "Elapsed time: " << timer.time() << std::endl;

  // Make a mesh out of the oriented bounding box
  Surface_mesh obb_sm;
  CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                        obb_points[4], obb_points[5], obb_points[6], obb_points[7], obb_sm);

  std::ofstream("obb.off") << obb_sm;

  return EXIT_SUCCESS;
}
