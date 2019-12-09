#define CGAL_OPTIMAL_BOUNDING_BOX_DEBUG

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/optimal_bounding_box.h>

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

  // test one
  Surface_mesh obb_sm;
  CGAL::optimal_bounding_box(sm, obb_sm);

  std::ofstream("obb.off") << obb_sm;

  return EXIT_SUCCESS;
}
