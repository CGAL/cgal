#include <CGAL/approximate_convex_segmentation.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <string>
#include "Utils.h" 

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

const double EPS = 1e-10;

template <typename Concurrency_tag>
bool test_on_convex_mesh(const std::string& path)
{
  Polyhedron mesh;
  if (!read_to_polyhedron(path.c_str(), mesh)) return false;    

  double concavity = CGAL::concavity_values<Concurrency_tag>(mesh);

  std::cout << "Concavity value: " << concavity << std::endl;
  if (concavity < 0 || concavity > EPS)
  {
    std::cerr << "Concavity value of a cube must be zero since it's convex" << std::endl;
    return false;
  }
  
  return true;
}

int main()
{
  expect_or_fail(test_on_convex_mesh<CGAL::Sequential_tag>("data/cube.off"));
#ifdef CGAL_LINKED_WITH_TBB
  expect_or_fail(test_on_convex_mesh<CGAL::Parallel_tag>("data/cube.off"));
#endif

  return EXIT_SUCCESS;
}

