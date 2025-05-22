#include <CGAL/approximate_convex_segmentation.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

#ifndef CGAL_LINKED_WITH_TBB
typedef CGAL::Sequential_tag Concurrency_tag;
#else
typedef CGAL::Parallel_tag Concurrency_tag;
#endif

int main(int argc, char* argv[])
{
  // read mesh
  Surface_mesh mesh;

  const char* filename = (argc > 1) ? argv[1] : "data/sword.off";
  std::ifstream input(filename);

  if (!input || !(input >> mesh))
  {
    std::cout << "Failed to read mesh" << std::endl;
    return EXIT_FAILURE;
  }

  if (CGAL::is_empty(mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cout << "Input mesh is invalid" << std::endl;
    return EXIT_FAILURE;
  }

  // compute concavity value
  double concavity = CGAL::concavity_values<Concurrency_tag>(mesh);

  // write result
  std::cout << "Concavity value: " << concavity << std::endl;

  return EXIT_SUCCESS;
}

