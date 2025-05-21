#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

int main(int argc, char** argv)
{
  if (argc==1)
  {
    std::cerr << "Nothing tested\n";
    return 1;
  }

  for(int i=0;i<(argc-1)/2; ++i)
  {
    std::cout << "Handling " << argv[2*i+1]
              << " expected res is " << argv[2*(i+1)] << "\n";

    std::ifstream input(argv[2*i+1]);
    assert(input);
    Surface_mesh sm;
    input >> sm;
    bool res = atoi(argv[2*(i+1)])>0;
    if (CGAL::Polygon_mesh_processing::does_bound_a_volume(sm)!=res)
      CGAL_error_msg("Result is not as expected (input orientation)");
    CGAL::Polygon_mesh_processing::reverse_face_orientations(sm);
    if (CGAL::Polygon_mesh_processing::does_bound_a_volume(sm)!=res)
      CGAL_error_msg("Result is not as expected (reversed orientation)");
  }
}
