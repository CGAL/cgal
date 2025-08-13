#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_octree_generation.h>
#include <iostream>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC_3;

int main()
{
  LCC_3 lcc;
  CGAL::compute_octree(lcc, "data/meshes/cube.off", 2, 5, false, false, true);
  std::cout<<"Octree generated with "<<lcc.number_of_darts()<<" darts"<<std::endl;
  return EXIT_SUCCESS;
}