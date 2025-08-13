/*!
  \ingroup PkgLinearCellComplexExamples
  \brief Example demonstrating octree generation from OFF mesh files.

  This example shows how to use the \cgal function `CGAL::compute_octree()` 
  to generate a hexahedral octree from an input OFF mesh file. The function 
  creates structured hexahedral subdivisions using AABB tree-based intersection 
  testing and supports configurable subdivision levels.
*/

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