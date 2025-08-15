/*!
  \ingroup PkgLinearCellComplexExamples
  \brief Example demonstrating octree generation from an in-memory FaceGraph
         and from a filename. The current 'regularized' flag is a stub
         (produces identical dart counts).
*/

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_octree_generation.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <fstream>
#include <iostream>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC_3;
typedef CGAL::Simple_cartesian<double>                     K;
typedef CGAL::Surface_mesh<K::Point_3>                     Surface_mesh;

int main()
{
  // Load mesh into a Surface_mesh (FaceGraph overload)
  std::ifstream in("data/meshes/cube.off");
  if(!in) {
    std::cerr << "Cannot open data/meshes/cube.off\n";
    return EXIT_FAILURE;
  }
  Surface_mesh sm;
  if(!(in >> sm)) {
    std::cerr << "Invalid OFF file\n";
    return EXIT_FAILURE;
  }
  if(!CGAL::is_triangle_mesh(sm))
    CGAL::Polygon_mesh_processing::triangulate_faces(sm);

  LCC_3 lcc_basic, lcc_reg;

  // Unregularized
  CGAL::compute_octree(lcc_basic, sm, 2, 5, false, false, false);
  // Regularized (stub: currently identical)
  CGAL::compute_octree(lcc_reg,   sm, 2, 5, false, false, true);

  std::cout << "Octree (basic) darts: " << lcc_basic.number_of_darts() << "\n";
  std::cout << "Octree (regularized) darts: " << lcc_reg.number_of_darts() << "\n";

  // Filename overload (const char*) (builds another LCC)
  LCC_3 lcc_from_file;
  CGAL::compute_octree(lcc_from_file, "data/meshes/cube.off", 2, 5, false, false, false);
  std::cout << "Octree (from filename) darts: " << lcc_from_file.number_of_darts() << "\n";

  return EXIT_SUCCESS;
}