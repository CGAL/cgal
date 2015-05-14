#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/repair.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/degtri_sliding.off";
  std::ifstream input(filename);

  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  std::size_t nb
    = CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh);

  std::cerr << "There were " << nb << " degenerate faces in this mesh" << std::endl;
  mesh.collect_garbage();
  std::cout << mesh << std::endl;
  return 0;
}
