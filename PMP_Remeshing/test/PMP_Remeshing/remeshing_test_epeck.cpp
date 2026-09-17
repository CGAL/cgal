#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/Surface_mesh.h>

#include <fstream>
#include <iostream>

namespace PMP = CGAL::Polygon_mesh_processing;

using Epeck = CGAL::Exact_predicates_exact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<Epeck::Point_3>;

int main(int argc, const char* argv[])
{
  const std::string filename = (argc < 2) ? CGAL::data_file_path("meshes/sphere.off") : argv[1];

  std::ifstream input(filename);
  Mesh mesh;
  if (!input || !(input >> mesh))
  {
    std::cerr << "Error: cannot read surface mesh : " << filename << "\n";
    assert(false);
  }

  double target_edge_length = 1.0;
  PMP::isotropic_remeshing(CGAL::faces(mesh), target_edge_length, mesh);

  return 0;
}
