#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

using K=CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh=CGAL::Surface_mesh<K::Point_3>;

namespace PMP=CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby.off");

  Mesh m;
  if (!CGAL::IO::read_polygon_mesh(fname, m)|| is_empty(m))
  {
    std::cerr << "ERROR: cannot read " << fname << "\n";
    exit(1);
  }

  Mesh kernel = PMP::kernel(m);;
  std::ofstream("kernel.off") << kernel;
}
