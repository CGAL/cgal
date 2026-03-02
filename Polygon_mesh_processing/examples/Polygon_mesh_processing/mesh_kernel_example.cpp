#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Real_timer.h>

using K=CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh=CGAL::Surface_mesh<K::Point_3>;

namespace PMP=CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby.off");
  const std::string oname = (argc > 2) ? argv[2] : "kernel.off";

  Mesh m;
  if (!CGAL::IO::read_polygon_mesh(fname, m)|| is_empty(m))
  {
    std::cerr << "ERROR: cannot read " << fname << "\n";
    exit(1);
  }

  CGAL::Real_timer timer;
  timer.start();
  Mesh kernel;
  PMP::kernel(m, kernel);
  timer.stop();

  std::cout << "Kernel computation of " << fname << " done in " << timer.time() << std::endl;
  std::cout << "Output containing " << vertices(kernel).size() << " vertices wrote in " << oname << std::endl;

  CGAL::IO::write_polygon_mesh("kernel.off", kernel, CGAL::parameters::stream_precision(17));

  return 0;
}
