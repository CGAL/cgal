#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_mesh_processing/clip.h>

#include <CGAL/Real_timer.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<K::Point_3>;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  const std::filesystem::path fname = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby.off");
  const std::filesystem::path ofname = (argc > 2) ? argv[2] : fname.stem().concat("-clipped").concat(fname.extension().native());

  Mesh mesh;
  if (!CGAL::IO::read_polygon_mesh(fname, mesh) || is_empty(mesh))
  {
    std::cerr << "ERROR: cannot read " << fname << "\n";
    return EXIT_FAILURE;
  }

  CGAL::Real_timer timer;
  timer.start();

  PMP::Orthogonal_cut_plane_traits<K> cut_traits;
  PMP::Orthogonal_cut_plane_traits<K>::Plane_3 pl(0, 0);
  PMP::clip(mesh, pl, CGAL::parameters::geom_traits(cut_traits));

  timer.stop();

  std::cout << "Clip computation of " << fname << " done in " << timer.time() << std::endl;
  std::cout << "Output containing " << vertices(mesh).size() << " vertices written in " << ofname << std::endl;

  CGAL::IO::write_polygon_mesh(ofname, mesh, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
