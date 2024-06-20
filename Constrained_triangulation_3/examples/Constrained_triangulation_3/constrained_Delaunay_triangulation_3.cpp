#include <CGAL/Constrained_Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/draw_triangulation_3.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using CDT = CGAL::Constrained_Delaunay_triangulation_3<K>;

int main(int argc, char* argv[])
{
  CGAL::Surface_mesh<K::Point_3> mesh;
  auto filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/mpi.off");
  std::ifstream in(filename);
  if(!in || !CGAL::IO::read_OFF(in, mesh)) {
    std::cerr << "Error: cannot read file " << filename << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Read " << mesh.number_of_vertices() << " vertices and "
            << mesh.number_of_faces() << " faces" << std::endl;
  CDT cdt(mesh);
  std::cout << "Number of vertices in the CDT: "
            << cdt.triangulation().number_of_vertices() << std::endl;
  // CGAL::draw(cdt.triangulation());
}
