#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_triangulation_3.h>

#include <algorithm>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

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

  auto ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(mesh);

  std::cout << "Number of vertices in the CDT: "
            << ccdt.triangulation().number_of_vertices() << '\n'
            << "Number of constrained facets in the CDT: "
            << ccdt.number_of_constrained_facets() << '\n';

  CGAL::draw(ccdt.triangulation());

  return EXIT_SUCCESS;
}
