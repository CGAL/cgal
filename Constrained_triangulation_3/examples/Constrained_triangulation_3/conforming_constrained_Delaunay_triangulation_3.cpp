#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/IO/write_MEDIT.h>

#include <cassert>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

int main(int argc, char* argv[])
{
  auto filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/mpi.off");

  CGAL::Surface_mesh<K::Point_3> mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "Error: cannot read file " << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Read " << mesh.number_of_vertices() << " vertices and "
            << mesh.number_of_faces() << " faces" << std::endl;

  auto ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(mesh);

  //! [use of ccdt.triangulation()]
  std::cout << "Number of vertices in the CDT: "
            << ccdt.triangulation().number_of_vertices() << '\n';
  //! [use of ccdt.triangulation()]
  std::cout << "Number of constrained facets in the CDT: "
            << ccdt.number_of_constrained_facets() << '\n';

  std::ofstream ofs(argc > 2 ? argv[2] : "out.mesh");
  ofs.precision(17);
  CGAL::IO::write_MEDIT(ofs, ccdt);

  //! [move ccdt to tr]
  auto tr = std::move(ccdt).triangulation();
  // Now `tr` is a valid `CGAL::Triangulation_3` object that can be used for further processing.
  // and the triangulation of `ccdt` is empty.
  std::cout << "Number of vertices in the triangulation `tr`: "
            << tr.number_of_vertices() << '\n';
  std::cout << "Number of vertices in `ccdt`: "
            << ccdt.triangulation().number_of_vertices() << '\n';
  assert(ccdt.triangulation().number_of_vertices() == 0);
  //! [move ccdt to tr]

  return EXIT_SUCCESS;
}
