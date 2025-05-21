#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/IO/write_MEDIT.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>

#include <CGAL/draw_constrained_triangulation_3.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point>;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const auto filename = (argc > 1) ? argv[1]
                       : CGAL::data_file_path("meshes/mpi_and_sphere.off");

  CGAL::Surface_mesh<K::Point_3> mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "Error: cannot read file " << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Number of facets in " << filename << ": "
            << mesh.number_of_faces() << "\n";

  auto ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(mesh,
                CGAL::parameters::return_empty_on_invalid_input(true));

  if(ccdt.number_of_constrained_facets() == 0)
  {
    std::cerr << "Error: no constrained facets in the CDT.\n";
    std::cerr << "Invalid input.\n";
    return EXIT_SUCCESS;
  }

  std::cout << "Number of constrained facets in the CDT: "
            << ccdt.number_of_constrained_facets() << '\n';

  std::ofstream ofs(argc > 2 ? argv[2] : "out.mesh");
  ofs.precision(17);
  CGAL::IO::write_MEDIT(ofs, ccdt);

  CGAL::draw(ccdt);

  return EXIT_SUCCESS;
}
