#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

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
  std::ifstream in(filename);
  if(!in || !(in >> mesh)) {
    std::cerr << "Error: cannot read file " << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Number of facets in " << filename << ": "
            << mesh.number_of_faces() << "\n";

  auto ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(mesh,
                CGAL::parameters::check_preconditions(true));

  if(ccdt.number_of_constrained_facets() == 0)
  {
    std::cerr << "Error: no constrained facets in the CDT.\n";
    std::cerr << "Checking preconditions has failed.\n";
    return EXIT_SUCCESS;
  }

  std::cout << "Number of constrained facets in the CDT: "
            << ccdt.number_of_constrained_facets() << '\n';

  CGAL::draw(ccdt);

  return EXIT_SUCCESS;
}
