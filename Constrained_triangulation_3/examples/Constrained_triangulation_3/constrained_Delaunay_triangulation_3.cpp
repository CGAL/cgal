#include <CGAL/make_constrained_Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/draw_triangulation_3.h>

#include <algorithm>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Tr = CGAL::default_constrained_triangulation_3_t<K>;

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

  Tr triangulation = CGAL::make_constrained_Delaunay_triangulation_3<Tr>(mesh);

  std::cout << "Number of vertices in the CDT: "
            << triangulation.number_of_vertices() << '\n'
            << "Number of constrained facets in the CDT: "
            << std::count_if(triangulation.finite_facets_begin(),
                             triangulation.finite_facets_end(),
                             [](const auto& f) {
                               auto [cell_handle, facet_index] = f;
                               return cell_handle->cdt_3_data().is_facet_constrained(facet_index);
                             }) << std::endl;

  // CGAL::draw(cdt.triangulation());
}
