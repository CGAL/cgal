#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <algorithm>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using CDT = CGAL::Conforming_constrained_Delaunay_triangulation_3<K>;

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

  auto cdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3<CDT>(mesh);
  static_assert(std::is_same_v<decltype(cdt), CDT>);
  CDT cdt2(mesh);
  const auto nb_cstr_facets = cdt2.number_of_constrained_facets();

  assert(cdt.triangulation().number_of_vertices() == cdt2.triangulation().number_of_vertices());
  assert(cdt.number_of_constrained_facets() == cdt2.number_of_constrained_facets());
  assert(cdt.number_of_constrained_facets() > mesh.num_faces());

  auto tr = std::move(cdt).triangulation();
  assert(0 == cdt.triangulation().number_of_vertices());
  assert(tr.number_of_vertices() == cdt2.triangulation().number_of_vertices());

  auto nb = 0u;
  for([[maybe_unused]] auto _ : cdt2.constrained_facets()) {
    ++nb;
  }
  assert(nb == nb_cstr_facets);
  assert(nb == std::distance(cdt2.constrained_facets_begin(), cdt2.constrained_facets_end()));

}
