#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_vertex_base_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_cell_base_3.h>

#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <algorithm>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Vbb = CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3<K>;
using Vb = CGAL::Conforming_constrained_Delaunay_triangulation_vertex_base_3<K, Vbb>;

using Cbb = CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3<K>;
using Cb = CGAL::Conforming_constrained_Delaunay_triangulation_cell_base_3<K, Cbb>;

using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using Tr = CGAL::Triangulation_3<K, Tds>;
using CDT = CGAL::Conforming_constrained_Delaunay_triangulation_3<K, Tr>;


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

  Tr tr = CGAL::convert_to_triangulation_3(std::move(cdt));

  CGAL::tetrahedral_isotropic_remeshing(tr, 2.,
                                        CGAL::parameters::number_of_iterations(3)
                                        .remesh_boundaries(false));

  std::cout << "Number of vertices in tr: " << tr.number_of_vertices() << std::endl;

  auto nb = 0u;
  for(auto f : tr.finite_facets())
  {
    const Tr::Cell_handle c = f.first;
    const int i = f.second;
    if(c->is_facet_on_surface(i))
      ++nb;
  }
  assert(nb == nb_cstr_facets);

  return EXIT_SUCCESS;
}
