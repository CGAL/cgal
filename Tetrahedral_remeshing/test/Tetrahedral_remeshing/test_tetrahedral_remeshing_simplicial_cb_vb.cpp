#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Simplicial_mesh_cell_base_3.h>
#include <CGAL/Simplicial_mesh_vertex_base_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/Tetrahedral_remeshing.h>

#include <CGAL/tags.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Subdomain_index = int;
using Surface_patch_index = std::pair<int, int>;
using Curve_index = char;
using Corner_index = short;

using Cb = CGAL::Simplicial_mesh_cell_base_3<Subdomain_index, Surface_patch_index>;
using Vb = CGAL::Simplicial_mesh_vertex_base_3<K, Subdomain_index, Surface_patch_index,
                                                  Curve_index, Corner_index>;

using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Sequential_tag>;
using Triangulation = CGAL::Triangulation_3<K, Tds>;

using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Triangulation>;

template<typename T3>
void generate_input_one_subdomain(const std::size_t nbv, T3& tr)
{
  CGAL::Random rng;

  typedef typename T3::Point Point;
  std::vector<Point> pts;
  while (pts.size() < nbv)
  {
    const float x = rng.uniform_real(-1.f, 1.f);
    const float y = rng.uniform_real(-1.f, 1.f);
    const float z = rng.uniform_real(-1.f, 1.f);

    pts.push_back(Point(x, y, z));
  }
  tr.insert(pts.begin(), pts.end());

  for (typename T3::Cell_handle c : tr.finite_cell_handles())
    c->set_subdomain_index(1);

  assert(tr.is_valid(true));

#ifdef CGAL_TETRAHEDRAL_REMESHING_GENERATE_INPUT_FILES
  std::ofstream out("data/triangulation_one_subdomain.binary.cgal",
                    std::ios_base::out | std::ios_base::binary);
  CGAL::save_binary_triangulation(out, tr);
  out.close();
#endif
}

int main(int argc, char* argv[])
{
  std::cout << "CGAL Random seed = " << CGAL::get_default_random().get_seed() << std::endl;

  Triangulation tr;
  generate_input_one_subdomain(100, tr);

  const double target_edge_length = (argc > 1) ? atof(argv[1]) : 0.1;

  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length);

  return EXIT_SUCCESS;
}
