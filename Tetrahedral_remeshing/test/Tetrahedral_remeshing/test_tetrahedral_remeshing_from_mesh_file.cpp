//#define CGAL_TETRAHEDRAL_REMESHING_DEBUG
//#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
//#define CGAL_DUMP_REMESHING_STEPS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/SMDS_3/tet_soup_to_c3t3.h>
#include <CGAL/Random.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;
typedef Remeshing_triangulation Tr;


int test(const std::string& filename, const bool perturb_vertices_positions = false)
{
  const double target_edge_length = 0.05;

  std::ifstream in(filename);

  Remeshing_triangulation tr;
  if(CGAL::SMDS_3::build_triangulation_from_file(in, tr,
    false /*verbose*/,
    false /*replace_domain_0*/,
    false /*allow_non_manifold*/,
    true /*allow_negative_orientation*/))
    std::cout << "build triangulation ok" << std::endl;
  else
    return EXIT_FAILURE;

  if (perturb_vertices_positions)
  {
    std::cout << "CGAL Random seed = " << CGAL::get_default_random().get_seed() << std::endl;
    using P = Remeshing_triangulation::Point;
    const int percentage = 10;// 1/percentage = fraction of perturbed vertices
    const double max_move = 0.1 * target_edge_length;

    CGAL::Random rng;
    int vid = 0;
    for (auto v : tr.finite_vertex_handles())
    {
      if(++vid % percentage != 0)
        continue;

      const P& pv = v->point();
      v->set_point(P{pv.x() + rng.get_double(-1., 1.) * max_move,
                     pv.y() + rng.get_double(-1., 1.) * max_move,
                     pv.z() + rng.get_double(-1., 1.) * max_move});
    }
  }

  auto count_negative_cells = [&](const Tr& tr) -> unsigned int
    {
      return std::count_if(tr.finite_cells_begin(),
                           tr.finite_cells_end(),
                           [](const auto& c)
                           {
                             return CGAL::POSITIVE != CGAL::orientation(c.vertex(0)->point(),
                                                                        c.vertex(1)->point(),
                                                                        c.vertex(2)->point(),
                                                                        c.vertex(3)->point());
                           });
      };

  unsigned int neg = count_negative_cells(tr);
  if(neg > 0)
    tr.may_have_badly_oriented_cells(true);
  std::cout << neg << " finite cells are badly oriented before remeshing" << std::endl;

  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length,
    CGAL::parameters::number_of_iterations(5)
    .remesh_boundaries(false));

  unsigned int neg2 = count_negative_cells(tr);
  if(neg2 > 0)
    std::cout << neg2 << " finite cells are badly oriented after remeshing" << std::endl;

  in.close();
  return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
  std::cout << "Run remeshing from elephant.mesh" << std::endl;
  const std::string file1 = CGAL::data_file_path("meshes/elephant.mesh");
  if(test(file1) != EXIT_SUCCESS)
    return EXIT_FAILURE;

  std::cout << "\n\nRun remeshing from perturbed elephant.mesh" << std::endl;
  const std::string file2 = CGAL::data_file_path("meshes/elephant.mesh");
  if(test(file2, true) != EXIT_SUCCESS)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
