#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
//#define CGAL_DUMP_REMESHING_STEPS
//#define CGAL_TETRAHEDRAL_REMESHING_DEBUG
//#define CGAL_TETRAHEDRAL_REMESHING_GENERATE_INPUT_FILES

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/tetrahedral_remeshing_io.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;

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

  CGAL_assertion(tr.is_valid(true));

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

  Remeshing_triangulation tr;
  generate_input_one_subdomain(1000, tr);

  const double target_edge_length = (argc > 1) ? atof(argv[1]) : 0.1;

  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length);

  return EXIT_SUCCESS;
}
