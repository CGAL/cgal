//#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
//#define CGAL_TETRAHEDRAL_REMESHING_DEBUG
//#define CGAL_TETRAHEDRAL_REMESHING_GENERATE_INPUT_FILES

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/tetrahedral_remeshing_io.h>
#include <CGAL/Random.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>

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
    const float x = rng.uniform_real(-10.f, 10.f);
    const float y = rng.uniform_real(-10.f, 10.f);
    const float z = rng.uniform_real(-10.f, 10.f);

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

  Remeshing_triangulation tr;
  generate_input_one_subdomain(1000, tr);

  const int target_edge_length = (argc > 1) ? atoi(argv[1]) : 1;

  std::ofstream ofs0("in.mesh");
  CGAL::IO::write_MEDIT(ofs0, tr);
  ofs0.close();

  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length);

  std::ofstream ofs("out.mesh");
  CGAL::IO::write_MEDIT(ofs, tr);
  ofs.close();

  return EXIT_SUCCESS;
}
