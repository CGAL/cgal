#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
//#define CGAL_TETRAHEDRAL_REMESHING_DEBUG

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

  // Collect options
  std::size_t nb_runs = 5;

  Remeshing_triangulation tr;
  generate_input_one_subdomain(1000, tr);

  const int target_edge_length = (argc > 1) ? atoi(argv[1]) : 1;

  std::ofstream ofs0("in.mesh");
  CGAL::IO::write_MEDIT(ofs0, tr);
  ofs0.close();

  const Remeshing_triangulation tr_ref = tr; // Keep a reference for comparison

  std::vector<std::string> output_tr;
  output_tr.reserve(nb_runs);

  for(std::size_t i = 0; i < nb_runs; ++i)
  {
    std::cout << "Run " << i << " of " << nb_runs << std::endl;

    tr = tr_ref; // Reset triangulation to reference
    CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length);

    std::ostringstream oss;
    CGAL::IO::write_MEDIT(oss, tr);
    output_tr.push_back(oss.str());
    oss.clear();

    if(i == 0)
      continue; // skip first run, it is the reference
    else
    {
      if(0 != output_tr[i-1].compare(output_tr[i]))
      {
        std::cerr << "******************************************" << std::endl;
        std::cerr << "*** Run " << i << " differs from run " << i-1 << std::endl;
        assert(false); // This should not happen
      }
      else
      {
        std::cout << "******************************************" << std::endl;
        std::cout << "*** Run " << i << " is identical to run " << i - 1 << std::endl;
        std::cout << "******************************************" << std::endl;
      }
    }
  }

  return EXIT_SUCCESS;
}
