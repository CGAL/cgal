//#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/tetrahedral_remeshing_io.h>

#include <iostream>
#include <fstream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;

template<typename T3>
bool generate_input_one_subdomain(const std::size_t nbv, T3& tr)
{
  CGAL::Random rng;

  typedef typename T3::Point Point;
  while (tr.number_of_vertices() < nbv)
    tr.insert(Point(rng.get_double(-1., 1.), rng.get_double(-1., 1.), rng.get_double(-1., 1.)));

  for (typename T3::Cell_handle c : tr.finite_cell_handles())
    c->set_subdomain_index(1);

  std::string filename("data/triangulation_one_subdomain.binary.cgal");
  std::ofstream out(filename, std::ios_base::out | std::ios_base::binary);
  out << "binary CGAL c3t3\n";
  CGAL::set_binary_mode(out);
  out << tr;

  return (!out.bad());
}

int main(int argc, char* argv[])
{
  std::cout << "CGAL Random seed = " << CGAL::get_default_random().get_seed() << std::endl;

  Remeshing_triangulation tr;
  generate_input_one_subdomain(1000, tr);

  const float target_edge_length = (argc > 1) ? atof(argv[1]) : 0.1f;

  CGAL::tetrahedral_adaptive_remeshing(tr, target_edge_length);

  return EXIT_SUCCESS;
}
