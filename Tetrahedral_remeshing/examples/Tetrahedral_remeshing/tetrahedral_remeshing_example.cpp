//#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include "tetrahedral_remeshing_io.h"

#include <iostream>
#include <fstream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;

template<typename T3>
bool generate_input_one_subdomain(const std::size_t nbv, T3& tr)
{
  CGAL::Random rng;
  std::cout << "CGAL Random seed = " << CGAL::get_default_random().get_seed() << std::endl;

  typedef typename T3::Point Point;
  while (tr.number_of_vertices() < nbv)
    tr.insert(Point(rng.get_double(-1., 1.), rng.get_double(-1., 1.), rng.get_double(-1., 1.)));

  for (typename T3::Cell_handle c : tr.finite_cell_handles())
    c->set_subdomain_index(1);

  std::string filename("data/triangulation_one_subdomain.binary.cgal");
  std::ofstream out(filename, std::ios_base::out | std::ios_base::binary);
  save_binary_triangulation(out, tr);

  return (!out.bad());
}

int main(int argc, char* argv[])
{
  Remeshing_triangulation tmp;
  generate_input_one_subdomain(1000, tmp);

  const char* filename     = (argc > 1) ? argv[1] : "data/triangulation_one_subdomain.binary.cgal";
  float target_edge_length = (argc > 2) ? atof(argv[2]) : 0.1f;

  std::ifstream input(filename, std::ios::in | std::ios::binary);

  Remeshing_triangulation t3;
  if (!input)
    return EXIT_FAILURE;

  if( !load_binary_triangulation(input, t3))
    return EXIT_FAILURE;

  CGAL::tetrahedral_adaptive_remeshing(t3, target_edge_length);

  // save output
  const std::string file_in(filename);

  // binary
  std::string file_out = file_in.substr(0, file_in.find_first_of("."));
  file_out.append("_out.binary.cgal");
  std::ofstream out(file_out.c_str(), std::ios_base::out | std::ios_base::binary);
  save_binary_triangulation(out, t3);

  return EXIT_SUCCESS;
}
