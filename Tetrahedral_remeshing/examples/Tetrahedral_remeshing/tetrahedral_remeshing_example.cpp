#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>

#include <CGAL/Triangulation_3.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/Random.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;


typedef CGAL::Triangulation_3<K>                                       T3;
typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K, int> Remeshing_triangulation;
//todo : add specialization for Cell_base without info
// (does not compile with `void` instead of `int`)

bool generate_input(const std::size_t& n,
                    const char* filename)
{
  T3 tr;
  CGAL::Random rng;

  while (tr.number_of_vertices() < n)
    tr.insert(T3::Point(rng.get_double(-1., 1.), rng.get_double(-1., 1.), rng.get_double(-1., 1.)));

  std::ofstream oFileT(filename, std::ios::out);
  // writing file output;
  oFileT << tr;

  return (!oFileT.bad());
}

int main(int argc, char* argv[])
{
  generate_input(1000, "data/random_sphere_triangulation.cgal");

  const char* filename     = (argc > 1) ? argv[1] : "data/random_sphere_triangulation.cgal";
  float target_edge_length = (argc > 2) ? atof(argv[2]) : 0.1f;

  std::ifstream input(filename, std::ios::in);
  if (!input)
  {
    std::cerr << "File " << filename << " could not be found" << std::endl;
    return EXIT_FAILURE;
  }

  T3 t3;
  input >> t3;
  CGAL_assertion(t3.is_valid());

  Remeshing_triangulation tr;
  CGAL::Tetrahedral_remeshing::build_remeshing_triangulation(t3, tr);
  
  CGAL::tetrahedral_adaptive_remeshing(tr, target_edge_length);

  std::ofstream oFileT("output.tr.cgal", std::ios::out);
  // writing file output;
  oFileT << tr;

  return EXIT_SUCCESS;
}
