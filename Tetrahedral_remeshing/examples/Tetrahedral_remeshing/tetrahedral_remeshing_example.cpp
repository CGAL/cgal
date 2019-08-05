#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K, int> Triangulation;
//todo : add specialization for Cell_base without info
// (does not compile with `void` instead of `int`)


int main(int argc, char* argv[])
{
  const char* filename     = (argc > 1) ? argv[1] : "data/pig.off";
  float target_edge_length = (argc > 2) ? atof(argv[2]) : 2.f;

  std::ifstream input(filename, std::ios::in);
  if (!input)
  {
    std::cerr << "File " << filename << " could not be found" << std::endl;
    return EXIT_FAILURE;
 }

  Triangulation tr;
  input >> tr;
  CGAL_assertion(tr.is_valid());

  CGAL::tetrahedral_adaptive_remeshing(tr, target_edge_length);

  std::ofstream oFileT("output.tr.cgal", std::ios::out);
  // writing file output;
  oFileT << tr;

  return EXIT_SUCCESS;
}
