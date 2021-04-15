#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include "tetrahedral_remeshing_generate_input.h"

#include <iostream>
#include <fstream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;

int main(int argc, char* argv[])
{
  const double target_edge_length = (argc > 1) ? atof(argv[1]) : 0.1;
  const std::size_t nbv = (argc > 2) ? atoi(argv[2]) : 1000;

  Remeshing_triangulation tr;
  CGAL::Tetrahedral_remeshing::generate_input_one_subdomain(nbv, tr);

  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length);

  return EXIT_SUCCESS;
}
