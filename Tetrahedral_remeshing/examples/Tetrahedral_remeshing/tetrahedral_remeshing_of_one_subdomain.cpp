#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>

#include <CGAL/Triangulation_3.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/Random.h>

#include "tetrahedral_remeshing_io.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;

struct Cells_of_subdomain
{
private:
  const int m_subdomain;

public:
  Cells_of_subdomain(const int& subdomain)
    : m_subdomain(subdomain)
  {}

  const bool operator()(Remeshing_triangulation::Cell_handle c) const
  {
    return m_subdomain == c->subdomain_index();
  }
};

int main(int argc, char* argv[])
{
  float target_edge_length = (argc > 1) ? atof(argv[1]) : 0.1f;

  Remeshing_triangulation tr;
  generate_input(2, 1000, tr);

  CGAL::tetrahedral_adaptive_remeshing(tr, target_edge_length,
    CGAL::parameters::cell_selector(Cells_of_subdomain(2)));

  std::ofstream oFileT("output.binary.cgal", std::ios::out);
  save_binary_triangulation(oFileT, tr);

  std::cout << "done" << std::endl;
  return EXIT_SUCCESS;
}

