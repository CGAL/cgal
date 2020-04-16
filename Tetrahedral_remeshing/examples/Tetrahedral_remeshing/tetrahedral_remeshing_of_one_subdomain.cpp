#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/Tetrahedral_remeshing/tetrahedral_remeshing_io.h>

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

  bool operator()(Remeshing_triangulation::Cell_handle c) const
  {
    return m_subdomain == c->subdomain_index();
  }
};

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/triangulation_two_subdomains.binary.cgal";
  const float target_edge_length = (argc > 2) ? atof(argv[2]) : 0.1f;

  std::ifstream input(filename, std::ios_base::in | std::ios_base::binary);
  if(!input)
    return EXIT_FAILURE;

  Remeshing_triangulation tr;
  load_binary_triangulation(input, tr);

  CGAL::tetrahedral_adaptive_remeshing(tr, target_edge_length,
      CGAL::parameters::cell_selector(Cells_of_subdomain(2)));

  std::ofstream ofile("output.binary.cgal", std::ios::out);
  save_binary_triangulation(ofile, tr);

  return EXIT_SUCCESS;
}

