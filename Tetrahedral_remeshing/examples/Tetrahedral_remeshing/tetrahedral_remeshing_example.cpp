#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include "tetrahedral_remeshing_generate_input.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;

int main(int argc, char* argv[])
{
  const double target_edge_length = (argc > 1) ? atof(argv[1]) : 0.1;
  const std::size_t nbv = (argc > 2) ? atoi(argv[2]) : 1000;

  Remeshing_triangulation tr;
  CGAL::Tetrahedral_remeshing::insert_random_points_in_cube(nbv, tr);

  /// A subdomain index 0 is considered outside and is not remeshed
  /// so we set finite cells to a non-zero `Subdomain_index`
  for (auto cell : tr.finite_cell_handles())
    cell->set_subdomain_index(1);

  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length);

  return EXIT_SUCCESS;
}
