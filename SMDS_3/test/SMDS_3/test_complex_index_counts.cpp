#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/IO/File_medit.h>

#include <cassert>
#include <sstream>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Tr = CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K>;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

int main()
{
  C3t3 c3t3;
  assert(c3t3.number_of_cells_in_subdomain(2) == 0);
  assert(c3t3.number_of_facets_in_patch(3) == 0);

  std::istringstream input(
    "MeshVersionFormatted 1\nDimension 3\n"
    "Vertices\n5\n"
    "0 0 0 0\n1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 -1 0\n"
    "Tetrahedra\n2\n1 2 3 4 2\n1 3 2 5 7\nEnd\n");
  const bool read = CGAL::IO::read_MEDIT(input, c3t3);
  assert(read);

  const C3t3& complex = c3t3;
  assert(complex.number_of_cells() == 2);
  assert(complex.number_of_cells_in_subdomain(2) == 1);
  assert(complex.number_of_cells_in_subdomain(7) == 1);
  assert(complex.number_of_cells_in_subdomain(0) == 0);
  assert(complex.number_of_cells_in_subdomain(9) == 0);

  Tr::Cell_handle cell = c3t3.cells_in_complex_begin();
  for(auto handle : c3t3.triangulation().finite_cell_handles())
  {
    c3t3.remove_from_complex(handle);
    c3t3.add_to_complex(handle, 2);
  }
  assert(complex.number_of_cells_in_subdomain(2) == 2);
  assert(complex.number_of_cells_in_subdomain(7) == 0);

  for(auto existing_facet : c3t3.triangulation().finite_facets())
    c3t3.remove_from_complex(existing_facet);
  auto facet = c3t3.triangulation().finite_facets_begin();
  const auto first_facet = *facet++;
  const auto second_facet = *facet++;
  const auto third_facet = *facet;
  c3t3.add_to_complex(first_facet, 3);
  c3t3.add_to_complex(second_facet, 3);
  c3t3.add_to_complex(third_facet, 5);
  assert(complex.number_of_facets() == 3);
  assert(complex.number_of_facets_in_patch(3) == 2);
  assert(complex.number_of_facets_in_patch(5) == 1);
  assert(complex.number_of_facets_in_patch(0) == 0);
  assert(complex.number_of_facets_in_patch(9) == 0);

  c3t3.remove_from_complex(second_facet);
  assert(complex.number_of_facets() == 2);
  assert(complex.number_of_facets_in_patch(3) == 1);
  c3t3.remove_from_complex(cell);
  assert(complex.number_of_cells() == 1);
  assert(complex.number_of_cells_in_subdomain(2) == 1);
}
