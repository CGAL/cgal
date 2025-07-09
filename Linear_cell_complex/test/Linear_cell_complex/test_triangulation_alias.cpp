#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_3_to_lcc.h>
#include <cassert>
#include <cstdlib>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC_3;
typedef CGAL::Delaunay_triangulation_3<LCC_3::Traits> Triangulation;

int main()
{
  LCC_3 lcc1, lcc2;
  Triangulation T;

  assert(T.dimension() == -1);

  auto d1 = CGAL::triangulation_3_to_lcc(lcc1, T);
  assert(d1 == LCC_3::null_descriptor);

  auto d2 = CGAL::import_from_triangulation_3(lcc2, T);
  assert(d2 == LCC_3::null_descriptor);

  return EXIT_SUCCESS;
}

