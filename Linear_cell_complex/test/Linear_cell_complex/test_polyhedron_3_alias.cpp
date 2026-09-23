#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/import_face_graph_to_lcc.h>
#ifndef CGAL_NO_DEPRECATED_CODE
#include <CGAL/Polyhedron_3_to_lcc.h>
#endif
#include <sstream>
#include <cassert>
#include <cstdlib>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC;
typedef CGAL::Polyhedron_3<LCC::Traits> Polyhedron;


int main()
{
#ifndef CGAL_NO_DEPRECATED_CODE
  std::stringstream ss("OFF\n0 0 0\n");

  Polyhedron P;
  ss >> P;

  LCC lcc0, lcc1, lcc2;

  auto d0 = CGAL::import_face_graph_to_lcc(P, lcc0);
  assert(d0 == LCC::null_descriptor);

  auto d1 = CGAL::polyhedron_3_to_lcc(lcc1, P);
  assert(d1 == LCC::null_descriptor);

  auto d2 = CGAL::import_from_polyhedron_3<LCC>(lcc2, P);
  assert(d2 == LCC::null_descriptor);
#endif // CGAL_NO_DEPRECATED_CODE
  return EXIT_SUCCESS;
}

