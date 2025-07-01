#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_3_to_lcc.h> 
#include <sstream>
#include <cassert>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC;
typedef CGAL::Polyhedron_3<LCC::Traits> Polyhedron;

int main()
{
  std::stringstream ss("OFF\n0 0 0\n");

  Polyhedron P;
  ss >> P;

  LCC lcc1, lcc2;

  auto d1 = CGAL::polyhedron_3_to_lcc(lcc1, P); 
  assert(d1 == LCC::null_descriptor);

  auto d2 = CGAL::import_from_polyhedron_3<LCC>(lcc2, P); 
  assert(d2 == LCC::null_descriptor);

  return 0;
}