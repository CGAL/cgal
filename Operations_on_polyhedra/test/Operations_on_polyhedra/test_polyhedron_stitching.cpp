#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_stitching.h>

#include <iostream>
#include <fstream>
#include <set>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;

void test(const char* fname)
{
  std::cout << "Testing " << fname << "..." << std::flush;
  std::ifstream input(fname);
  Polyhedron P;
  input >> P;

  CGAL::polyhedron_stitching(P);

  std::ofstream output("output.off");
  output << P;
  output.close();

  assert(P.is_valid(false, 5));
  std::cout << "OK\n";
}

int main()
{
  test("F17.off");
  test("full_border.off");
  test("full_border_quads.off");
  test("half_border.off");
  test("mid_border.off");
  test("multiple_incidence.off");
}

