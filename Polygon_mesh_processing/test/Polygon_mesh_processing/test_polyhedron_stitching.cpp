#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/stitch_polyhedron.h>

#include <iostream>
#include <fstream>
#include <set>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;

void test(const char* fname)
{
  std::cout << "Testing " << fname << "..." << std::flush;
  std::ifstream input(fname);
  assert(input);
  Polyhedron P;
  input >> P;

  assert(P.size_of_vertices()!=0);
  
  CGAL::Polygon_mesh_processing::stitch_polyhedron(P);

  std::ofstream output("output.off");
  output << P;
  output.close();

  P.normalize_border();
  assert(P.is_valid(false, 5));
  std::cout << "OK\n";
}

int main()
{
  test("data_stitching/full_border.off");
  test("data_stitching/full_border_quads.off");
  test("data_stitching/half_border.off");
  test("data_stitching/mid_border.off");
  test("data_stitching/multiple_incidence.off");
  test("data_stitching/incidence_3.off");
}

