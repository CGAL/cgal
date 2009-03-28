#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/minkowski_sum_3.h>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;

int main() {

  Nef_polyhedron N0, N1;
  std::ifstream in("cube.nef3");
  in >> N0;
  std::cin >> N1;
  Nef_polyhedron result = CGAL::minkowski_sum_3(N0, N1);
}
