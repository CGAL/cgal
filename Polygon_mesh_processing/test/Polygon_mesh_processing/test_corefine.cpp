#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;


void test(const char* f1, const char* f2)
{
  std::cout << "Corefining " << f1
            << " and " << f2 << "\n";

  std::cout << "  with Surface_mesh\n";
  Surface_mesh sm1, sm2;
  std::ifstream input(f1);
  assert(input);
  input >> sm1;
  input.close();
  input.open(f2);
  assert(input);
  input >> sm2;
  input.close();

  CGAL::Polygon_mesh_processing::corefine(sm1, sm2);

  assert(sm1.is_valid());
  assert(sm2.is_valid());

  std::cout << "  with Polyhedron_3\n";
  Polyhedron_3 P, Q;
  input.open(f1);
  assert(input);
  input >> P;
  input.close();
  input.open(f2);
  assert(input);
  input >> Q;

  CGAL::Polygon_mesh_processing::corefine(P, Q);

  assert(P.is_valid());
  assert(Q.is_valid());
}
int main(int argc, char** argv)
{
  for(int i=0; i< (argc-1)/2;++i)
  {
    test(argv[2*i+1], argv[2*(i+1)]);
    test(argv[2*(i+1)], argv[2*i+1]);
  }
}
