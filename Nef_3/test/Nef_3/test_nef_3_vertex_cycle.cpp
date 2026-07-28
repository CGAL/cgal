#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_items.h>
#include <vector>
#include <cassert>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Nef_polyhedron_3<K, CGAL::SNC_items> Nef_polyhedron;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;

int main()
{
  std::vector<Point_3> cycle;
  cycle.push_back(Point_3(0,0,0));
  cycle.push_back(Point_3(1,0,0));
  cycle.push_back(Point_3(1,1,0));
  cycle.push_back(Point_3(0,1,0));

  Vector_3 normal(0,0,1);

  Nef_polyhedron nef(cycle.begin(), cycle.end(), normal);

  assert(nef.number_of_vertices() == 4);
  assert(nef.number_of_halfedges() == 10);

  return 0;
}
