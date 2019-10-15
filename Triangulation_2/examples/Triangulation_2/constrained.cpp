#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <cassert>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Exact_predicates_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
typedef CDT::Point Point;
typedef CDT::Edge  Edge;
int
main( )
{
  CDT cdt;
  std::cout << "Inserting a grid of 5x5 constraints " << std::endl;
  for (int i = 1; i < 6; ++i)
    cdt.insert_constraint( Point(0,i), Point(6,i));
  for (int j = 1; j < 6; ++j)
    cdt.insert_constraint( Point(j,0), Point(j,6));

  assert(cdt.is_valid());
  int count = 0;
  for (const Edge& e : cdt.finite_edges())
    if (cdt.is_constrained(e))
      ++count;
  std::cout << "The number of resulting constrained edges is  ";
  std::cout <<  count << std::endl;
  return 0;
}
