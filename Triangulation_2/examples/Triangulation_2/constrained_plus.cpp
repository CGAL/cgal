#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <cassert>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;

typedef CGAL::Exact_intersections_tag                     Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,CGAL::Default,Itag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT>       CDTplus;
typedef CDTplus::Point                                    Point;

int
main( )
{
  CDTplus cdt;
  std::cout  << "Inserting a grid 5 x 5 of constraints " << std::endl;
  for (int i = 1; i < 6; ++i)
    cdt.insert_constraint( Point(0,i), Point(6,i));
  for (int j = 1; j < 6; ++j)
    cdt.insert_constraint( Point(j,0), Point(j,6));

  assert(cdt.is_valid());
  int count = 0;
  for (CDTplus::Subconstraint_iterator scit = cdt.subconstraints_begin();
       scit != cdt.subconstraints_end();
       ++scit)  ++count;
  std::cout << "The number of resulting constrained edges is  "
            <<  count << std::endl;
  return 0;
}
