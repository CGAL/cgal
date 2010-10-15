#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <cassert>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_2<K>             Vbb;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>   Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>      TDS;
typedef CGAL::Exact_predicates_tag                       Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS,Itag> CDT;
typedef CGAL::Triangulation_hierarchy_2<CDT>             CDTH;
typedef CGAL::Constrained_triangulation_plus_2<CDTH>     Triangulation;

typedef Triangulation::Point                             Point;

int
main( )
{
  Triangulation cdt;
  std::cout << "Inserting a grid 5 x 5 of  constraints " << std::endl;
    for (int i = 1; i < 6; ++i)
    cdt.insert_constraint( Point(0,i), Point(6,i));
    for (int j = 1; j < 6; ++j)
    cdt.insert_constraint( Point(j,0), Point(j,6));

  int count = 0;
  for (Triangulation::Subconstraint_iterator scit = cdt.subconstraints_begin();
       scit != cdt.subconstraints_end();
       ++scit)  ++count;
  std::cout << "The number of resulting constrained edges is  ";
  std::cout <<  count << std::endl;

  //verbose mode of is_valid ; shows the number of vertices at each  level
  std::cout << "The number of vertices at successive levels" << std::endl;
  assert(cdt.is_valid(true));

  return 0;
}
