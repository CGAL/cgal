// file examples/Triangulation_2/constrained_plus.C
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

struct K : CGAL::Exact_predicates_exact_constructions_kernel {};

typedef CGAL::Triangulation_vertex_base_2<K>              Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>    Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>       TDS;
typedef CGAL::Exact_intersections_tag                     Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS,Itag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT>       CDTplus;
typedef CDTplus::Point                                    Point;

int
main( )
{
  CDTplus cdt;
  std::cerr << "Inserting a grid of constraints " << std::endl;
  std::cerr << "Inserting five horizontal constraints  " << std::endl;
  for (int i = 1; i < 6; ++i) 
    cdt.insert_constraint( Point(0,i), Point(6,i));
  std::cerr << "Inserting five vertical constraints   " << std::endl;
  for (int j = 1; j < 6; ++j) 
    cdt.insert_constraint( Point(j,0), Point(j,6));
  
  assert(cdt.is_valid());
  int count = 0;
  for (CDTplus::Subconstraint_iterator scit = cdt.subconstraints_begin();
       scit != cdt.subconstraints_end();
       ++scit)  ++count;
  std::cerr << "The number of resulting constrained edges is  ";
  std::cerr <<  count << std::endl;
  return 0;
}
