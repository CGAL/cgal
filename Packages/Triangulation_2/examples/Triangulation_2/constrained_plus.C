// file examples/Triangulation_2/constrained_plus.C
#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/intersections.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

typedef CGAL::Simple_cartesian<double>  K1;
typedef CGAL::Filtered_kernel<K1>       K2;
struct K : public K2 {};

typedef CGAL::Triangulation_vertex_base_2<K>              Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>    Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>       TDS;
typedef CGAL::Exact_predicates_tag                        Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS,Itag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT> CDTplus;

typedef CDTplus::Point                                     Point;

int
main( )
{
   Point pt[12] = {
     Point(0,0), Point(0,4), 
     Point(1,0), Point(1,4), 
     Point(2,0), Point(2,4),
     Point(-1,1), Point(3,1),
     Point(-1,2), Point(3,2),
     Point(0.5,1), Point(2.5,1)
   };

  CDTplus  tr;
  int  nc = 0, nsubc=0;
  for(int j=0; j<11; j+=2){
    ++nc;
    tr.insert(pt[j],pt[j+1]);
  }

  nsubc = std::distance(tr.subconstraints_begin(), tr.subconstraints_end());
  std::cerr << nc << " constraints inserted" << std::endl;
  std::cerr << nsubc << " resulting subconstraints" << std::endl;
  return 0;
}
