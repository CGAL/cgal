#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel                 K;
typedef CGAL::Polygon_2<K>                                                Polygon_2;
typedef CGAL::Exact_intersections_tag                                     Itag_;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,CGAL::Default, Itag_> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT>                       CDTP;

typedef CDTP::Point                                                       Point;
typedef CDTP::Constraint_id                                               Cid;
typedef CDTP::Vertex_handle                                               Vertex_handle;


int main()
{
  CDTP cdtp;
  
  Vertex_handle handle1=cdtp.insert(Point(0,0));
  Vertex_handle handle2=cdtp.insert(Point(0,10));
  Vertex_handle handle3=cdtp.insert(Point(10,10));
  Vertex_handle handle4=cdtp.insert(Point(10,0));
  
  auto id = cdtp.insert_constraint(handle1,handle3);
  
  cdtp.remove_constraint(handle1,handle3);
  
  auto lastConstraintId=cdtp.insert_constraint(handle1,handle3);
  
  if (lastConstraintId == nullptr){
    std::cout << "problem" << std::endl;
    assert(false);
  }
  return 0;
}
