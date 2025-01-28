#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>

#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                K;
typedef CGAL::Projection_traits_3<K>                                       GT;

typedef CGAL::Exact_predicates_tag                                         Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<GT, CGAL::Default, Itag> CDT;
typedef CDT::Point                                                         Point;
typedef CDT::Edge                                                          Edge;

int main()
{
  //4 points on the diagonal plane of a cube
  std::vector<Point> ps(4);
  ps[0] = Point(0,0,0);
  ps[1] = Point(3,1,-1);
  ps[2] = Point(-1, 3, -3);
  ps[3] = Point(1,0.5,-0.5);

  GT gt{ { 0, 1, 1} };
  CDT cdt(gt);
  for(int i = 0; i< 4; ++i)
    cdt.insert(ps[i]);

  for(int i = 1; i < 3; ++i)
    cdt.insert_constraint(ps[i], ps[i+1]);

  for(CDT::Face_handle f : cdt.all_face_handles())
  {
   for(int i=0; i<3; ++i)
     std::cout << f->vertex(i)->point() << " ";
   std::cout << std::endl;
  }
}
