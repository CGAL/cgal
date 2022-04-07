#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include <boost/unordered_map.hpp>
#include <unordered_map>


typedef int                                           Number_type;
typedef CGAL::Simple_cartesian<Number_type>           Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Vertex_handle                  Vertex_handle;
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;

int main()
{
  Arrangement_2   arr;

  Segment_2       s1(Point_2(1, 3), Point_2(3, 5));
  Segment_2       s2(Point_2(3, 5), Point_2(5, 3));
  Segment_2       s3(Point_2(5, 3), Point_2(3, 1));
  Segment_2       s4(Point_2(3, 1), Point_2(1, 3));
  Segment_2       s5(Point_2(1, 3), Point_2(5, 3));

  Halfedge_handle e1 = arr.insert_in_face_interior(s1, arr.unbounded_face());
  Vertex_handle   v1 = e1->source();



  std::map<Vertex_handle,int> sm;
  sm[v1]= 1;
  std::unordered_map<Vertex_handle,int> sum;
  sum[v1]= 1;
  boost::unordered_map<Vertex_handle,int> bum;
  bum[v1]= 1;
  return 0;
}
