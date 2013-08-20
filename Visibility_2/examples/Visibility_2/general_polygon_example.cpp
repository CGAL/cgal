#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/test_model_methods.h>
#include <CGAL/test_utils.h>
#include <CGAL/Triangular_expansion_visibility_2_.h>
#include <CGAL/test_utils.h>

#include <iostream>

typedef CGAL::Gmpq                                              Number_type;
typedef CGAL::Cartesian<Number_type>                            Kernel;
typedef Kernel::Point_2                                         Point_2;
typedef Kernel::Segment_2                                       Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                           Arrangement_2;
typedef Arrangement_2::Vertex_const_handle                      Vertex_const_handle;
typedef Arrangement_2::Halfedge_const_handle                    Halfedge_const_handle;
typedef Arrangement_2::Face_handle                              Face_handle;
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2, CGAL::Tag_true>
                                                                TEV;
int main() {
  Point_2 p1(1, 2), p2(12, 3), p3(19, -2), p4(12, 6), p5(14, 14), p6(9, 5);
  Point_2 h1(8,3), h2(10, 3), h3(8, 4), h4(10, 6), h5(11, 6), h6(11, 7);
  Segment_2 s[12];
  s[0] = Segment_2(p1, p2);
  s[1] = Segment_2(p2, p3);
  s[2] = Segment_2(p3, p4);
  s[3] = Segment_2(p4, p5);
  s[4] = Segment_2(p5, p6);
  s[5] = Segment_2(p6, p1);

  s[6] = Segment_2(h1, h2);
  s[7] = Segment_2(h2, h3);
  s[8] = Segment_2(h3, h1);
  s[9] = Segment_2(h4, h5);
  s[10] = Segment_2(h5, h6);
  s[11] = Segment_2(h6, h4);
  Arrangement_2 env;
  CGAL::insert_curves(env, &s[0], &s[12]);
  //find the halfedge whose target is the query point.
  Point_2 query_point = p4;
  Halfedge_const_handle he = env.halfedges_begin();
  while (he->source()->point() != p3 || he->target()->point() != p4)
    he++;
  //visibility query
  Arrangement_2 output_arr;
  TEV tev(env);
  Face_handle fh = tev.compute_visibility(query_point, he, output_arr);
  //print out the visibility region.
  std::cout << "Regularized visibility region of q has "
            << output_arr.number_of_edges()
            << " edges." << std::endl;
  Arrangement_2::Ccb_halfedge_circulator curr = fh->outer_ccb();
  std::cout << "Traverse the face of the visibility region." << std::endl;
  std::cout << "[" << curr->curve() << "]"<< std::endl;
  while (++curr != fh->outer_ccb())
    std::cout << "[" << curr->curve() << "]" << std::endl;
  return 0;
}
