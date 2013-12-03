#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <istream>

typedef CGAL::Exact_predicates_exact_constructions_kernel               Kernel;
typedef Kernel::Point_2                                                 Point_2;
typedef Kernel::Segment_2                                               Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                              Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                                   Arrangement_2;
typedef Arrangement_2::Edge_const_iterator                              Edge_const_iterator;
typedef Arrangement_2::Face_handle                                      Face_handle;
typedef Arrangement_2::Ccb_halfedge_circulator                          Ccb_halfedge_circulator;
typedef CGAL::Simple_polygon_visibility_2<Arrangement_2, CGAL::Tag_false>
                                                                        NSPV;
typedef CGAL::Simple_polygon_visibility_2<Arrangement_2, CGAL::Tag_true>
                                                                        RSPV;

int main() {
  //create environment
  Point_2 p1(0, 4), p2(0, 0), p3(3, 2), p4(4, 0), p5(4, 4), p6(1, 2);
  Segment_2 s[6];
  s[0] = Segment_2(p1, p2);
  s[1] = Segment_2(p2, p3);
  s[2] = Segment_2(p3, p4);
  s[3] = Segment_2(p4, p5);
  s[4] = Segment_2(p5, p6);
  s[5] = Segment_2(p6, p1);
  Arrangement_2 env;
  CGAL::insert_non_intersecting_curves(env, &s[0], &s[6]);
  //locate the query point in the arrangement
  Point_2 query_point(0.5, 2);
  Arrangement_2::Face_const_handle face;
  CGAL::Arr_naive_point_location<Arrangement_2> pl(env);
  CGAL::Object obj = pl.locate(query_point);
  CGAL::assign(face, obj);
  //visibility query
  Arrangement_2 non_regular_output;
  NSPV non_regular_visibility(env);
  Face_handle non_regular_fh = non_regular_visibility.compute_visibility(query_point, face, non_regular_output);
  std::cout << "Non-regularized visibility region of q has "
            << non_regular_output.number_of_edges()
            << " edges:" << std::endl;
  for (Edge_const_iterator eit = non_regular_output.edges_begin(); eit != non_regular_output.edges_end(); ++eit)
    std::cout << "[" << eit->curve() << "]" << std::endl;
  Arrangement_2 regular_output;
  RSPV regular_visibility(env);
  Face_handle regular_fh = regular_visibility.compute_visibility(query_point, face, regular_output);
  Ccb_halfedge_circulator curr = regular_fh->outer_ccb();
  std::cout << "Regularized visibility region of q has "
            << regular_output.number_of_edges()
            << " edges:" << std::endl;
  for (Edge_const_iterator eit = regular_output.edges_begin(); eit != regular_output.edges_end(); ++eit)
    std::cout << "[" << eit->curve() << "]" << std::endl;

  //For a regular face, we can also get its whole boundary by traversing its outer CCB.
  std::cout << "Traverse the face of the regularized visibility region:" << std::endl;
  std::cout << "[" << curr->curve() << "]"<< std::endl;
  while (++curr != regular_fh->outer_ccb())
    std::cout << "[" << curr->curve() << "]" << std::endl;
  return 0;
}

