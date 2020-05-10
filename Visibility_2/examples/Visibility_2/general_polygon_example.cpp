
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <iostream>
#include <vector>

// Define the used kernel and arrangement
typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef Kernel::Point_2                                         Point_2;
typedef Kernel::Segment_2                                       Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                           Arrangement_2;
typedef Arrangement_2::Halfedge_const_handle                    Halfedge_const_handle;
typedef Arrangement_2::Face_handle                              Face_handle;

// Define the used visibility class
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2>  TEV;

int main() {
  // Defining the input geometry
  Point_2 p1(1,2), p2(12, 3), p3(19,-2), p4(12,6), p5(14,14), p6( 9,5);
  Point_2 h1(8,3), h2(10, 3), h3( 8, 4), h4(10,6), h5(11, 6), h6(11,7);
  std::vector<Segment_2> segments;
  segments.push_back(Segment_2(p1,p2));
  segments.push_back(Segment_2(p2,p3));
  segments.push_back(Segment_2(p3,p4));
  segments.push_back(Segment_2(p4,p5));
  segments.push_back(Segment_2(p5,p6));
  segments.push_back(Segment_2(p6,p1));

  segments.push_back(Segment_2(h1,h2));
  segments.push_back(Segment_2(h2,h3));
  segments.push_back(Segment_2(h3,h1));
  segments.push_back(Segment_2(h4,h5));
  segments.push_back(Segment_2(h5,h6));
  segments.push_back(Segment_2(h6,h4));

  // insert geometry into the arrangement
  Arrangement_2 env;
  CGAL::insert_non_intersecting_curves(env,segments.begin(),segments.end());

  //Find the halfedge whose target is the query point.
  //(usually you may know that already by other means)
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

  std::cout << "Boundary edges of the visibility region:" << std::endl;
  Arrangement_2::Ccb_halfedge_circulator curr = fh->outer_ccb();
  std::cout << "[" << curr->source()->point() << " -> " << curr->target()->point() << "]" << std::endl;
  while (++curr != fh->outer_ccb())
    std::cout << "[" << curr->source()->point() << " -> " << curr->target()->point() << "]"<< std::endl;
  return 0;
}
