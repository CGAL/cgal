// file: examples/Arrangement_2/ex_arr_hist_1.C

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>

#include "point_location_utils.h"

typedef CGAL::Gmpq                                    Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Curve_2                             Segment_2;
typedef CGAL::Arrangement_with_history_2<Traits_2>    Arr_with_hist_2;
typedef Arr_with_hist_2::Curve_handle                 Curve_handle;
typedef CGAL::Arr_trapezoid_ric_point_location<Arr_with_hist_2>  
                                                      Point_location;

int main ()
{
  Arr_with_hist_2   arr;

  Segment_2         s1 (Point_2 (0, 3), Point_2 (4, 3));
  Curve_handle      c1 = arr.insert (s1);
  Segment_2         s2 (Point_2 (3, 1), Point_2 (3, 5));
  Curve_handle      c2 = arr.insert (s2);
  Segment_2         s3 (Point_2 (2, 3), Point_2 (5, 3));
  Curve_handle      c3 = arr.insert (s3);
  Segment_2         segs[3];

  segs[0] = Segment_2 (Point_2 (2, 6), Point_2 (7, 1));
  segs[1] = Segment_2 (Point_2 (0, 0), Point_2 (2, 6));
  segs[2] = Segment_2 (Point_2 (3, 4), Point_2 (3, 6));
  arr.insert (segs, segs + 3);

  std::cout << "Removing [" << s1 << "] : ";
  std::cout << arr.remove (c1) << " edges have been removed." << std::endl;

  // Print the arrangement.
  Arr_with_hist_2::Vertex_iterator  vit;

  std::cout << arr.number_of_vertices() << " vertices:" << std::endl;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
    std::cout << "(" << vit->point() << ")" << std::endl;

  // Print the arrangement edges.
  Arr_with_hist_2::Edge_iterator             eit;
  Arr_with_hist_2::Origin_curve_iterator     ocit;

  std::cout << arr.number_of_edges() << " edges:" << std::endl;
  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit)
  {
    std::cout << "[" << eit->curve() << "]. Origin: ";
    for (ocit = arr.origin_curves_begin (eit->handle());
	 ocit != arr.origin_curves_end (eit->handle()); ++ocit)
    {
      std::cout << " [" << *ocit << "]" << std::flush;
    }
    std::cout << std::endl;
  }

  // Test the split and merge functions:
  Segment_2         s4 (Point_2 (5, 6), Point_2 (7, 4));
  Curve_handle      c4 = arr.insert (s4);
  Arr_with_hist_2::Halfedge_handle   h = arr.split_edge (*(c4->edges_begin()), 
							 Point_2 (6, 5));

  std::cout << "V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;

  arr.merge_edge (h, h->next());
  std::cout << "V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;


  Point_location  pl (arr);

  point_location_query (pl, Point_2 (Number_type (7, 2), 4));
  point_location_query (pl, Point_2 (6, 2));
  point_location_query (pl, Point_2 (2, 4));

  return (0);
}

