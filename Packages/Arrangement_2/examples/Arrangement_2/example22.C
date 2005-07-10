// file: examples/Arrangement_2/example9.C


//#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>

typedef CGAL::Quotient<int>                           Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef CGAL::Arr_naive_point_location<Arrangement_2> Point_location;
typedef Arrangement_2::Vertex_handle                  Vertex_handle;

int main ()
{
  // Construct the arrangement containing one face.
  Arrangement_2  arr;
  Point_location pl (arr);

  Segment_2      cv1 (Point_2 (1, 1), Point_2 (5, 2));
  Segment_2      cv2 (Point_2 (5, 2), Point_2 (9, 1));
  Segment_2      cv3 (Point_2 (9, 1), Point_2 (9, 6));
  Segment_2      cv4 (Point_2 (9, 6), Point_2 (5, 5));
  Segment_2      cv5 (Point_2 (5, 5), Point_2 (1, 6));
  Segment_2      cv6 (Point_2 (1, 6), Point_2 (1, 1));

  insert_non_intersecting (arr, pl, cv1);
  insert_non_intersecting (arr, pl, cv2);
  insert_non_intersecting (arr, pl, cv3);
  insert_non_intersecting (arr, pl, cv4);
  insert_non_intersecting (arr, pl, cv5);
  insert_non_intersecting (arr, pl, cv6);

  std::cout << "Before inserting isolated vertices:" 
	    << "   V = " << arr.number_of_vertices()
	    << ",  E = " << arr.number_of_edges() 
	    << ",  F = " << arr.number_of_faces() << std::endl;

  // Get the bounded face, which forms the single hole inside the unbounded
  // face.
  CGAL_assertion (arr.number_of_faces() == 2);
  Arrangement_2::Face_handle    uf = arr.unbounded_face();
  Arrangement_2::Holes_iterator hole = uf->holes_begin();
  Arrangement_2::Face_handle    f = (*hole)->twin()->face();

  // Insert several vertices in the interior of this face.
  Vertex_handle  v1 = arr.insert_isolated_vertex (Point_2 (2, 4), f);
  Vertex_handle  v2 = arr.insert_isolated_vertex (Point_2 (3, 4), f);
  Vertex_handle  v3 = arr.insert_isolated_vertex (Point_2 (6, 4), f);
  Vertex_handle  v4 = arr.insert_isolated_vertex (Point_2 (3, 3), f);
  Vertex_handle  v5 = arr.insert_isolated_vertex (Point_2 (6, 3), f);

  std::cout << "After inserting isolated vertices:" 
	    << "   V = " << arr.number_of_vertices()
	    << ",  E = " << arr.number_of_edges() 
	    << ",  F = " << arr.number_of_faces() << std::endl;

  Arrangement_2::Isolated_vertices_iterator   iv_iter;

  std::cout << "Isolated vertices in the face: " << std::flush;
  for (iv_iter = f->isolated_vertices_begin();
       iv_iter != f->isolated_vertices_end(); ++iv_iter)
  {
    std::cout << "(" << iv_iter->point() << ") " << std::flush;
  }
  std::cout << std::endl;

  // Insert segments whose endpoints correspond to isolated vertices.
  Segment_2      seg1 (Point_2 (2, 4), Point_2 (3, 4));
  Segment_2      seg2 (Point_2 (6, 4), Point_2 (7, 4));

  arr.insert_at_vertices (seg1, v1, v2);
  arr.insert_from_left_vertex (seg2, v3);

  std::cout << "After connecting isolated vertices:" 
	    << "   V = " << arr.number_of_vertices()
	    << ",  E = " << arr.number_of_edges() 
	    << ",  F = " << arr.number_of_faces() << std::endl;

  std::cout << "Isolated vertices in the face: ";
  for (iv_iter = f->isolated_vertices_begin();
       iv_iter != f->isolated_vertices_end(); ++iv_iter)
    std::cout << "(" << iv_iter->point() << ") ";
  std::cout << std::endl;

  // Split the face into two and print the remaining isolated vertices.
  Segment_2       split_seg (Point_2 (5, 2), Point_2 (5, 5));

  insert_non_intersecting (arr, pl, split_seg);

  std::cout << "After Splitting the face:" 
	    << "   V = " << arr.number_of_vertices()
	    << ",  E = " << arr.number_of_edges() 
	    << ",  F = " << arr.number_of_faces() << std::endl;

  Arrangement_2::Face_iterator   fit;
  int                            i = 1;

  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit, ++i)
  {
    if (fit->isolated_vertices_begin() != fit->isolated_vertices_end())
    {
      std::cout << "Isolated vertices in face no. " << i << ": ";
      for (iv_iter = fit->isolated_vertices_begin();
	   iv_iter != fit->isolated_vertices_end(); ++iv_iter)
      {
	std::cout << "(" << iv_iter->point() << ") ";
	CGAL_assertion (arr.incident_face (iv_iter->handle()) == fit);
      }
      std::cout << std::endl;
    }
  }

  arr.remove_isolated_vertex (v4);

  return (0);
}

