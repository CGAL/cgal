// file: examples/Arrangement_2/example1.C


//#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_observer.h>

typedef CGAL::Gmpq                                    Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;
typedef Arrangement_2::Halfedge_const_handle          Halfedge_const_handle;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Point_location;

int main ()
{
  // Construct the arrangement using the specialized insertion functions.
  Arrangement_2  arr;
  Point_location pl(arr);
  
  Segment_2       cv1 (Point_2(0, 1), Point_2(6, 9));
  Segment_2       cv2 (Point_2(6, 9), Point_2(6, 1));
  Segment_2       cv3 (Point_2(0, 1), Point_2(6, 1));
  
  Segment_2       cv4 (Point_2(4, 4), Point_2(7, 7));
  Segment_2       cv5 (Point_2(7, 7), Point_2(8, 0));
  Segment_2       cv6 (Point_2(8, 0), Point_2(4, 4));
  
  Segment_2       cv7 (Point_2(5, 3), Point_2(5, 5));

  insert(arr, pl, cv1);
  insert(arr, pl, cv2);
  insert(arr, pl, cv3);
  insert(arr, pl, cv4);
  insert(arr, pl, cv5);
  insert(arr, pl, cv6);
  insert(arr, pl, cv7);

  CGAL::Object obj = pl.locate (Point_2 (Number_type (11, 2), 
					 Number_type (11, 2)));
  Halfedge_const_handle h;
  if (CGAL::assign(h, obj))
  {
    arr.remove_edge (arr.non_const_handle(h));
  }

  obj = pl.locate (Point_2 (Number_type (11, 2), 
			    Number_type (5, 2)));
  if (CGAL::assign(h, obj))
  {
    arr.remove_edge (arr.non_const_handle(h));
  }
  
  // Print the arrangement vertices.
  Arrangement_2::Vertex_const_iterator  vit;
  Arrangement_2::Vertex_const_handle    vh;
  int                                   i, j;

  std::cout << arr.number_of_vertices() << " vertices:" << std::endl;
  for (i = 1, vit = arr.vertices_begin();
       vit != arr.vertices_end(); vit++, i++)
  {
    vh = vit;
    std::cout << '\t' << i << ": " << vh->point() << std::endl;
  }
  std::cout << std::endl;

  // Print the arrangement edges.
  Arrangement_2::Edge_const_iterator    eit;
  Arrangement_2::Halfedge_const_handle  hh;

  std::cout << arr.number_of_edges() << " edges:" << std::endl;
  for (i = 1, eit = arr.edges_begin(); eit != arr.edges_end(); eit++, i++)
  {
    hh = eit->handle();
    std::cout << '\t' << i << ": " << hh->curve() << std::endl;
  }
  std::cout << std::endl;

  // Print the arrangement faces.
  Arrangement_2::Face_const_iterator           fit;
  Arrangement_2::Face_const_handle             fh;
  Arrangement_2::Ccb_halfedge_const_circulator ccb;
  Arrangement_2::Holes_const_iterator          hoit;

  std::cout << arr.number_of_faces() << " faces." << std::endl;
  for (i = 1, fit = arr.faces_begin(); fit != arr.faces_end(); fit++, i++)
  {
    // Print the outer boundary of the face.
    fh = fit;
    std::cout << '\t' << i << ": ";
    if (fh->is_unbounded())
    {
      std::cout << "Unbounded face." << std::endl;
    }
    else
    {
      ccb = fh->outer_ccb();
      std::cout << ccb->source()->point();
      do
      {
        std::cout << " -> " << ccb->target()->point();
        ccb++;
      } while (ccb != fh->outer_ccb());
      std::cout << std::endl;
    }

    // Print the holes.
    for (j = 1, hoit = fh->holes_begin(); hoit != fh->holes_end(); hoit++, j++)
    {
      std::cout << "\t\tHole " << i << ": ";
      ccb = *hoit;
      std::cout << ccb->source()->point();
      do
      {
        std::cout << " -> " << ccb->target()->point();
        ccb++;
      } while (ccb != *hoit);
      std::cout << std::endl;
    }
  }
  std::cout << std::endl;


  
  
  return (0);
}

