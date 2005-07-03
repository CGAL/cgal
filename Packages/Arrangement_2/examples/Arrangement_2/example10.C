// file: examples/Arrangement_2/example10.C

//#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_walk_along_line_point_location.h>

typedef CGAL::Quotient<int>                                     Number_type;
typedef CGAL::Cartesian<Number_type>                            Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef Traits_2::Point_2                                       Point_2;
typedef Traits_2::X_monotone_curve_2                            Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                           Arrangement_2;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Point_location;
typedef Arrangement_2::Vertex_handle                            Vertex_handle;

int main ()
{
  // Construct the arrangement containing one face.
  Arrangement_2  arr;
  Point_location pl (arr);

  Segment_2      cv1 (Point_2 (1, 1), Point_2 (9, 1));
  Segment_2      cv2 (Point_2 (9, 1), Point_2 (5, 7));
  Segment_2      cv3 (Point_2 (5, 7), Point_2 (1, 1));

  insert_non_intersecting (arr, pl, cv1);
  insert_non_intersecting (arr, pl, cv2);
  insert_non_intersecting (arr, pl, cv3);

  // Insert several vertices in the interior of this face.
  insert_vertex (arr, pl, Point_2 (5, 7));
  insert_vertex (arr, pl, Point_2 (5, 1));
  insert_vertex (arr, pl, Point_2 (4, 3));
  insert_vertex (arr, pl, Point_2 (4, 5));
  insert_vertex (arr, pl, Point_2 (8, 5));
  insert_vertex (arr, pl, Point_2 (6, 3));

  std::cout << "V = " << arr.number_of_vertices()
	    << ",  E = " << arr.number_of_edges() 
	    << ",  F = " << arr.number_of_faces() << std::endl;

  Arrangement_2::Isolated_vertices_iterator   iv_iter;

  // Print the remaining isolated vertices.
  Arrangement_2::Face_iterator   fcit;
  int                            k = 1;

  for (fcit = arr.faces_begin(); fcit != arr.faces_end(); ++fcit, ++k)
  {
    if (fcit->isolated_vertices_begin() != fcit->isolated_vertices_end())
    {
      std::cout << "Isolated vertices in face no. " << k << ": ";
      for (iv_iter = fcit->isolated_vertices_begin();
	   iv_iter != fcit->isolated_vertices_end(); ++iv_iter)
      {
	std::cout << "(" << iv_iter->point() << ") ";
	CGAL_assertion (arr.incident_face (iv_iter->handle()) == fcit);
      }
      std::cout << std::endl;
    }
  }

  // Insert two curve going through some isolated vertices.
  Segment_2      cv4 (Point_2 (3, 0), Point_2 (9, 6));
  Segment_2      cv5 (Point_2 (4, 5), Point_2 (6, 3));

  std::cout << "cv4:"<< std::endl;
  insert (arr, pl, cv4);
  std::cout << "cv5:"<< std::endl;
  insert (arr, pl, cv5);
  
  // Remove all removable vertices.
  Arrangement_2::Vertex_iterator   vtit, next;
  
  for (vtit = arr.vertices_begin(); vtit != arr.vertices_end(); vtit = next)
  {
    next = vtit;
    ++next;

    remove_vertex (arr, vtit->handle());
  }

  std::cout << "V = " << arr.number_of_vertices()
	    << ",  E = " << arr.number_of_edges() 
	    << ",  F = " << arr.number_of_faces() << std::endl;

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
	++ccb;
      } while (ccb != fh->outer_ccb());
      std::cout << std::endl;
    }

    // Print the holes.
    for (j = 1, hoit = fh->holes_begin(); hoit != fh->holes_end(); ++hoit, j++)
    {
      std::cout << "\t\tHole " << i << ": ";
      ccb = *hoit;
      std::cout << ccb->source()->point();
      do
      {
	std::cout << " -> " << ccb->target()->point();
	++ccb;
      } while (ccb != *hoit);
      std::cout << std::endl;
    }
  }
  std::cout << std::endl;

  return (0);
}

