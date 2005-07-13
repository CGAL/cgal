// file: examples/Arrangement_2/example14.C


//#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_observer.h>

typedef CGAL::Quotient<CGAL::MP_Float>                Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;
typedef CGAL::Arr_naive_point_location<Arrangement_2> Point_location;

class My_observer : public CGAL::Arr_observer<Arrangement_2>
{
public:

  My_observer (Arrangement_2& arr) :
    CGAL::Arr_observer<Arrangement_2> (arr)
  {}

  virtual void before_split_face (Face_handle,
                                  Halfedge_handle e)
  {
    std::cout << "-> New face, caused by: " << e->curve() << std::endl;
  }

  virtual void after_modify_edge (Halfedge_handle e)
  {
    std::cout << "-> Existing edge " << e->curve() 
	      << " has just been modified." << std::endl;
  }

};

int main ()
{
  // Construct the arrangement using the specialized insertion functions.
  Arrangement_2  arr;
  My_observer    obs (arr);

  Segment_2       cv1 (Point_2(1.0, 0.0), Point_2(3.0, 2.0));
  Segment_2       cv2 (Point_2(4.0, -1.0), Point_2(3.0, -2.0));
  Segment_2       cv3 (Point_2(4.0, -1.0), Point_2(1.0, 0.0));
  Segment_2       cv4 (Point_2(1.0, 0.0), Point_2(4.0, 1.0));
  Segment_2       cv5 (Point_2(3.0, 2.0), Point_2(4.0, 1.0));
  Segment_2       cv6 (Point_2(6.0, 0.0), Point_2(4.0, -1.0));
  Segment_2       cv7 (Point_2(4.0, 1.0), Point_2(6.0, 0.0));

  Halfedge_handle h1 = arr.insert_in_face_interior (cv1, arr.unbounded_face());
  Halfedge_handle h2 = arr.insert_in_face_interior (cv2, arr.unbounded_face());
  Halfedge_handle h3 = arr.insert_at_vertices (cv3, h2, h1->twin());
  Halfedge_handle h4 = arr.insert_from_left_vertex (cv4, h1->twin());
  Halfedge_handle h5 = arr.insert_at_vertices (cv5, h1, h4);
  Halfedge_handle h6 = arr.insert_from_left_vertex (cv6, h3->twin());
  Halfedge_handle h7 = arr.insert_at_vertices(cv7, h5, h6);

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
    std::cout << '\t' << i << ": " << eit->curve() << std::endl;
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

  // Perform point location.
  Point_location   pl (arr);
  Point_2          q (0, 0);
  CGAL::Object     obj;

  obj = pl.locate (q);
  if (CGAL::assign (fh, obj))
  {
    if (fh->is_unbounded())
      std::cout << "Inside unbounded face." << std::endl;
    else
      std::cout << "Inside face." << std::endl;
  }
  else if (CGAL::assign (hh, obj))
  {
    std::cout << "On halfedge: " << hh->curve() << std::endl;
  }
  else if (CGAL::assign (vh, obj))
  {
    std::cout << "On vertex: " << vh->point() << std::endl;
  }
  else
  {
    std::cout << "Illegal point-location result." << std::endl;    
  }

  // Insert additional segments.
  Segment_2       s1 (Point_2(-1.0, 0.0), Point_2(0.0, 2.0));
  Segment_2       s2 (Point_2(-1.0, 0.0), Point_2(-1.0, -2.0));
  Segment_2       s3 (Point_2(0.0, 2.0), Point_2(1.0, 0.0));
  Segment_2       s4 (Point_2(-1.0, -2.0), Point_2(3.0, -2.0));

  insert_non_intersecting (arr, pl, s1);
  insert_non_intersecting (arr, pl, s2);
  insert_non_intersecting (arr, pl, s3);
  insert_non_intersecting (arr, pl, s4);

  // Perform point location again.
  obj = pl.locate (q);
  if (CGAL::assign (fh, obj))
  {
    if (fh->is_unbounded())
      std::cout << "Inside unbounded face." << std::endl;
    else
      std::cout << "Inside face." << std::endl;
  }
  else if (CGAL::assign (hh, obj))
  {
    std::cout << "On halfedge: " << hh->curve() << std::endl;
  }
  else if (CGAL::assign (vh, obj))
  {
    std::cout << "On vertex: " << vh->point() << std::endl;
  }
  else
  {
    std::cout << "Illegal point-location result." << std::endl;    
  }

  // Test the insertion function (iis2 and iis3 cause some overlaps).
  Segment_2       iis1 (Point_2(0.0, -3.0), Point_2(5.0, 2.0));

  insert (arr, pl, iis1);

  std::cout << "V = " << arr.number_of_vertices()
	    << ",  E = " << arr.number_of_edges() 
	    << ",  F = " << arr.number_of_faces() << std::endl;
  
  Segment_2       iis2 (Point_2(-0.0, -2.0), Point_2(5.0, -2.0));

  insert (arr, pl, iis2);

  std::cout << "V = " << arr.number_of_vertices()
	    << ",  E = " << arr.number_of_edges() 
	    << ",  F = " << arr.number_of_faces() << std::endl;

  Segment_2       iis3 (Point_2(2.0, 2.0), Point_2(5.0, 0.5));

  insert (arr, pl, iis3);

  std::cout << "V = " << arr.number_of_vertices()
	    << ",  E = " << arr.number_of_edges() 
	    << ",  F = " << arr.number_of_faces() << std::endl;
 
  return (0);
}

