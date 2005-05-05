// file: examples/Arrangement_2/example2.C


//#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_observer.h>

typedef CGAL::Quotient<int>                           Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef CGAL::Arr_naive_point_location<Arrangement_2> Point_location;

class My_observer : public CGAL::Arr_observer<Arrangement_2>
{
public:

  My_observer (Arrangement_2& arr) :
    CGAL::Arr_observer<Arrangement_2> (arr)
  {}

  virtual void before_merge_face (Face_handle,
                                  Face_handle,
                                  Halfedge_handle e)
  {
    std::cout << "-> The removal of :  [ " << e.curve()
	      << " ]  causes two faces to merge." << std::endl;
  }

  virtual void before_merge_edge (Halfedge_handle e1,
                                  Halfedge_handle e2,
                                  const X_monotone_curve_2& c)
  {
    std::cout << "-> Merging  [ " << e1.curve()
	      << " ]  and  [ " << e2.curve()
	      << " ]  to form  [ " << c << " ]  ." << std::endl;
  }
};

int main ()
{
  // Construct the arrangement containing one diamond-shaped face.
  Arrangement_2  arr;
  My_observer    obs (arr);
  Point_location pl (arr);

  Segment_2      cv1 (Point_2(-1, 0), Point_2(0, 1));
  Segment_2      cv2 (Point_2(0, 1), Point_2(1, 0));
  Segment_2      cv3 (Point_2(1, 0), Point_2(0, -1));
  Segment_2      cv4 (Point_2(0, -1), Point_2(-1, 0));

  arr_insert_non_intersecting (arr, pl, cv1);
  arr_insert_non_intersecting (arr, pl, cv2);
  arr_insert_non_intersecting (arr, pl, cv3);
  arr_insert_non_intersecting (arr, pl, cv4);

  std::cout << "V = " << arr.number_of_vertices()
	    << ",  E = " << arr.number_of_edges() 
	    << ",  F = " << arr.number_of_faces() << std::endl;

  // Insert a vertical segment dividing the diamond into two:
  Segment_2      cv_vert (Point_2(0, -1), Point_2(0, 1));
  Arrangement_2::Halfedge_handle
                 he_vert = arr_insert_non_intersecting (arr, pl, cv_vert);

  std::cout << "V = " << arr.number_of_vertices()
	    << ",  E = " << arr.number_of_edges() 
	    << ",  F = " << arr.number_of_faces() << std::endl;

  // Insert a horizontal segment dividing the diamond into four:
  Segment_2      cv_horiz (Point_2(-1, 0), Point_2(1, 0));

  arr_insert (arr, pl, cv_horiz);

  std::cout << "V = " << arr.number_of_vertices()
	    << ",  E = " << arr.number_of_edges() 
	    << ",  F = " << arr.number_of_faces() << std::endl;

  // Now remove a portion of the vertical segment.
  arr_remove_edge (arr, he_vert);
 
  std::cout << "V = " << arr.number_of_vertices()
	    << ",  E = " << arr.number_of_edges() 
	    << ",  F = " << arr.number_of_faces() << std::endl;

  // Remove the other vertical segment.
  Arrangement_2::Edge_iterator    eit;
  Arrangement_2::Traits_2         traits;

  for (eit = arr.edges_begin(); eit != arr.edges_end(); eit++)
  {
    if (traits.is_vertical_2_object() ((*eit).curve()))
    {
      arr_remove_edge (arr, *eit);

      std::cout << "V = " << arr.number_of_vertices()
		<< ",  E = " << arr.number_of_edges() 
		<< ",  F = " << arr.number_of_faces() << std::endl;
      break;
    }
  }

  return (0);
}

