//! \file examples/Arrangement_on_surface_2/face_extension_overlay.cpp
// A face overlay of two arrangements with extended face records.

#include "arr_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>

typedef CGAL::Cartesian<Number_type>                            Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef Traits_2::Point_2                                       Point_2;
typedef Traits_2::X_monotone_curve_2                            Segment_2;
typedef CGAL::Arr_face_extended_dcel<Traits_2, bool>            Dcel;
typedef CGAL::Arrangement_2<Traits_2, Dcel>                     Arrangement_2;
typedef CGAL::Arr_face_overlay_traits<Arrangement_2,
                                      Arrangement_2,
                                      Arrangement_2,
                                      std::logical_and<bool> >  Overlay_traits;

int main ()
{
  // Construct the first arrangement, containing a square-shaped face.
  Arrangement_2          arr1;

  Segment_2      s1 (Point_2(2, 2), Point_2(6, 2));
  Segment_2      s2 (Point_2(6, 2), Point_2(6, 6));
  Segment_2      s3 (Point_2(6, 6), Point_2(2, 6));
  Segment_2      s4 (Point_2(2, 6), Point_2(2, 2));

  insert_non_intersecting_curve (arr1, s1);
  insert_non_intersecting_curve (arr1, s2);
  insert_non_intersecting_curve (arr1, s3);
  insert_non_intersecting_curve (arr1, s4);

  // Mark just the bounded face.
  Arrangement_2::Face_iterator   fit;

  CGAL_assertion (arr1.number_of_faces() == 2);
  for (fit = arr1.faces_begin(); fit != arr1.faces_end(); ++fit)
    fit->set_data (fit != arr1.unbounded_face());

  // Construct the second arrangement, containing a rhombus-shaped face.
  Arrangement_2          arr2;

  Segment_2      t1 (Point_2(4, 1), Point_2(7, 4));
  Segment_2      t2 (Point_2(7, 4), Point_2(4, 7));
  Segment_2      t3 (Point_2(4, 7), Point_2(1, 4));
  Segment_2      t4 (Point_2(1, 4), Point_2(4, 1));

  insert_non_intersecting_curve (arr2, t1);
  insert_non_intersecting_curve (arr2, t2);
  insert_non_intersecting_curve (arr2, t3);
  insert_non_intersecting_curve (arr2, t4);

  // Mark just the bounded face.
  CGAL_assertion (arr2.number_of_faces() == 2);
  for (fit = arr2.faces_begin(); fit != arr2.faces_end(); ++fit)
    fit->set_data (fit != arr2.unbounded_face());

  // Compute the overlay of the two arrangements, marking only the faces that
  // are intersections of two marked faces in arr1 and arr2, respectively.
  Arrangement_2          overlay_arr;
  Overlay_traits         overlay_traits;

  overlay (arr1, arr2, overlay_arr, overlay_traits);

  // Go over the faces of the overlaid arrangement and print just the marked
  // ones.
  Arrangement_2::Ccb_halfedge_circulator    curr;

  std::cout << "The union is: ";
  for (fit = overlay_arr.faces_begin(); fit != overlay_arr.faces_end(); ++fit) {
    if (! fit->data())
      continue;

    curr = fit->outer_ccb();
    std::cout << curr->source()->point();
    do {
      std::cout << " --> " << curr->target()->point();
      ++curr;
    } while (curr != fit->outer_ccb());
    std::cout << std::endl;
  }

  return 0;
}
