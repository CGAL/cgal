//! \file examples/Arrangement_on_surface_2/overlay_unbounded.cpp
// A face overlay of two arrangements with unbounded faces.

#include <string>
#include <boost/lexical_cast.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>

// Define a functor for creating a label from a character and an integer.
struct Overlay_label
{
  std::string operator() (char c, int i) const
  {
    return boost::lexical_cast<std::string>(c) +
      boost::lexical_cast<std::string>(i);
  }
};

typedef CGAL::Exact_predicates_exact_constructions_kernel   Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>                   Traits_2;
typedef Traits_2::Point_2                                   Point_2;
typedef Traits_2::Segment_2                                 Segment_2;
typedef Traits_2::Ray_2                                     Ray_2;
typedef Traits_2::Line_2                                    Line_2;
typedef Traits_2::X_monotone_curve_2                        X_monotone_curves_2;
typedef CGAL::Arr_face_extended_dcel<Traits_2, char>        DcelA;
typedef CGAL::Arrangement_2<Traits_2, DcelA>                ArrangementA_2;
typedef CGAL::Arr_face_extended_dcel<Traits_2, int>         DcelB;
typedef CGAL::Arrangement_2<Traits_2, DcelB>                ArrangementB_2;
typedef CGAL::Arr_face_extended_dcel<Traits_2, std::string> DcelRes;
typedef CGAL::Arrangement_2<Traits_2, DcelRes>              ArrangementRes_2;
typedef CGAL::Arr_face_overlay_traits<ArrangementA_2,
                                      ArrangementB_2,
                                      ArrangementRes_2,
                                      Overlay_label>        Overlay_traits;

int main ()
{
  // Construct the first arrangement, induced by two line y = x and y = -x.
  ArrangementA_2          arr1;

  insert (arr1, Line_2 (Point_2(0, 0), Point_2(1, 1)));
  insert (arr1, Line_2 (Point_2(0, 0), Point_2(1, -1)));

  // Label the four (unbounded) face of the arrangement as 'A' to 'D'.
  // We do so by traversing the incident faces to the halfedges around the
  // single arrangement vertex (0, 0).
  CGAL_assertion (arr1.number_of_vertices() == 1);

  ArrangementA_2::Halfedge_around_vertex_circulator  first, curr;
  char                                               clabel = 'A';

  curr = first = arr1.vertices_begin()->incident_halfedges();
  do {
    curr->face()->set_data (clabel);
    ++clabel;
    ++curr;
  } while (curr != first);
  std::cout << "Done with arr1." << std::endl;

  // Construct the second arrangement, containing a single square-shaped face.
  ArrangementB_2          arr2;

  insert (arr2, Segment_2 (Point_2(-4, -4), Point_2(4, -4)));
  insert (arr2, Segment_2 (Point_2(4, -4), Point_2(4, 4)));
  insert (arr2, Segment_2 (Point_2(4, 4), Point_2(-4, 4)));
  insert (arr2, Segment_2 (Point_2(-4, 4), Point_2(-4, -4)));

  // Give the unbounded face the index 1, and the bounded face the index 2.
  CGAL_assertion (arr2.number_of_faces() == 2);

  ArrangementB_2::Face_iterator    fit;

  for (fit = arr2.faces_begin(); fit != arr2.faces_end(); ++fit)
    fit->set_data ((fit == arr2.unbounded_face()) ? 1 : 2);
  std::cout << "Done with arr2." << std::endl;

  // Compute the overlay of the two arrangements.
  ArrangementRes_2       overlay_arr;
  Overlay_traits         overlay_traits;

  overlay (arr1, arr2, overlay_arr, overlay_traits);

  // Go over the faces of the overlaid arrangement and their labels.
  ArrangementRes_2::Face_iterator  res_fit;

  std::cout << "The overlay faces are: ";
  for (res_fit = overlay_arr.faces_begin();
       res_fit != overlay_arr.faces_end(); ++res_fit)
  {
    std::cout << res_fit->data() << " ("
              << (res_fit->is_unbounded() ? "unbounded" : "bounded")
              << ")." << std::endl;
  }

  return 0;
}
