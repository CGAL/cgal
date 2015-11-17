//! \file examples/Arrangement_on_surface_2/face_extension.cpp
// Extending the arrangement-face records.

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_observer.h>

typedef CGAL::Cartesian<CGAL::Exact_rational>          Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>             Traits_2;
typedef Traits_2::Point_2                              Point_2;
typedef Traits_2::X_monotone_curve_2                   Segment_2;
typedef CGAL::Arr_face_extended_dcel<Traits_2, int>    Dcel;
typedef CGAL::Arrangement_2<Traits_2, Dcel>            Arrangement_2;

// An arrangement observer, used to receive notifications of face splits and
// to update the indices of the newly created faces.
class Face_index_observer : public CGAL::Arr_observer<Arrangement_2>
{
private:
  int     n_faces;          // The current number of faces.

public:

  Face_index_observer (Arrangement_2& arr) :
    CGAL::Arr_observer<Arrangement_2> (arr),
    n_faces (0)
  {
    CGAL_precondition (arr.is_empty());

    arr.unbounded_face()->set_data (0);
    n_faces++;
  }

  virtual void after_split_face (Face_handle /* old_face */,
                                 Face_handle new_face, bool )
  {
    // Assign index to the new face.
    new_face->set_data (n_faces);
    n_faces++;
  }
};

int main ()
{
  // Construct the arrangement containing two intersecting triangles.
  Arrangement_2          arr;
  Face_index_observer    obs (arr);

  Segment_2      s1 (Point_2(4, 1), Point_2(7, 6));
  Segment_2      s2 (Point_2(1, 6), Point_2(7, 6));
  Segment_2      s3 (Point_2(4, 1), Point_2(1, 6));
  Segment_2      s4 (Point_2(1, 3), Point_2(7, 3));
  Segment_2      s5 (Point_2(1, 3), Point_2(4, 8));
  Segment_2      s6 (Point_2(4, 8), Point_2(7, 3));

  insert_non_intersecting_curve (arr, s1);
  insert_non_intersecting_curve (arr, s2);
  insert_non_intersecting_curve (arr, s3);
  insert (arr, s4);
  insert (arr, s5);
  insert (arr, s6);

  // Go over all arrangement faces and print the index of each face and it
  // outer boundary. The face index is stored in its data field in our case.
  Arrangement_2::Face_const_iterator            fit;
  Arrangement_2::Ccb_halfedge_const_circulator  curr;

  std::cout << arr.number_of_faces() << " faces:" << std::endl;
  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    std::cout << "Face no. " << fit->data() << ": ";
    if (fit->is_unbounded())
      std::cout << "Unbounded." << std::endl;
    else {
      curr = fit->outer_ccb();
      std::cout << curr->source()->point();
      do {
        std::cout << " --> " << curr->target()->point();
        ++curr;
      } while (curr != fit->outer_ccb());
      std::cout << std::endl;
    }
  }

  return 0;
}
