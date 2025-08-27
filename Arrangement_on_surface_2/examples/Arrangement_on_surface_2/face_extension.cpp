//! \file examples/Arrangement_on_surface_2/face_extension.cpp
// Extending the arrangement-face records.

#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_observer.h>

#include "arr_exact_construction_segments.h"

using Dcel = CGAL::Arr_face_extended_dcel<Traits, size_t>;
using Ex_arrangement = CGAL::Arrangement_2<Traits, Dcel>;

// An arrangement observer, used to receive notifications of face splits and
// to update the indices of the newly created faces.
class Face_index_observer : public CGAL::Arr_observer<Ex_arrangement> {
private:
  size_t n_faces;                       // the current number of faces

public:
  Face_index_observer(Ex_arrangement& arr) :
    CGAL::Arr_observer<Ex_arrangement>(arr),
    n_faces(0)
  {
    CGAL_precondition(arr.is_empty());
    arr.unbounded_face()->set_data (0);
  }

  virtual void after_split_face(Face_handle, Face_handle new_face, bool)
  {
    new_face->set_data(++n_faces);        // assign index to the new face
  }
};

int main() {
  // Construct the arrangement containing two intersecting triangles.
  Ex_arrangement arr;
  Face_index_observer obs(arr);
  insert_non_intersecting_curve(arr, Segment(Point(4, 1), Point(7, 6)));
  insert_non_intersecting_curve(arr, Segment(Point(1, 6), Point(7, 6)));
  insert_non_intersecting_curve(arr, Segment(Point(4, 1), Point(1, 6)));
  insert(arr, Segment(Point(1, 3), Point(7, 3)));
  insert(arr, Segment(Point(1, 3), Point(4, 8)));
  insert(arr, Segment(Point(4, 8), Point(7, 3)));

  // Go over all arrangement faces and print the index of each face and its
  // outer boundary. The face index is stored in the data field.
  std::cout << arr.number_of_faces() << " faces:\n";
  for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    std::cout << "Face no. " << fit->data() << ": ";
    if (fit->is_unbounded()) std::cout << "Unbounded.\n";
    else {
      auto curr = fit->outer_ccb();
      std::cout << curr->source()->point();
      do std::cout << " --> " << curr->target()->point();
      while (++curr != fit->outer_ccb());
      std::cout << std::endl;
    }
  }

  return 0;
}
