// file: examples/Arrangement_2/example19.C

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
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
typedef CGAL::Arr_naive_point_location<Arrangement_2> Point_location;

// An arrangement observer, used to receive notifications of face splits and
// face mergers.
class My_observer : public CGAL::Arr_observer<Arrangement_2>
{
public:

  My_observer (Arrangement_2& arr) :
    CGAL::Arr_observer<Arrangement_2> (arr)
  {}

  virtual void before_split_face (Face_handle,
                                  Halfedge_handle e)
  {
    std::cout << "-> The insertion of :  [ " << e->curve()
              << " ]  causes a face to split." << std::endl;
  }

  virtual void before_merge_face (Face_handle,
                                  Face_handle,
                                  Halfedge_handle e)
  {
    std::cout << "-> The removal of :  [ " << e->curve()
              << " ]  causes two faces to merge." << std::endl;
  }

};

int main ()
{
  // Construct the arrangement containing one diamond-shaped face.
  Arrangement_2  arr;
  My_observer    obs (arr);
  Point_location pl (arr);

  Segment_2      s1 (Point_2(-1, 0), Point_2(0, 1));
  Segment_2      s2 (Point_2(0, 1), Point_2(1, 0));
  Segment_2      s3 (Point_2(1, 0), Point_2(0, -1));
  Segment_2      s4 (Point_2(0, -1), Point_2(-1, 0));

  insert_non_intersecting (arr, s1, pl);
  insert_non_intersecting (arr, s2, pl);
  insert_non_intersecting (arr, s3, pl);
  insert_non_intersecting (arr, s4, pl);

  // Insert a vertical segment dividing the diamond into two, and a
  // a horizontal segment dividing the diamond into four:
  Segment_2      s_vert (Point_2(0, -1), Point_2(0, 1));
  Arrangement_2::Halfedge_handle
                 e_vert = insert_non_intersecting (arr, s_vert, pl);

  Segment_2      s_horiz (Point_2(-1, 0), Point_2(1, 0));

  insert (arr, s_horiz, pl);

  std::cout << "V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;

  // Now remove a portion of the vertical segment.
  remove_edge (arr, e_vert);
 
  std::cout << "V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;

  return (0);
}

