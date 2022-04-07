// Testing the removal function, covering all possible scenarios.

#include <vector>
#include <iostream>

#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Exact_rational                            Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>               Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Segment_2                             Segment_2;
typedef Traits_2::Ray_2                                 Ray_2;
typedef Traits_2::Line_2                                Line_2;
typedef Traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;
typedef Arrangement_2::Vertex_handle                    Vertex_handle;
typedef Arrangement_2::Face_handle                      Face_handle;
typedef Arrangement_2::Halfedge_handle                  Halfedge_handle;

bool test(const Arrangement_2& arr)
{
  if (arr.number_of_vertices() != 0) {
    std::cerr << "(A) Number of vertices (" << arr.number_of_vertices()
              << ") not 0!" << std::endl;
    return false;
  }

  // Check the validity more thoroughly.
  if (! CGAL::is_valid(arr)) {
    std::cerr << "The arrangement is NOT valid!" << std::endl;
    return false;
  }

  return true;
}

bool test_ray(Arrangement_2& arr, Face_handle f)
{
  Segment_2 s1(Point_2(0, 0), Point_2(1, 0));
  Halfedge_handle eh1 =
    arr.insert_in_face_interior(X_monotone_curve_2(s1), f);
  Vertex_handle vh = eh1->target();
  Ray_2 ray(Point_2(1, 0), Point_2(2, 0));
  Halfedge_handle eh2 =
    arr.insert_from_left_vertex(X_monotone_curve_2(ray), vh);

  // Remove the edges
  arr.remove_edge(eh2);
  arr.remove_edge(eh1);

  if (!::test(arr)) return false;

  return true;
}

int main()
{
  Arrangement_2 arr;

  // Construct an arrangement of five linear objects.
  std::vector<X_monotone_curve_2>  curves;
  curves.push_back(Ray_2(Point_2(0, 0), Point_2(0, 1)));
  curves.push_back(Ray_2(Point_2(0, 0), Point_2(1, 0)));
  curves.push_back(Ray_2(Point_2(0, 0), Point_2(0, -1)));
  curves.push_back(Ray_2(Point_2(0, 0), Point_2(-1, 0)));
  curves.push_back(Ray_2(Point_2(0, 0), Point_2(1, 1)));
  curves.push_back(Ray_2(Point_2(0, 0), Point_2(1, -1)));
  curves.push_back(Ray_2(Point_2(0, 0), Point_2(-1, 1)));
  curves.push_back(Ray_2(Point_2(0, 0), Point_2(-1, -1)));

  std::vector<Halfedge_handle> hhs(curves.size());
  for (size_t i = 0; i < curves.size(); i++)
    hhs[i] = insert_non_intersecting_curve(arr, curves[i]);
  bool valid = arr.is_valid();
  std::cout << "Arrangement size:"
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;
  std::cout << "The arrangement is "
            << (valid ? "valid." : "NOT valid!") << std::endl;
  if (! valid) return (1);

  // Remove the edges.
  for (size_t i = 0; i < hhs.size(); i++) {
    arr.remove_edge(hhs[i]);
    bool valid = arr.is_valid();
    std::cout << "  Removed " << i+1 << " curve(s), arrangement is "
              << (valid ? "valid." : "NOT valid!") << std::endl;

    if (! valid) return (1);
  }

  std::cout << "Final arrangement size:"
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  if (!::test(arr)) return 1;

  /* Construct another arrangement of a segment connected to a ray.
   * First remove the ray, then remove the segment.
   *
   *  o---------------o
   *  |               |
   *  |               |
   *  |               |
   *  |       o--o----|
   *  |               |
   *  |               |
   *  |               |
   *  o---------------o
   */
  if (!test_ray(arr, arr.unbounded_face())) return 1;

  /* Construct another arrangement of a segment connected to a ray.
   * First remove the ray, then remove the segment.
   *
   *  o---------------o
   *  |               |
   *  |               |
   *  |               |
   *  |       o--o----|
   *  |               |
   *  o---------------o
   *  |               |
   *  o---------------o
   */
  Line_2 l1(Point_2(0, -1), Point_2(1, -1));
  Halfedge_handle eh1 =
    arr.insert_in_face_interior(X_monotone_curve_2(l1), arr.unbounded_face());

  if (!test_ray(arr, eh1->face())) return 1;
  arr.remove_edge(eh1);

  /* Construct another arrangement of a segment connected to a ray.
   * First remove the ray, then remove the segment.
   *
   *  o-----o---------o
   *  |     |         |
   *  |     |         |
   *  |     |         |
   *  |     | o--o----|
   *  |     |         |
   *  |     |         |
   *  |     |         |
   *  o-----o---------o
   */
  Line_2 l2(Point_2(-1, 0), Point_2(-1, 1));
  Halfedge_handle eh2 =
    arr.insert_in_face_interior(X_monotone_curve_2(l2), arr.unbounded_face());

  if (!test_ray(arr, eh2->twin()->face())) return 1;
  arr.remove_edge(eh2);

  return 0;
}
