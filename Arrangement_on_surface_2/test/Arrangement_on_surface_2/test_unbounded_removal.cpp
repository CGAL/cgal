// Testing the removal function, covering all possible scenarios.

#include <vector>

#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/IO/Arr_iostream.h>

typedef CGAL::Gmpq                                      Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>               Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Segment_2                             Segment_2;
typedef Traits_2::Ray_2                                 Ray_2;
typedef Traits_2::Line_2                                Line_2;
typedef Traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;
typedef Arrangement_2::Halfedge_handle                  Halfedge_handle;

int main()
{
  // Construct an arrangement of five linear objects.
  std::vector<X_monotone_curve_2>  curves;
  curves.push_back(Ray_2(Point_2 (0, 0), Point_2(0, 1)));
  curves.push_back(Ray_2(Point_2 (0, 0), Point_2(1, 0)));
  // curves.push_back(Ray_2(Point_2 (0, 0), Point_2(0, -1)));
  // curves.push_back(Ray_2(Point_2 (0, 0), Point_2(-1, 0)));
  // curves.push_back(Ray_2(Point_2 (0, 0), Point_2(1, 1)));
  // curves.push_back(Ray_2(Point_2 (0, 0), Point_2(1, -1)));
  // curves.push_back(Ray_2(Point_2 (0, 0), Point_2(-1, 1)));
  // curves.push_back(Ray_2(Point_2 (0, 0), Point_2(-1, -1)));

  Arrangement_2 arr;
  std::vector<Halfedge_handle> hhs(curves.size());
  for (unsigned int i = 0; i < curves.size(); i++)
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
  for (unsigned int i = 0; i < hhs.size(); i++) {
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

  // Check the validity more thoroughly.
  valid = is_valid(arr);
  std::cout << "Arrangement is "
              << (valid ? "valid." : "NOT valid!") << std::endl;

  return (0);
}

