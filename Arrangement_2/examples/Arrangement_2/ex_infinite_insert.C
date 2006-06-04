//! \file examples/Arrangement_2/ex_infinite_non_intersecting.C
// Constructing an arrangement of unbounded linear objects using the insertion
// function for non-intersecting curves.

#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>

typedef CGAL::Gmpq                                    Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>             Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Segment_2                           Segment_2;
typedef Traits_2::Ray_2                               Ray_2;
typedef Traits_2::Line_2                              Line_2;
typedef Traits_2::X_monotone_curve_2                  X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef CGAL::Arr_naive_point_location<Arrangement_2> Naive_pl;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Walk_pl;

int main ()
{
  // Construct 3 lines and 2 rays.
  const int          n_curves = 5;
  X_monotone_curve_2 curves[5];

  curves[0] = Line_2 (Point_2 (0, 0), Point_2 (2, 1));
  curves[1] = Line_2 (Point_2 (0, 0), Point_2 (-2, 1));
  curves[2] = Ray_2 (Point_2 (0, 2), Point_2 (-2, 0));
  curves[3] = Ray_2 (Point_2 (0, 2), Point_2 (2, 0));
  curves[4] = Line_2 (Point_2 (-1, 4), Point_2 (1, 4));

  // Construct the arrangement by inserting the curves incermentally.
  Arrangement_2      arr;
  //  Naive_pl           naive_pl (arr);
  Walk_pl            walk_pl (arr);
  int                k;

  for (k = 0; k < n_curves; k++)
  {
    std::cout << "Inserting curve no. " << k + 1 << std::endl;
    insert_x_monotone_curve (arr, curves[k], walk_pl);
  }

  // Print out the size of the resulting arrangement.
  Arrangement_2::Face_const_iterator   fit;
  unsigned int                         n_unb_faces = 0;

  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
  {
    if (fit->is_unbounded())
      n_unb_faces++;
  }

  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << " (" << arr.number_of_vertices_at_infinity()
            << " at infinity)"
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces()
            << " (" << n_unb_faces << " unbounded)" << std::endl;

  // Print the vertices.
  Arrangement_2::Vertex_const_iterator   vit;

  std::cout << "The vertices:" << std::endl;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
    std::cout << "  (" << vit->point() << ")" << std::endl;

  return (0);
}
