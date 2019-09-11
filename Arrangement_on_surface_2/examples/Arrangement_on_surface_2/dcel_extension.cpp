//! \file examples/Arrangement_on_surface_2/dcel_extension.cpp
// Extending all DCEL records (vertices, edges and faces).

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>

enum Color {BLUE, RED, WHITE};

typedef CGAL::Cartesian<CGAL::Exact_rational>                   Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef Traits_2::Point_2                                       Point_2;
typedef Traits_2::X_monotone_curve_2                            Segment_2;
typedef CGAL::Arr_extended_dcel<Traits_2,Color, bool, int>      Dcel;
typedef CGAL::Arrangement_2<Traits_2, Dcel>                     Arrangement_2;

int main ()
{
  // Construct the arrangement containing two intersecting triangles.
  Arrangement_2          arr;

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

  // Go over all arrangement vertices and set their colors according to our
  // coloring convention.
  Arrangement_2::Vertex_iterator            vit;
  std::size_t                               degree;

  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
  {
    degree = vit->degree();
    if (degree == 0)
      vit->set_data (BLUE);       // Isolated vertex.
    else if (degree <= 2)
      vit->set_data (RED);        // Vertex represents an endpoint.
    else
      vit->set_data (WHITE);      // Vertex represents an intersection point.
  }

  // Go over all arrangement edges and set their flags.
  Arrangement_2::Edge_iterator              eit;
  bool                                      flag;

  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    // Check if the halfedge has the same direction as its associated
    // segment. Note that its twin always has an opposite direction.
    flag = (eit->source()->point() == eit->curve().source());
    eit->set_data (flag);
    eit->twin()->set_data (!flag);
  }

  // For each arrangement face, print the outer boundary and its size.
  Arrangement_2::Face_iterator              fit;
  Arrangement_2::Ccb_halfedge_circulator    curr;
  int                                       boundary_size;

  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    boundary_size = 0;
    if (! fit->is_unbounded()) {
      curr = fit->outer_ccb();
      do {
        ++boundary_size;
        ++curr;
      } while (curr != fit->outer_ccb());
    }
    fit->set_data (boundary_size);
  }

  // Copy the arrangement and print the vertices.
  Arrangement_2    arr2 = arr;

  std::cout << "The arrangement vertices:" << std::endl;
  for (vit = arr2.vertices_begin(); vit != arr2.vertices_end(); ++vit) {
    std::cout << '(' << vit->point() << ") - ";
    switch (vit->data()) {
     case BLUE  : std::cout << "BLUE."  << std::endl; break;
     case RED   : std::cout << "RED."   << std::endl; break;
     case WHITE : std::cout << "WHITE." << std::endl; break;
     default: break;
    }
  }

  return 0;
}
