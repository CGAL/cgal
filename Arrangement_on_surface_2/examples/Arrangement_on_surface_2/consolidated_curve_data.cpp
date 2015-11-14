//! \file examples/Arrangement_on_surface_2/consolidated_curve_data.cpp
// Associating a color attribute with segments using the consolidated
// curve-data traits.

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_landmarks_point_location.h>

enum Segment_color {
  RED,
  BLUE
};

typedef CGAL::Cartesian<CGAL::Exact_rational>             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                Segment_traits_2;
typedef Segment_traits_2::Curve_2                         Segment_2;
typedef CGAL::Arr_consolidated_curve_data_traits_2
                   <Segment_traits_2, Segment_color>      Traits_2;
typedef Traits_2::Point_2                                 Point_2;
typedef Traits_2::Curve_2                                 Colored_segment_2;
typedef CGAL::Arrangement_2<Traits_2>                     Arrangement_2;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2> Landmarks_pl;

int main ()
{
  // Construct an arrangement containing three RED line segments.
  Arrangement_2     arr;
  Landmarks_pl      pl (arr);

  Segment_2         s1 (Point_2(-1, -1), Point_2(1, 3));
  Segment_2         s2 (Point_2(2, 0), Point_2(3, 3));
  Segment_2         s3 (Point_2(0, 3), Point_2(2, 5));

  insert (arr, Colored_segment_2 (s1, RED), pl);
  insert (arr, Colored_segment_2 (s2, RED), pl);
  insert (arr, Colored_segment_2 (s3, RED), pl);

  // Insert three BLUE line segments.
  Segment_2         s4 (Point_2(-1, 3), Point_2(4, 1));
  Segment_2         s5 (Point_2(-1, 0), Point_2(4, 1));
  Segment_2         s6 (Point_2(-2, 1), Point_2(1, 4));

  insert (arr, Colored_segment_2 (s4, BLUE), pl);
  insert (arr, Colored_segment_2 (s5, BLUE), pl);
  insert (arr, Colored_segment_2 (s6, BLUE), pl);

  // Go over all vertices and print just the ones corresponding to intersection
  // points between RED segments and BLUE segments. Note that we skip endpoints
  // of overlapping sections.
  Arrangement_2::Vertex_const_iterator   vit;
  Segment_color                          color;

  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    // Go over the incident halfedges of the current vertex and examine their
    // colors.
    bool       has_red = false;
    bool       has_blue = false;

    Arrangement_2::Halfedge_around_vertex_const_circulator  eit, first;

    eit = first = vit->incident_halfedges();
    do {
      // Get the color of the current half-edge.
      if (eit->curve().data().size() == 1) {
        color = eit->curve().data().front();

        if (color == RED)
          has_red = true;
        else if (color == BLUE)
          has_blue = true;
      }

      ++eit;
    } while (eit != first);

    // Print the vertex only if incident RED and BLUE edges were found.
    if (has_red && has_blue)
    {
      std::cout << "Red-blue intersection at (" << vit->point() << ")"
                << std::endl;
    }
  }

  // Locate the edges that correspond to a red-blue overlap.
  Arrangement_2::Edge_iterator   eit;

  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit)
  {
    // Go over the incident edges of the current vertex and examine their
    // colors.
    bool       has_red = false;
    bool       has_blue = false;

    Traits_2::Data_container::const_iterator       dit;

    for (dit = eit->curve().data().begin(); dit != eit->curve().data().end();
         ++dit)
    {
      if (*dit == RED)
        has_red = true;
      else if (*dit == BLUE)
        has_blue = true;
    }

    // Print the edge only if it corresponds to a red-blue overlap.
    if (has_red && has_blue)
      std::cout << "Red-blue overlap at [" << eit->curve() << "]"  << std::endl;
  }
  return 0;
}
