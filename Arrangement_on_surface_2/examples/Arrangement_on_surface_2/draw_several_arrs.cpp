//! \file examples/Arrangement_on_surface_2/draw_several_arrs.cpp
// Drawing two arrangements in the same window.

#include <cstdlib>
#include <vector>

#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Basic_viewer.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/draw_arrangement_2.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Number_type = Kernel::FT;

using Segment_traits = CGAL::Arr_segment_traits_2<Kernel>;
using Segment_arrangement = CGAL::Arrangement_2<Segment_traits>;

using Circle_traits = CGAL::Arr_circle_segment_traits_2<Kernel>;
using Circle_arrangement = CGAL::Arrangement_2<Circle_traits>;

int main() {
  // An arrangement of segments that forms a grid.
  Segment_arrangement seg_arr;
  using Point = Kernel::Point_2;
  using Segment = Kernel::Segment_2;
  std::vector<Segment> segments;
  for (int i = 0; i <= 4; ++i) {
    segments.emplace_back(Point(0, 2 * i), Point(8, 2 * i));
    segments.emplace_back(Point(2 * i, 0), Point(2 * i, 8));
  }
  CGAL::insert(seg_arr, segments.begin(), segments.end());

  // An arrangement of circles.
  Circle_arrangement circ_arr;
  using Circle = Kernel::Circle_2;
  using Curve = Circle_traits::Curve_2;
  CGAL::insert(circ_arr, Curve(Circle(Point(4, 4), Number_type(16))));
  CGAL::insert(circ_arr, Curve(Circle(Point(2, 6), Number_type(4))));
  CGAL::insert(circ_arr, Curve(Circle(Point(6, 2), Number_type(4))));

  // Draw the circles in blue on top of the grid; do not draw their faces, so
  // that they do not hide the faces of the grid.
  using Circle_gso = CGAL::Graphics_scene_options<Circle_arrangement, Circle_arrangement::Vertex_const_handle,
                                                  Circle_arrangement::Halfedge_const_handle,
                                                  Circle_arrangement::Face_const_handle>;
  Circle_gso circ_gso;
  circ_gso.disable_faces();
  circ_gso.colored_edge = [](const Circle_arrangement&, Circle_arrangement::Halfedge_const_handle) { return true; };
  circ_gso.edge_color = [](const Circle_arrangement&, Circle_arrangement::Halfedge_const_handle)
  { return CGAL::IO::Color(0, 0, 255); };

  CGAL::Graphics_scene scene;
  CGAL::add_to_graphics_scene(seg_arr, scene);
  CGAL::add_to_graphics_scene(circ_arr, scene, circ_gso);
  CGAL::draw_graphics_scene(scene, "a grid and circles");
  return EXIT_SUCCESS;
}
