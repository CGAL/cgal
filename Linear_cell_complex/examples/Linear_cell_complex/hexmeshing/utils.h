#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/draw_linear_cell_complex.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT       FT;
typedef Kernel::Point_3  Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3> LCC;
typedef typename LCC::Dart_handle Dart_handle;
typedef typename LCC::Vertex_attribute_handle Vertex_handle;
typedef typename LCC::size_type   size_type;

void mark_face(LCC& lcc, Dart_handle dart, size_type mark){
  auto it = lcc.darts_of_orbit<1>(dart).begin();
  auto end = lcc.darts_of_orbit<1>(dart).end();

  for (; it != end; it++){
    if (!lcc.is_whole_cell_marked<0>(it, mark)){
      lcc.mark_cell<0>(it, mark);
    }
  }
}

void mark_vol(LCC& lcc, Dart_handle dart, size_type mark){
  auto it = lcc.darts_of_orbit<1,2>(dart).begin();
  auto end = lcc.darts_of_orbit<1,2>(dart).end();

  for (; it != end; it++){
    if (!lcc.is_whole_cell_marked<0>(it, mark)){
      lcc.mark_cell<0>(it, mark);
    }
  }
}

void mark_edge(LCC& lcc, Dart_handle dart, size_type mark){
  if (!lcc.is_whole_cell_marked<0>(dart, mark)){
    lcc.mark_cell<0>(dart, mark);
  }

  auto ext = lcc.other_extremity(dart);
  if (!lcc.is_whole_cell_marked<0>(ext, mark)){
    lcc.mark_cell<0>(ext, mark);
  }
}

void render(LCC &lcc, size_type marked_cell_0, size_type marked_cell_1)
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(lcc, buffer);

  // Hacking a bit here, we are going to draw colored points for marked vertices and edges just for debugging purposes
  // All of this is temporary of course and very ugly

  (const_cast<CGAL::Buffer_for_vao &>(buffer.get_buffer_for_mono_points())).clear();
  (const_cast<CGAL::Buffer_for_vao &>(buffer.get_buffer_for_mono_segments())).clear();

  // Marked vertices
  for (auto it = lcc.one_dart_per_cell<0>().begin(), end = lcc.one_dart_per_cell<0>().end(); it != end; it++)
  {

    if (lcc.is_marked(it, marked_cell_0))
    {
      buffer.add_point(
          lcc.point(it),
          CGAL::Color(0, 255, 0));
      continue;
    }

    buffer.add_point(lcc.point(it));
  }

  // Marked edges
  for (auto it = lcc.one_dart_per_cell<1>().begin(), end = lcc.one_dart_per_cell<1>().end(); it != end; it++){
    // lcc.is_whole_cell_marked<1>()
    if (lcc.is_whole_cell_marked<1>(it, marked_cell_1)){
      buffer.add_segment(lcc.point(it), lcc.point(lcc.other_extremity(it)), CGAL::Color(255, 0, 0));
      continue;
    }

    buffer.add_segment(lcc.point(it), lcc.point(lcc.other_extremity(it)));
  }


  draw_graphics_scene(buffer);
}

