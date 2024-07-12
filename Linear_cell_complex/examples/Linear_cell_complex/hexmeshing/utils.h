#pragma once

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/draw_linear_cell_complex.h>

/**
 * 
 * Will be later removed
 * Used for debugging stuff
 * 
 */

template <typename LCC>
void mark_face(LCC& lcc, typename LCC::Dart_handle dart, typename LCC::size_type mark){
  auto it = lcc.template darts_of_orbit<1>(dart).begin();
  auto end = lcc.template darts_of_orbit<1>(dart).end();

  for (; it != end; it++){
    if (!lcc.template is_whole_cell_marked<0>(it, mark)){
      lcc.template mark_cell<0>(it, mark);
    }
  }
}

template <typename LCC>
void mark_vol(LCC& lcc, typename LCC::Dart_handle dart, typename LCC::size_type mark){
  auto it = lcc.template darts_of_orbit<1,2>(dart).begin();
  auto end = lcc.template darts_of_orbit<1,2>(dart).end();

  for (; it != end; it++){
    if (!lcc.template is_whole_cell_marked<0>(it, mark)){
      lcc.template mark_cell<0>(it, mark);
    }
  }
}

template <typename LCC>
void mark_edge(LCC& lcc, typename LCC::Dart_handle dart, typename LCC::size_type mark){
  if (!lcc.template is_whole_cell_marked<0>(dart, mark)){
    lcc.template mark_cell<0>(dart, mark);
  }

  auto ext = lcc.other_extremity(dart);
  if (!lcc.template is_whole_cell_marked<0>(ext, mark)){
    lcc.template mark_cell<0>(ext, mark);
  }
}

template <typename LCC>
void render(LCC &lcc, typename LCC::size_type marked_cell_0, typename LCC::size_type marked_cell_1)
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(lcc, buffer);

  // Hacking a bit here, we are going to draw colored points for marked vertices and edges just for debugging purposes
  // All of this is temporary of course and very ugly

  (const_cast<CGAL::Buffer_for_vao &>(buffer.get_buffer_for_mono_points())).clear();
  (const_cast<CGAL::Buffer_for_vao &>(buffer.get_buffer_for_mono_segments())).clear();

  // Marked vertices
  for (auto it = lcc.template one_dart_per_cell<0>().begin(), end = lcc.template one_dart_per_cell<0>().end(); it != end; it++)
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
  for (auto it = lcc.template one_dart_per_cell<1>().begin(), 
      end = lcc.template one_dart_per_cell<1>().end(); 
      it != end; 
      it++){
    // lcc.is_whole_cell_marked<1>()
    if (lcc.template is_whole_cell_marked<1>(it, marked_cell_1)){
      buffer.add_segment(lcc.point(it), lcc.point(lcc.other_extremity(it)), CGAL::Color(0, 0, 255));
      continue;
    }

    buffer.add_segment(lcc.point(it), lcc.point(lcc.other_extremity(it)));
  }

  draw_graphics_scene(buffer);
}
