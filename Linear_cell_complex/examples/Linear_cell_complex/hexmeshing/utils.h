#pragma once

#include <CGAL/Graphics_scene.h>
#include <CGAL/Handle.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Union_find.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <unordered_set>

template <typename T>
// TODO make union_find arg const
std::vector<typename CGAL::Union_find<T>::handle> get_partitions(CGAL::Union_find<T>& union_find){
  using Handle = typename CGAL::Union_find<T>::handle;
  std::vector<Handle> result;

  if (union_find.number_of_sets() == 0) return result;

  for (auto it = union_find.begin(), end = union_find.end(); it != end; it++ ){
    if (it.ptr()->up == nullptr) {
      result.push_back(it);

      if (result.size() == union_find.number_of_sets())
        return result;
    }
  }

  CGAL_assertion_msg(false, "get_partions(Union_find<T>) function did not find all sets. Bug?");
  return result;
};
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
using LCCSceneOptions = CGAL::Graphics_scene_options<LCC,
                          typename LCC::Dart_const_handle,
                          typename LCC::Dart_const_handle,
                          typename LCC::Dart_const_handle,
                          typename LCC::Dart_const_handle>;

template <typename LCC>
CGAL::IO::Color rand_color_from_dart(const LCC& lcc, typename LCC::Dart_const_handle dart){
  CGAL::Random random((unsigned int)(lcc.darts().index(dart)));
  return CGAL::get_random_color(random);
}

template <typename LCC>
void render(LCC &lcc, CGAL::Graphics_scene& buffer, typename LCC::size_type marked_cell_0, typename LCC::size_type marked_cell_1){
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


template <typename LCC>
void render(LCC &lcc, typename LCC::size_type marked_cell_0, typename LCC::size_type marked_cell_1)
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(lcc, buffer);
  render(lcc, buffer, marked_cell_0, marked_cell_1);
}


template <typename LCC>
void render(LCC &lcc, LCCSceneOptions<LCC> options, typename LCC::size_type marked_cell_0, typename LCC::size_type marked_cell_1)
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(lcc, buffer, options);
  render(lcc, buffer, marked_cell_0, marked_cell_1);
}
