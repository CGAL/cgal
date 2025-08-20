
#ifndef HEXMESHING_THREADS_UTILS_H
#define HEXMESHING_THREADS_UTILS_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>

namespace CGAL::Hexmeshing {
  template <typename HexData>
  void thread_number_vertex_in_edge(HexData& hdata,
    Dart_handle node, Dart_handle extremity0, Dart_handle extremity1){}

  template <typename HexData>
  void thread_number_vertex_in_1t_face(HexData& hdata, Dart_handle node) {}

  template <typename HexData>
  void thread_number_vertex_in_1t_vol(HexData& hdata, Dart_handle v_signature_start) {}

  template <typename HexData>
  void thread_join_3_template_vertex__pair(HexData& hdata, Dart_handle edge) {}

  template <typename HexData>
  void thread_join_3_template_vertex__pairpair(HexData& hdata, Dart_handle edge) {}

  template <typename HexData>
  void thread_communicate_marked_nodes(HexData&, RefinementData&, size_type) {}

  template <typename HexData>
  void thread_communicate_cells_id_and_3t(HexData&, RefinementData&){}

  template <typename HexData>
  void thread_remove_ghosts(HexData& hdata) {}
}





#endif