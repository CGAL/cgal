// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Soichiro Yamazaki <soichiro19998@gmail.com>, Théo Bénard <benard320@gmail.com>
//

#ifndef CGAL_HEXMESHING_THREADS_UTILS_H
#define CGAL_HEXMESHING_THREADS_UTILS_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>

namespace CGAL::internal::Hexmeshing
{
  template <typename HexData>
  void thread_number_vertex_in_edge(HexData& hdata,
    Dart_descriptor node, Dart_descriptor extremity0, Dart_descriptor extremity1){}

  template <typename HexData>
  void thread_number_vertex_in_1t_face(HexData& hdata, Dart_descriptor node) {}

  template <typename HexData>
  void thread_number_vertex_in_1t_vol(HexData& hdata, Dart_descriptor v_signature_start) {}

  template <typename HexData>
  void thread_join_3_template_vertex__pair(HexData& hdata, Dart_descriptor edge) {}

  template <typename HexData>
  void thread_join_3_template_vertex__pairpair(HexData& hdata, Dart_descriptor edge) {}

  template <typename HexData>
  void thread_communicate_marked_nodes(HexData&, RefinementData&, size_type) {}

  template <typename HexData>
  void thread_communicate_cells_id_and_3t(HexData&, RefinementData&){}

  template <typename HexData>
  void thread_remove_ghosts(HexData& hdata) {}
}

#endif
