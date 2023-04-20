// Copyright (c) 2021  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau
//

#ifndef CGAL_TEST_T23_MOVE_SEMANTIC_C
#define CGAL_TEST_T23_MOVE_SEMANTIC_C

#include <cassert>

namespace CGAL {
  namespace Testsuite {
    namespace Triangulation_23 {
      template <typename Tr>
      void test_move_semantic(Tr source_tr) {
        const auto dimension = source_tr.dimension();
        const auto nb_of_vertices = source_tr.number_of_vertices();
        auto check_triangulation_validity = [&](const Tr& tr) {
          assert(tr.is_valid());
          assert(tr.number_of_vertices() == nb_of_vertices);
          assert(tr.dimension() == dimension);
        };
        auto check_moved_from_triangulation = [](const Tr& tr_copy) {
          assert(tr_copy.dimension() == -2);
          assert(tr_copy.number_of_vertices() + 1 == 0);
        };
        auto check_empty_triangulation = [](const Tr& tr_copy2) {
          assert(tr_copy2.dimension() == -1);
          assert(tr_copy2.number_of_vertices() == 0);
        };
        // move constructor
        {
          Tr tr_copy(source_tr);
          check_triangulation_validity(tr_copy);

          Tr tr_move_constructed(std::move(tr_copy));
          check_triangulation_validity(tr_move_constructed);
          check_moved_from_triangulation(tr_copy);

          Tr tr_copy2(source_tr);
          Tr tr_move_constructed2(std::move(tr_copy2));
          check_moved_from_triangulation(tr_copy2);
          tr_copy2.clear();
          check_empty_triangulation(tr_copy2);

          Tr tr_copy3(source_tr);
          Tr tr_move_constructed3(std::move(tr_copy3));
          check_moved_from_triangulation(tr_copy3);
          tr_copy3 = source_tr;
          check_triangulation_validity(tr_copy3);
        }
        // move-assignment
        {
          Tr tr_copy4(source_tr);
          Tr tr_move_assigned;
          tr_move_assigned = std::move(tr_copy4);
          check_triangulation_validity(tr_move_assigned);
          check_moved_from_triangulation(tr_copy4);
        }
      };
    }
  }
}

#endif // CGAL_TEST_T23_MOVE_SEMANTIC_C
