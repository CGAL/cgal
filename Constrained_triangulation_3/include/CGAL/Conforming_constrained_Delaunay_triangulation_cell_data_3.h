// Copyright (c) 2019-2024  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_CELL_DATA_3_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_CELL_DATA_3_H

#include <CGAL/license/Constrained_triangulation_3.h>

#include <CGAL/Constrained_triangulation_3/internal/config.h>

#include <array>
#include <bitset>

namespace CGAL {
enum class CDT_3_cell_marker {
  CLEAR = 0,
  IN_REGION = 1,
  VISITED = 1,
  ON_REGION_BOUNDARY = 2,
  nb_of_markers
};

/*!
 * @ingroup PkgConstrainedTriangulation3Classes
 * @brief Internal per-cell data for \cgal 3D conforming constrained Delaunay triangulations
 *
 * This class is an internal detail of the implementation of \cgal 3D conforming constrained Delaunay triangulations.
 *
 * Any model of the `ConformingConstrainedDelaunayTriangulationCellBase_3` concept must include one object of this type
 * as a non-static data member.
 */
class Conforming_constrained_Delaunay_triangulation_cell_data_3 {
  /// @cond SKIP_IN_MANUAL
  template <typename Tr> friend class Conforming_constrained_Delaunay_triangulation_3_impl;
  /// @endcond

  std::array<CDT_3_signed_index, 4> face_id = { -1, -1, -1, -1 };
  std::array<void*, 4> facet_2d = {nullptr, nullptr, nullptr, nullptr};
  std::bitset<static_cast<unsigned>(CDT_3_cell_marker::nb_of_markers)> markers;

  bool is_marked() const { return markers.any(); }
  bool is_marked(CDT_3_cell_marker m) const { return markers.test(static_cast<unsigned>(m)); }
  void set_mark(CDT_3_cell_marker m) { markers.set(static_cast<unsigned>(m)); }
  void clear_mark(CDT_3_cell_marker m) { markers.reset(static_cast<unsigned>(m)); }
  void clear_marks() { markers.reset(); }

  template <typename Facet_handle>
  void set_facet_constraint(int i, CDT_3_signed_index face_id,
                            Facet_handle facet_2d)
  {
    this->face_id[unsigned(i)] = face_id;
    this->facet_2d[unsigned(i)] = static_cast<void*>(facet_2d == Facet_handle{} ?  nullptr : std::addressof(*facet_2d));
  }

  template <typename CDT_2>
  auto face_2 (const CDT_2& cdt, int i) const {
    using Face = typename CDT_2::Face;
    auto ptr = static_cast<Face*>(facet_2d[unsigned(i)]);
    return cdt.tds().faces().iterator_to(*ptr);
  }
public:
  /// @brief Returns if the i-th facet of the cell is constrained.
  bool is_facet_constrained(int i) const { return face_id[unsigned(i)] >= 0; }

  /// @brief Returns the index of the constraint that constrains the
  /// i-th facet of the cell.
  /// @pre `is_facet_constrained(i)`
  CDT_3_signed_index face_constraint_index(int i) const {
    return face_id[unsigned(i)];
  }
  /// @brief Set the index of the constraint that constrains the i-th facet of the cell.
  void set_face_constraint_index(int i, CDT_3_signed_index index) {
    face_id[unsigned(i)] = index;
  }
};

} // namespace CGAL

#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_CELL_DATA_3_H
