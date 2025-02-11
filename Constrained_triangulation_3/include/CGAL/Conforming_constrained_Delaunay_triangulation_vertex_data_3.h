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

#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_VERTEX_DATA_3_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_VERTEX_DATA_3_H

#include <CGAL/license/Constrained_triangulation_3.h>

#include <CGAL/assertions.h>
#include <CGAL/Constrained_triangulation_3/internal/config.h>

#include <bitset>

namespace CGAL {

#ifdef DOXYGEN_RUNNING
/*!
 * @ingroup PkgCT_3Classes
 * @brief Internal per-vertex data for \cgal 3D conforming constrained Delaunay triangulations
 *
 * This class is an internal detail of the implementation of \cgal 3D conforming constrained Delaunay triangulations.
 *
 * Any model of the `ConformingConstrainedDelaunayTriangulationVertexBase_3` concept must include one object of this type
 * as a non-static data member.
 */
struct Conforming_constrained_Delaunay_triangulation_vertex_data_3 {};
#else // DOXYGEN_RUNNING

enum class CDT_3_vertex_type { FREE, CORNER, STEINER_ON_EDGE, STEINER_IN_FACE };

enum class CDT_3_vertex_marker {
  CLEAR = 0,
  REGION_BORDER,
  REGION_INSIDE,
  CAVITY,
  CAVITY_ABOVE,
  CAVITY_BELOW,
  nb_of_markers
};

class Conforming_constrained_Delaunay_triangulation_vertex_data_3 {
protected:
  // TODO: check and improve the compactness of this class
  CDT_3_vertex_type m_vertex_type = CDT_3_vertex_type::FREE;
  std::bitset<static_cast<int>(CDT_3_vertex_marker::nb_of_markers)> mark{};
  struct C_id {
    void* ptr = nullptr;
    std::size_t id = 0;
    friend bool operator==(const C_id& lhs, const C_id& rhs) {
      return lhs.ptr == rhs.ptr && lhs.id == rhs.id;
    }
  };
  union U {
    struct On_edge {
      int nb_of_incident_constraints = 0;
      C_id c_id{};
    } on_edge;
    struct On_face{
      CDT_3_signed_index face_index = 0;
    } on_face;
  } u {U::On_edge{}};

public:
  friend bool operator==(const Conforming_constrained_Delaunay_triangulation_vertex_data_3& lhs,
                         const Conforming_constrained_Delaunay_triangulation_vertex_data_3& rhs) {
    return lhs.m_vertex_type == rhs.m_vertex_type && lhs.mark == rhs.mark && std::invoke([&]() {
             if(lhs.m_vertex_type == CDT_3_vertex_type::STEINER_ON_EDGE) {
               return lhs.u.on_edge.nb_of_incident_constraints == rhs.u.on_edge.nb_of_incident_constraints &&
                      lhs.u.on_edge.c_id == rhs.u.on_edge.c_id;
             } else if(lhs.m_vertex_type == CDT_3_vertex_type::STEINER_IN_FACE) {
               return lhs.u.on_face.face_index == rhs.u.on_face.face_index;
             }
             return true;
           });
  }

  template <typename T>
  void set_on_constraint(T constraint_id) {
    ++u.on_edge.nb_of_incident_constraints;
    u.on_edge.c_id.ptr = constraint_id.vl_with_info_pointer();
    u.on_edge.c_id.id = static_cast<std::size_t>(constraint_id.index());
  }

  int number_of_incident_constraints() const {
    CGAL_assertion(u.on_edge.nb_of_incident_constraints >= 0);
    return u.on_edge.nb_of_incident_constraints;
  }

  void set_mark(CDT_3_vertex_marker marker) {
    mark.set(static_cast<unsigned int>(marker));
  }

  void clear_marks() {
    mark.reset();
  }

  void clear_mark(CDT_3_vertex_marker marker) {
    mark.reset(static_cast<unsigned int>(marker));
  }

  bool is_marked(CDT_3_vertex_marker marker) const {
    return mark.test(static_cast<unsigned int>(marker));
  }

  bool is_marked() const {
    return mark.any();
  }

  template<typename Triangulation>
  auto constraint_id(const Triangulation&) const {
    CGAL_assertion(m_vertex_type != CDT_3_vertex_type::STEINER_IN_FACE);
    using C_id = typename Triangulation::Constraint_id;
    using size_type = typename Triangulation::size_type;
    using Vertex_list_w_info_ptr = decltype(std::declval<C_id>().vl_with_info_pointer());
    auto ptr = static_cast<Vertex_list_w_info_ptr>(u.on_edge.c_id.ptr);
    auto id = static_cast<size_type>(u.on_edge.c_id.id);
    return C_id{ptr, id};
  }

  void set_Steiner_vertex_in_face(CDT_3_signed_index face_index) {
    m_vertex_type = CDT_3_vertex_type::STEINER_IN_FACE;
    u.on_face = typename U::On_face{face_index};
  }

  CDT_3_signed_index face_index() const {
    CGAL_assertion(m_vertex_type == CDT_3_vertex_type::STEINER_IN_FACE);
    return u.on_face.face_index;
  }

  CDT_3_vertex_type vertex_type() const { return m_vertex_type; }
  void set_vertex_type(CDT_3_vertex_type type) { m_vertex_type = type; }
  bool is_Steiner_vertex_on_edge() const { return m_vertex_type == CDT_3_vertex_type::STEINER_ON_EDGE; }
  bool is_Steiner_vertex_in_face() const { return m_vertex_type == CDT_3_vertex_type::STEINER_IN_FACE; }
};

#endif // DOXYGEN_RUNNING

} // namespace CGAL

#if CGAL_CXX20
#  include <concepts>

  static_assert(std::regular<CGAL::Conforming_constrained_Delaunay_triangulation_vertex_data_3>);

#endif

#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_VERTEX_DATA_3_H