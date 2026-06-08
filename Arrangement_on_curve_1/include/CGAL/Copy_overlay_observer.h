// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_COPY_OVERLAY_OBSERVER_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_COPY_OVERLAY_OBSERVER_H

#include <type_traits>
#include <boost/property_map/property_map.hpp>

namespace CGAL {
namespace Arrangement_on_curve_1 {

template <typename ArrangementA, typename ArrangementB, typename ArrangementR>
class Copy_overlay_observer {
public:
  using Vertex_const_descriptor_a = typename ArrangementA::Vertex_const_descriptor;
  using Vertex_const_descriptor_b = typename ArrangementB::Vertex_const_descriptor;
  using Vertex_descriptor_r = typename ArrangementR::Vertex_descriptor;

  using Edge_const_descriptor_a = typename ArrangementA::Edge_const_descriptor;
  using Edge_const_descriptor_b = typename ArrangementB::Edge_const_descriptor;
  using Edge_descriptor_r = typename ArrangementR::Edge_descriptor;

  using Vertex_data_map_a = typename ArrangementA::Topology_traits::Vertex_data_map;
  using Vertex_data_map_b = typename ArrangementB::Topology_traits::Vertex_data_map;
  using Vertex_data_map_r = typename ArrangementR::Topology_traits::Vertex_data_map;

  using Edge_data_map_a = typename ArrangementA::Topology_traits::Edge_data_map;
  using Edge_data_map_b = typename ArrangementB::Topology_traits::Edge_data_map;
  using Edge_data_map_r = typename ArrangementR::Topology_traits::Edge_data_map;

private:
  Vertex_data_map_r m_v_map_r;
  Edge_data_map_r m_e_map_r;

  Vertex_data_map_a m_v_map_a;
  Edge_data_map_a m_e_map_a;

  Vertex_data_map_b m_v_map_b;
  Edge_data_map_b m_e_map_b;

public:
  Copy_overlay_observer(const ArrangementA& arr_a, const ArrangementB& arr_b, ArrangementR& arr_r) :
    m_v_map_r(arr_r.vertex_data_map()),
    m_e_map_r(arr_r.edge_data_map()),
    m_v_map_a(arr_a.vertex_data_map()),
    m_e_map_a(arr_a.edge_data_map()),
    m_v_map_b(arr_b.vertex_data_map()),
    m_e_map_b(arr_b.edge_data_map())
  {}

  // 1. Two vertices coincide: defaults to copying data from vertex A
  void create_vertex(Vertex_const_descriptor_a v_a, Vertex_const_descriptor_b, Vertex_descriptor_r v_res) {
    using Vertex_data_r = typename Vertex_data_map_r::value_type;
    if constexpr (! std::is_void_v<Vertex_data_r>) {
      using Vertex_data_a = typename Vertex_data_map_a::value_type;
      if constexpr (std::is_convertible_v<Vertex_data_a, Vertex_data_r>) {
        put(m_v_map_r, v_res, get(m_v_map_a, v_a));
      }
    }
  }

  // 2. Vertex A splits Edge B: copies data from vertex A
  void create_vertex(Vertex_const_descriptor_a v_a, Edge_const_descriptor_b, Vertex_descriptor_r v_res) {
    using Vertex_data_r = typename Vertex_data_map_r::value_type;
    if constexpr (! std::is_void_v<Vertex_data_r>) {
      using Vertex_data_a = typename Vertex_data_map_a::value_type;
      if constexpr (std::is_convertible_v<Vertex_data_a, Vertex_data_r>) {
        put(m_v_map_r, v_res, get(m_v_map_a, v_a));
      }
    }
  }

  // 3. Edge A is split by Vertex B: copies data from vertex B
  void create_vertex(Edge_const_descriptor_a, Vertex_const_descriptor_b v_b, Vertex_descriptor_r v_res) {
    using Vertex_data_r = typename Vertex_data_map_r::value_type;
    if constexpr (! std::is_void_v<Vertex_data_r>) {
      using Vertex_data_b = typename Vertex_data_map_b::value_type;
      if constexpr (std::is_convertible_v<Vertex_data_b, Vertex_data_r>) {
        put(m_v_map_r, v_res, get(m_v_map_b, v_b));
      }
    }
  }

  // 4. Two edges overlap: defaults to copying data from edge A
  void create_edge(Edge_const_descriptor_a e_a, Edge_const_descriptor_b, Edge_descriptor_r e_res) {
    using Edge_data_r = typename Edge_data_map_r::value_type;
    if constexpr (! std::is_void_v<Edge_data_r>) {
      using Edge_data_a = typename Edge_data_map_a::value_type;
      if constexpr (std::is_convertible_v<Edge_data_a, Edge_data_r>) {
        put(m_e_map_r, e_res, get(m_e_map_a, e_a));
      }
    }
  }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
