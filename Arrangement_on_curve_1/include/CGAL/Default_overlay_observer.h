// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_DEFAULT_OVERLAY_OBSERVER_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_DEFAULT_OVERLAY_OBSERVER_H

namespace CGAL {
namespace Arrangement_on_curve_1 {

template <typename ArrangementA, typename ArrangementB, typename ArrangementR>
class Default_overlay_observer {
public:
  using Vertex_descriptor_a = typename ArrangementA::Vertex_descriptor;
  using Vertex_descriptor_b = typename ArrangementB::Vertex_descriptor;
  using Vertex_descriptor_r = typename ArrangementR::Vertex_descriptor;

  using Edge_descriptor_a   = typename ArrangementA::Edge_descriptor;
  using Edge_descriptor_b   = typename ArrangementB::Edge_descriptor;
  using Edge_descriptor_r   = typename ArrangementR::Edge_descriptor;

  using Vertex_data_map_a = typename ArrangementA::Topology_traits::Vertex_data_map;
  using Vertex_data_map_b = typename ArrangementB::Topology_traits::Vertex_data_map;
  using Vertex_data_map_r = typename ArrangementR::Topology_traits::Vertex_data_map;

  using Edge_data_map_a   = typename ArrangementA::Topology_traits::Edge_data_map;
  using Edge_data_map_b   = typename ArrangementB::Topology_traits::Edge_data_map;
  using Edge_data_map_r   = typename ArrangementR::Topology_traits::Edge_data_map;

public:
  Default_overlay_observer(const ArrangementA&, const ArrangementB&, ArrangementR&) {}

  // All notification methods are intentional no-ops
  void create_vertex(Vertex_descriptor_a, Vertex_descriptor_b, Vertex_descriptor_r) {}
  void create_vertex_from_a(Vertex_descriptor_a, Vertex_descriptor_r) {}
  void create_vertex_from_b(Vertex_descriptor_b, Vertex_descriptor_r) {}

  void create_edge(Edge_descriptor_a, Edge_descriptor_b, Edge_descriptor_r) {}
  void create_edge_from_a(Edge_descriptor_a, Edge_descriptor_r) {}
  void create_edge_from_b(Edge_descriptor_b, Edge_descriptor_r) {}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
