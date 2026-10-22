// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_CLOSED_TOPOLOGY_TRAITS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_CLOSED_TOPOLOGY_TRAITS_H

#include <CGAL/Arrangement_on_curve_1/Topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/Parameter_space.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1TopologyTraitsClasses
 *
 * \brief Topology traits for an arrangement on a closed master curve.
 *
 * Derives from `Topology_traits` with the left side open and the right side
 * closed, and overrides `create_vertex` to register the vertex as the right
 * boundary vertex when the point lies on the closed right boundary, using
 * `GeometryTraits_1::Parameter_space_in_x_1`.
 */
template <typename GeometryTraits_1,
          typename VertexData = void,
          typename EdgeData = void,
          bool UseVector = false,
          typename Allocator_ = std::allocator<char>>
class Closed_topology_traits :
  public Topology_traits<GeometryTraits_1, Closed_side_tag, Closed_side_tag, VertexData, EdgeData, UseVector,
                         Allocator_> {
public:
  using Geometry_traits_1 = GeometryTraits_1;
  using Vertex_data = VertexData;
  using Edge_data = EdgeData;
  using Allocator = Allocator_;

  // The standard stl type
  using allocator_type = Allocator_;

private:
  using Base =
    Topology_traits<Geometry_traits_1, Closed_side_tag, Closed_side_tag, Vertex_data, Edge_data, UseVector, Allocator>;

public:
  using Base::Base;
  using typename Base::Point_1;
  using typename Base::Vertex_descriptor;

  //! creates a new vertex ssociated with point `p`.
  Vertex_descriptor create_vertex(const Point_1& p, const GeometryTraits_1& gt) {
    Vertex_descriptor v = Base::create_vertex(p);
    auto ps = gt.parameter_space_in_x_1_object()(p);
    if (ps == Parameter_space::LEFT_BOUNDARY) this->set_left_boundary_vertex(v);
    else if (ps == Parameter_space::RIGHT_BOUNDARY) this->set_right_boundary_vertex(v);
    return v;
  }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
