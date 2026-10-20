// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_RIGHT_BOUNDED_TOPOLOGY_TRAITS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_RIGHT_BOUNDED_TOPOLOGY_TRAITS_H

/*! \file
 * Topology traits for a 1D arrangement on a right-bounded master curve: (−∞, b].
 * The left boundary is open; the right boundary is closed.
 */

#include <CGAL/Arrangement_on_curve_1/Topology_traits.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1TopologyTraitsClasses
 *
 * \brief Topology traits for an arrangement on a right-bounded master curve (−∞, b].
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
          typename Allocator = std::allocator<char>>
class Right_bounded_topology_traits :
  public Topology_traits<GeometryTraits_1, Open_side_tag, Closed_side_tag, VertexData, EdgeData, UseVector, Allocator> {
  using Base =
    Topology_traits<GeometryTraits_1, Open_side_tag, Closed_side_tag, VertexData, EdgeData, UseVector, Allocator>;
public:
  using typename Base::Point_1;
  using typename Base::Vertex_descriptor;

  Vertex_descriptor create_vertex(const Point_1& p, const GeometryTraits_1& gt) {
    Vertex_descriptor v = Base::create_vertex(p);
    if (gt.parameter_space_in_x_1_object()(p) == ARR_RIGHT_BOUNDARY) this->set_right_boundary_vertex(v);
    return v;
  }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
