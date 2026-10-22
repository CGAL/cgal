// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_OPEN_TOPOLOGY_TRAITS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_OPEN_TOPOLOGY_TRAITS_H

/*! \file
 * Topology traits for a 1D arrangement on an open master curve.
 * Both boundaries are open; no boundary vertices are ever registered.
 */

#include <CGAL/Arrangement_on_curve_1/Topology_traits_base.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1TopologyTraitsClasses
 *
 * \brief Topology traits for an arrangement on an open master curve.
 *
 * Both boundaries of the parameter space are open: the master curve excludes
 * its endpoints, if any; e.g., a line or an open segment.  No boundary vertices are ever registered.
 *
 * Derives from `Topology_traits_base`; `create_vertex` ignores the geometry
 * traits, since no point lies on a boundary.
 *
 * \tparam GeometryTraits_1 a model of `AocTraits_1` with
 *         `Left_side_category = Open_side_tag` and
 *         `Right_side_category = Open_side_tag`.
 * \tparam VertexData user data attached to each vertex (default: `void`).
 * \tparam EdgeData   user data attached to each edge   (default: `void`).
 * \tparam UseVector  if `true`, containers use `std::vector` with index
 *                    descriptors; otherwise `std::list` with iterator
 *                    descriptors (default: `false`).
 * \tparam Allocator  standard allocator (default: `std::allocator<char>`).
 */
template <typename GeometryTraits_1,
          typename VertexData = void,
          typename EdgeData = void,
          bool UseVector = false,
          typename Allocator_ = std::allocator<char>>
class Open_topology_traits :
  public Topology_traits_base<typename GeometryTraits_1::Point_1, VertexData, EdgeData, UseVector, Allocator_> {
public:
  using Geometry_traits_1 = GeometryTraits_1;
  using Vertex_data = VertexData;
  using Edge_data = EdgeData;
  using Allocator = Allocator_;

  //! The standard stl type.
  using allocator_type = Allocator;

private:
  using Base = Topology_traits_base<typename GeometryTraits_1::Point_1, VertexData, EdgeData, UseVector, Allocator>;

public:
  using Base::Base;
  using typename Base::Point_1;
  using typename Base::Vertex_descriptor;

  //! creates a new vertex associated with point `p`.
  Vertex_descriptor create_vertex(const Point_1& p, const GeometryTraits_1& /* gt */) { return Base::create_vertex(p); }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
