// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_UNBOUNDED_TOPOLOGY_TRAITS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_UNBOUNDED_TOPOLOGY_TRAITS_H

/*! \file
 * Topology traits for a 1D arrangement on an unbounded master curve: (−∞, +∞).
 * Both boundaries are open; no boundary vertices are ever registered.
 */

#include <CGAL/Arrangement_on_curve_1/Topology_traits.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1TopologyTraitsClasses
 *
 * \brief Topology traits for an arrangement on an unbounded master curve \f$(−\infty, +\infty)\f$.
 *
 * Both boundaries of the parameter space are open: neither endpoint of the
 * master curve is in the domain.  No boundary vertices are ever registered.
 *
 * Derives from `Topology_traits` with both sides set to `Open_side_tag`.
 * The base-class `create_vertex` already handles the open/open case correctly
 * (it ignores the geometry traits), so no override is needed here.
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
          typename Allocator = std::allocator<char>>
class Unbounded_topology_traits :
  public Topology_traits<GeometryTraits_1, Open_side_tag, Open_side_tag, VertexData, EdgeData, UseVector, Allocator>
{};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
