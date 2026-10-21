// Copyright (c) 2026 Tel-Aviv University (Israel). All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1TopologyTraitsClasses
 *
 * \brief `Identified_topology_traits` is a topology traits model of
 * `AocTopologyTraits` for 1D arrangements on a closed loop (identified
 * boundary).
 *
 * When both sides of the parameter space are identified (tag:
 * `Identified_side_tag`), the curve is a closed loop, for example a circle,
 * a direction-space \f$\mathbb{S}_1\f$, or any periodic 1D domain.
 *
 * Key topological differences from `Open_topology_traits`:
 *
 * - There is **no** sentinel edge.  An empty arrangement does not contain vertices.
 *
 * - Vertices and edges form a cyclic sequence:
 *   \f$\ldots \leftarrow e_0 \rightarrow v_0 \leftarrow e_1 \rightarrow v_1 \ldots\f$
 *   where the last edge wraps back to the first vertex.
 *
 * - A edge incident on the left to a vertex that maps to the identification point
 *   provides the canonical entry point for cyclic traversal (locate, overlay).
 *
 * - `insert_empty(p)` creates the very first vertex and a single loop edge
 *   (left and right both pointing to that vertex).
 *
 * All shared combinatorial infrastructure is inherited from
 * `Topology_traits_base`.  This class adds only `m_start_edge` / `m_has_start`
 * and the `start_edge()` / `set_start_edge()` API.
 *
 * \tparam Point_1_ the geometric point type stored at each vertex.
 * \tparam VertexData optional user data per vertex (`void` = none, default).
 * \tparam EdgeData optional user data per edge (`void` = none, default).
 * \tparam UseVector if `true`, use `std::vector` with `std::size_t` index
 *   descriptors; if `false` (default), use `std::list` with stable iterator
 *   descriptors.
 * \tparam Allocator the allocator used for internal containers.
 *
 * \cgalModels{AocTopologyTraits}
 *
 * \sa `Topology_traits_base`
 * \sa `Open_topology_traits`
 * \sa `Closed_topology_traits`
 */
template <typename GeometryTraits_1,
          typename VertexData = void,
          typename EdgeData = void,
          bool UseVector = false,
          typename Allocator = std::allocator<char>>
class Identified_topology_traits {
public:
  /// \name Creation
  /// @{

  //! constructs an empty topology traits with an explicit allocator object.
  explicit Identified_topology_traits(const Allocator& alloc);

  /// @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL
