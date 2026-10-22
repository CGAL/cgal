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
 * `AocTopologyTraits` for 1D arrangements on a cyclic curve (identified
 * boundary).
 *
 * When both sides of the parameter space are identified (tag:
 * `Identified_side_tag`), the curve is cyclic, for example a circle,
 * a direction-space \f$\mathbb{S}^1\f$, or any periodic 1D domain.
 *
 * Key topological differences from `Open_topology_traits`:
 *
 * - An empty arrangement consists of a single edge that spans the entire
 *   curve and no vertices, as in all topology traits; however, every edge of a
 *   non-empty arrangement has a vertex on both sides.
 *
 * - Vertices and edges form a cyclic sequence:
 *   \f$\ldots \leftarrow e_0 \rightarrow v_0 \leftarrow e_1 \rightarrow v_1 \ldots\f$
 *   where the last edge wraps back to the first vertex.
 *
 * - A edge incident on the left to a vertex that maps to the identification point
 *   provides the canonical entry point for cyclic traversal (locate, overlay).
 *
 * - \link ArrangementOnCurve_1::insert_empty() `insert_empty(p)`\endlink creates the very first
 *   vertex and a single loop edge (left and right both pointing to that vertex).
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
 * \sa `Open_topology_traits`
 * \sa `Closed_topology_traits`
 * \sa `Left_closed_topology_traits`
 * \sa `Right_closed_topology_traits`
 */
template <typename GeometryTraits_1,
          typename VertexData = void,
          typename EdgeData = void,
          bool UseVector = false,
          typename Allocator = std::allocator<char>>
class Identified_topology_traits {
public:
  //! \name Creation
  //! @{

  //! constructs an empty topology traits with an explicit allocator object.
  explicit Identified_topology_traits(const Allocator& alloc);

  //! @}

  //! \name Identification Vertex
  //! @{

  /*! determines whether the identification vertex exists.  It does not exist
   * until the identification point is inserted into the arrangement, and does
   * not exist again after `clear()` or after the vertex is erased.
   */
  bool has_identification_vertex() const;

  /*! obtains the identification vertex.
   * \pre `has_identification_vertex()` is `true`.
   */
  Vertex_descriptor identification_vertex();

  /*! obtains the identification vertex (const version).
   * \pre `has_identification_vertex()` is `true`.
   */
  Vertex_const_descriptor identification_vertex() const;

  /*! sets the identification vertex to `v`.  It is called by `create_vertex()`
   * when the point of the new vertex is the identification point.
   */
  void set_identification_vertex(Vertex_descriptor v);

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL
