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
 * \brief The class template `Closed_topology_traits` is a model of the
 * `AocTopologyTraits` for 1D arrangements embedded on a closed master curve,
 * that is, a curve that includes its two endpoints, parameterized over a closed
 * interval parameter space.
 *
 * It maintains the arrangement as an alternating sequence of vertices and
 * edges. The class template supports both list-based and vector-based
 * underlying storage structures, and permits optional user data extension on
 * vertices and edges via property maps.
 *
 * A *closed* side means that the curve includes its endpoint on that side.
 * Once inserted, the endpoint is the *boundary vertex* on that side.
 *
 * \tparam GeometryTraits_1 a model of `AocTraits_1` with
 *         `Left_side_category = Closed_side_tag` and `Right_side_category = Closed_side_tag`.
 * \tparam VertexData optional user-defined data attached to each vertex (defaults to `void`).
 * \tparam EdgeData optional user-defined data attached to each edge (defaults to `void`).
 * \tparam UseVector when set to `true`, stores vertices and edges in `std::vector`
 *   to enable \f$O(\log n)\f$ binary search point location. Descriptors become `std::size_t`
 *   indices. When `false` (default), uses `std::list` where descriptors are stable
 *   bidirectional iterators.
 * \tparam Allocator allocator used for internal node and container allocations (defaults to `std::allocator<char>`).
 *
 * \cgalModels{AocTopologyTraits}
 *
 * \sa `CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>`
 */
template <typename GeometryTraits_1,
          typename VertexData = void,
          typename EdgeData = void,
          bool UseVector = false,
          typename Allocator = std::allocator<char>>
class Closed_topology_traits {
public:
  //! \name Creation
  //! @{

  //! constructs an empty topology traits with an explicit allocator object.
  explicit Closed_topology_traits(const Allocator& alloc);

  //! @}

  //! \name Boundary Vertices
  //! @{

  /*! determines whether the left boundary vertex exists.  It does not exist
   * until a point on the closed left boundary is inserted into the arrangement,
   * and does not exist again after `clear()` or after the vertex is erased.
   */
  bool has_left_boundary_vertex() const;

  /*! obtains the left boundary vertex.
   * \pre `has_left_boundary_vertex()` is `true`.
   */
  Vertex_descriptor left_boundary_vertex();

  /*! obtains the left boundary vertex (const version).
   * \pre `has_left_boundary_vertex()` is `true`.
   */
  Vertex_const_descriptor left_boundary_vertex() const;

  /*! sets the left boundary vertex to `v`.  It is called when a vertex is created
   * whose point lies on the closed left boundary, either by `create_vertex(p, gt)`
   * or by the arrangement.
   */
  void set_left_boundary_vertex(Vertex_descriptor v);

  /*! determines whether the right boundary vertex exists.  It does not exist
   * until a point on the closed right boundary is inserted into the arrangement,
   * and does not exist again after `clear()` or after the vertex is erased.
   */
  bool has_right_boundary_vertex() const;

  /*! obtains the right boundary vertex.
   * \pre `has_right_boundary_vertex()` is `true`.
   */
  Vertex_descriptor right_boundary_vertex();

  /*! obtains the right boundary vertex (const version).
   * \pre `has_right_boundary_vertex()` is `true`.
   */
  Vertex_const_descriptor right_boundary_vertex() const;

  /*! sets the right boundary vertex to `v`.  It is called when a vertex is created
   * whose point lies on the closed right boundary, either by `create_vertex(p, gt)`
   * or by the arrangement.
   */
  void set_right_boundary_vertex(Vertex_descriptor v);

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL
