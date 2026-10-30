// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1Ref
 *
 * \brief `Arrangement_on_curve_1` represents a subdivision of a continuous 1D
 * curve, referred to as the *master curve*, into maximally connected subcurves.
 * The subdivision is represented as a linked list of alternating vertices and
 * edges, referred to as cells, where a vertex maps to a point and an edge maps
 * to an open subcurve.  Thus, the operations that obtain incident cells
 * provided by the class template take constant time, and the space needed to
 * store the arrangement is linear in the complexity of the arrangement.
 *
 * ### Geometry and Topology Traits
 *
 * The class is parameterized by a geometry traits (`GeometryTraits_1`) and a
 * topology traits (`TopologyTraits`).  The geometry traits governs the
 * geometric predicates and boundary-category tags; the topology traits governs
 * the underlying combinatorial data structure (storage, descriptors, and
 * sentinel/boundary book-keeping).
 *
 * ### Parameter Space and Boundary Categories
 *
 * The master curve \f$C = \phi_C(X)\f$ is given by a continuous function
 * \f$\phi_C: X \rightarrow \mathbb{R}^d\f$; it is defined over a
 * one-dimensional interval \f$X\f$ called the parameter space. The parameter
 * space can be open, closed, left-closed, right-closed, or identified. `GeometryTraits_1`
 * declares the types `Left_side_category` and `Right_side_category`, each of
 * which is one of `Open_side_tag`, `Closed_side_tag`, or `Identified_side_tag`.
 * The static data member, for example, `is_identified` is `true` if and only if
 * both categories are `Identified_side_tag`; that is, the curve is cyclic,
 * such as \f$\mathbb{S}^1\f$. When `TopologyTraits` is omitted, it
 * defaults to `Default_topology_traits<GeometryTraits_1>::type`, the topology
 * traits that matches the two categories, e.g.,
 * `Identified_topology_traits<Circle_traits_1<Kernel>>` for
 * `Circle_traits_1<Kernel>`.  Supply the topology traits explicitly to attach
 * user data to vertices or edges, or to use vector storage, e.g.:
 *
 * \code
 * using Gt = Circle_traits_1<Kernel>;
 * using Arr = Arrangement_on_curve_1<Gt, Identified_topology_traits<Gt, Vertex_data, void, true>>;
 * \endcode
 *
 * ### Traits Object Storing

 * An arrangement object maintains a smart pointer to a shared immutable
 * geometry traits object, optimizing lifecycle tracking and enabling efficient
 * data synchronization across topological overlays.
 *
 * ### Binary-search mode
 *
 * When `BinarySearch = true` and `TopologyTraits::UseVector = true`, the
 * `locate()` function uses binary search over the sorted vertex vector,
 * reducing point location from \f$O(n)\f$ to \f$O(log n)\f$. Binary-search mode
 * is not supported for an identified parameter space.
 *
 * \tparam GeometryTraits_1 a model of `AocTraits_1` and optionally
 *   `AocClosedBoundaryTraits_1` or `AocIdentifiedBoundaryTraits_1`.
 * \tparam TopologyTraits a model of `AocTopologyTraits`; defaults to
 *   `Default_topology_traits<GeometryTraits_1>::type`.
 * \tparam BinarySearch enable \f$O(log n)\f$ binary-search locate
 *   (requires `TopologyTraits::UseVector = true`; default: `false`).
 *
 * \cgalModels{ArrangementOnCurve_1}
 */
template <typename GeometryTraits_1,
          typename TopologyTraits = typename Default_topology_traits<GeometryTraits_1>::type,
          bool BinarySearch = false>
class Arrangement_on_curve_1 {
public:
  //! \name Types
  //! @{

  //! the geometry traits type.
  typedef GeometryTraits_1 Geometry_traits_1;

  //! the topology traits type.
  typedef TopologyTraits Topology_traits;

  //! a shared-ownership pointer to a const geometry traits object.
  typedef std::shared_ptr<const Geometry_traits_1> Shared_geometry_traits;

  //! the allocator type of the topology traits.
  typedef typename Topology_traits::Allocator Allocator;

  //! the standard STL allocator type.
  typedef Allocator allocator_type;

  //! @}

  //! \name Static Predicates
  //! @{

  //! is `true` if binary-search point location is enabled, that is, `BinarySearch` is `true`.
  static constexpr bool binary_search_enabled = BinarySearch;

  //! is `true` if `Left_side_category` is `Closed_side_tag`.
  static constexpr bool left_is_closed;

  //! is `true` if `Right_side_category` is `Closed_side_tag`.
  static constexpr bool right_is_closed;

  //! @}

  //! \name Creation
  //! @{

  /*! constructs an empty arrangement with a default-constructed geometry
   * traits object.
   */
  Arrangement_on_curve_1();

  /*! constructs an empty arrangement that shares the geometry traits object
   * pointed to by `sgt`.
   * \pre `sgt` is not null.
   */
  Arrangement_on_curve_1(Shared_geometry_traits sgt);

  /*! constructs an arrangement that shares the geometry traits object pointed
   * to by `sgt` and has the topology traits `tt`.
   * \pre `sgt` is not null.
   */
  Arrangement_on_curve_1(Shared_geometry_traits sgt, Topology_traits tt);

  /*! constructs an empty arrangement that shares the geometry traits object
   * pointed to by `sgt`, and whose topology traits uses the allocator `alloc`.
   * \pre `sgt` is not null.
   */
  Arrangement_on_curve_1(Shared_geometry_traits sgt, const Allocator& alloc);

  //! @}

  //! \name Access Functions
  //! @{

  /*! obtains the geometry traits.
   */
  const Geometry_traits_1& geometry_traits_1() const;

  /*! obtains the shared-ownership pointer to the geometry traits.
   */
  Shared_geometry_traits shared_geometry_traits_1() const;

  /*! obtains the topology traits.
   */
  const Topology_traits& topology_traits() const;

  /*! obtains the topology traits (non-const version).
   */
  Topology_traits& topology_traits();

  /*! obtains the allocator of the topology traits.
   */
  Allocator get_allocator() const;

  //! @}

  //! \name Modifiers
  //! @{

  /*! replaces the shared geometry traits pointer with `sgt`.  The previous
   * geometry traits object is released, and destroyed if no other owner holds
   * a reference to it.
   * \pre `empty()` is `true`.
   * \pre `sgt` is not null.
   */
  void reset_shared_geometry_traits(Shared_geometry_traits sgt);

  /*! inserts a vertex at the left boundary.  The vertex becomes the left
   * vertex of \link ArrangementOnCurve_1::left_boundary_edge() `left_boundary_edge()`\endlink; no edge is created.
   * \pre `left_is_closed` is `true`.
   * \pre \link ArrangementOnCurve_1::left_boundary_edge() `left_boundary_edge()`\endlink has no left vertex.
   * \pre `p` lies on the left boundary.
   */
  Vertex_descriptor insert_left_boundary_vertex(const Point_1& p);

  /*! inserts a vertex at the right boundary.  The vertex becomes the right
   * vertex of `e`; no edge is created.
   * \pre `right_is_closed` is `true`.
   * \pre `e` is \link ArrangementOnCurve_1::right_boundary_edge() `right_boundary_edge()`\endlink and has no right
   * vertex.
   * \pre `p` lies on the right boundary.
   */
  Vertex_descriptor insert_right_boundary_vertex(Edge_descriptor e, const Point_1& p);

  /*! inserts a vertex at the identification point.  The vertex becomes the
   * left vertex of
   * \link ArrangementOnCurve_1::left_boundary_edge() `left_boundary_edge()`\endlink,
   * and a new edge is created between the last vertex and the new vertex.
   * \pre `is_identified` is `true`.
   * \pre The arrangement is not empty and does not contain the identification vertex.
   * \pre `p` is the identification point.
   */
  Vertex_descriptor insert_identification_vertex(const Point_1& p);

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL
