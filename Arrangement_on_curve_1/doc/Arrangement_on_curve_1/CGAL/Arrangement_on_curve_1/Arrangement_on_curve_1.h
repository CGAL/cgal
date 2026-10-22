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
 * space can be open, half-open, closed, or identified. `GeometryTraits_1`
 * declares the types `Left_side_category` and `Right_side_category`, each of
 * which is one of `Open_side_tag`, `Closed_side_tag`, or `Identified_side_tag`.
 * The static data member, for example, `is_identified` is `true` if and only if
 * both categories are `Identified_side_tag`; that is, the curve is a closed
 * loop such as \f$\mathbb{S}^1\f$. When `TopologyTraits` is omitted, it
 * defaults to `Open_topology_traits<GeometryTraits_1>`, which is correct for
 * any geometry traits with `Left_side_category = Right_side_category =
 * Open_side_tag`.  For other combinations supply the matching traits
 * explicitly, e.g.:
 *
 * \code
 * using Arr = Arrangement_on_curve_1<Circle_traits_1<Kernel>>,
 *   Identified_topology_traits<Circle_traits_1<Kernel>>;
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
 *   `Open_topology_traits<GeometryTraits_1::Point_1>`.
 * \tparam BinarySearch enable \f$O(log n)\f$ binary-search locate
 *   (requires `TopologyTraits::UseVector = true`; default: `false`).
 *
 * \cgalModels{ArrangementOnCurve_1}
 */
template <typename GeometryTraits_1,
          typename TopologyTraits = Open_topology_traits<GeometryTraits_1>,
          bool BinarySearch = false>
class Arrangement_on_curve_1 {
public:
  //! \name Types
  //! @{

  //! The geometry traits type.
  typedef GeometryTraits_1 Geometry_traits_1;

  //! The topology traits type.
  typedef TopologyTraits Topology_traits;

  //! A shared-ownership pointer to a const geometry traits object.
  typedef std::shared_ptr<const Geometry_traits_1> Shared_geometry_traits;

  //! Whether \f$O(log n)\f$ binary-search locate is enabled.
  static constexpr bool binary_search_enabled = BinarySearch;

  typedef typename Topology_traits::Point_1                Point_1;
  typedef typename Topology_traits::Vertex_descriptor      Vertex_descriptor;
  typedef typename Topology_traits::Edge_descriptor        Edge_descriptor;
  typedef typename Topology_traits::Vertex_const_descriptor Vertex_const_descriptor;
  typedef typename Topology_traits::Edge_const_descriptor  Edge_const_descriptor;
  typedef typename Topology_traits::Size                   Size;

  /*! A point-location result: a vertex or edge descriptor.
   */
  typedef std::variant<Vertex_descriptor, Edge_descriptor> Location_result;

  /*! A const point-location result.
   */
  typedef std::variant<Vertex_const_descriptor, Edge_const_descriptor> Const_location_result;

  //! Tag type describing the topological behaviour of the left boundary.
  //! Equals `GeometryTraits_1::Left_side_category`.
  typedef typename Geometry_traits_1::Left_side_category  Left_side_category;

  //! Tag type describing the topological behaviour of the right boundary.
  //! Equals `GeometryTraits_1::Right_side_category`.
  typedef typename Geometry_traits_1::Right_side_category Right_side_category;

  //! @}

  //! \name Static Predicates
  //! @{

  /*! is `true` if both `Left_side_category` and `Right_side_category` are
   * `Identified_side_tag`.
   */
  static constexpr bool is_identified;

  //! @}

  //! \name Creation
  //! @{

  /*! default constructor; creates an empty arrangement.
   */
  Arrangement_on_curve_1();

  /*! constructs with a given geometry traits object.
   * \param gt the geometry traits to use for comparisons.
   */
  explicit Arrangement_on_curve_1(const Geometry_traits_1& gt);

  /*! constructs with a shared-ownership pointer to the geometry traits.
   * \param sgt a shared pointer to the geometry traits.
   */
  explicit Arrangement_on_curve_1(Shared_geometry_traits sgt);

  /*! constructs with a shared-ownership pointer to the geometry traits and an
   * explicit topology traits object.
   * \param sgt a shared pointer to the geometry traits.
   * \param tt the topology traits to use for storage and navigation.
   */
  Arrangement_on_curve_1(Shared_geometry_traits sgt, Topology_traits tt);

  /*! replaces the shared geometry traits pointer.
   *
   * Transfers shared ownership of a new geometry traits object into the
   * arrangement.  The arrangement immediately uses the new traits for all
   * subsequent geometric operations.  The previous traits object is released
   * (and destroyed if no other owner holds a reference to it).
   *
   * \param new_sgt the new shared-ownership pointer to the geometry traits.
   */
  void reset_shared_geometry_traits(Shared_geometry_traits new_sgt);

  //! @}

  //! \name Access
  //! @{

  /*! obtains a pointer to the underlying geometry traits.
   */
  const Geometry_traits_1* geometry_traits() const;

  /*! obtains the shared-ownership pointer to the geometry traits.
   */
  Shared_geometry_traits shared_geometry_traits_1() const;

  /*! obtains a const reference to the internal topology traits object.
   */
  const Topology_traits& topology_traits() const;

  /*! obtains a mutable reference to the internal topology traits object.
   */
  Topology_traits& topology_traits();

  //! @}

  //! \name Modifiers
  //! @{

  /*! resets the shared smart pointer allowing empty arrangements to bind to a separate existing traits memory frame.
   * \pre `is_empty() == true`
   */
  void reset_shared_geometry_traits(Shared_geometry_traits new_shared_traits_traits);

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL
