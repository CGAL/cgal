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
 * `AocTopologyTraits` for 1D arrangements embedded on a closed master curve
 * parameterized over a closed interval parameter space with (finite) boundary
 * endpoints.
 *
 * It maintains the arrangement as an alternating sequence of vertices and
 * edges. The class template supports both list-based and vector-based
 * underlying storage structures, and permits optional user data extension on
 * vertices and edges via property maps.
 *
 * A *closed* boundary side means the carrier curve has a known, finite
 * endpoint on that side.  The arrangement stores a *boundary vertex* at each
 * closed end rather than an unbounded sentinel edge.
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

  //! constructs an empty topology traits initialized with a single unbounded edge.
  Closed_topology_traits();

  //! constructs an empty topology traits with an explicit allocator object.
  explicit Closed_topology_traits(const Allocator_type& alloc);

  //! copy constructor.
  Closed_topology_traits(const Closed_topology_traits& other);

  //! copy assignment.
  Closed_topology_traits& operator=(const Closed_topology_traits& other);

  //! move constructor.
  Closed_topology_traits(Closed_topology_traits&& other) noexcept;

  //! move assignment.
  Closed_topology_traits& operator=(Closed_topology_traits&& other) noexcept;

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL
