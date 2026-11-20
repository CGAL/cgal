// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1TopologyTraitsClasses
 *
 * \brief The class template `Open_topology_traits` is a model of the
 * `AocTopologyTraits` concept for 1D arrangements embedded on an open
 * master curve, that is, a curve that excludes its endpoints, if any,
 * parameterized over an open interval parameter space.
 *
 * It maintains the arrangement as an alternating sequence of vertices and
 * edges. The class template supports both list-based and vector-based
 * underlying storage structures, and permits optional user data extension on
 * vertices and edges via property maps.
 *
 * \tparam GeometryTraits_1 a model of `AocTraits_1` with
 *         `Left_side_category = Open_side_tag` and `Right_side_category = Open_side_tag`.
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
template <typename GeometryTraits_1, typename VertexData = void, typename EdgeData = void,
          bool UseVector = false, typename Allocator = std::allocator<char>>
class Open_topology_traits {
public:
  //! \name Creation
  //! @{

  //! constructs an empty topology traits with an explicit allocator object.
  explicit Open_topology_traits(const Allocator& alloc);

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL
