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
 * \brief The class template `Default_topology_traits` selects the topology
 * traits that matches the side categories of the geometry traits.
 *
 * It is the default of the template parameter `TopologyTraits` of
 * `Arrangement_on_curve_1`.  The selection is:
 *
 * | `Left_side_category`  | `Right_side_category` | `type`                         |
 * |-----------------------|-----------------------|--------------------------------|
 * | `Open_side_tag`       | `Open_side_tag`       | `Open_topology_traits`         |
 * | `Closed_side_tag`     | `Closed_side_tag`     | `Closed_topology_traits`       |
 * | `Closed_side_tag`     | `Open_side_tag`       | `Left_closed_topology_traits`  |
 * | `Open_side_tag`       | `Closed_side_tag`     | `Right_closed_topology_traits` |
 * | `Identified_side_tag` | `Identified_side_tag` | `Identified_topology_traits`   |
 *
 * Any other combination is rejected at compile time.
 *
 * \tparam GeometryTraits_1 a model of `AocTraits_1`.
 * \tparam VertexData optional user-defined data attached to each vertex (defaults to `void`).
 * \tparam EdgeData optional user-defined data attached to each edge (defaults to `void`).
 * \tparam UseVector when set to `true`, stores vertices and edges in `std::vector`; when `false` (default), uses
 *   `std::list`.
 * \tparam Allocator allocator used for internal node and container allocations (defaults to `std::allocator<char>`).
 *
 * \sa `CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>`
 */
template <typename GeometryTraits_1,
          typename VertexData = void,
          typename EdgeData = void,
          bool UseVector = false,
          typename Allocator = std::allocator<char>>
struct Default_topology_traits {
  //! \name Types
  //! @{

  /*! the selected topology traits, instantiated with `GeometryTraits_1`,
   * `VertexData`, `EdgeData`, `UseVector`, and `Allocator`; a model of
   * `AocTopologyTraits`.
   */
  typedef unspecified_type type;

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL
