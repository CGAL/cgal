// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_DEFAULT_TOPOLOGY_TRAITS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_DEFAULT_TOPOLOGY_TRAITS_H

/*! \file
 * Selects the topology traits that matches the side categories of a geometry
 * traits class.
 */

#include <memory>
#include <type_traits>

#include <CGAL/Arrangement_on_curve_1/Side_tags.h>
#include <CGAL/Arrangement_on_curve_1/Open_topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/Closed_topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/Left_closed_topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/Right_closed_topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/Identified_topology_traits.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {
namespace internal {

template <typename T> struct Always_false : std::false_type {};

// The primary template handles the unsupported combinations of side
// categories, e.g., a parameter space identified on one side only.
template <typename GeometryTraits_1, typename VertexData, typename EdgeData, bool UseVector, typename Allocator,
          typename LeftSideCategory, typename RightSideCategory>
struct Default_topology_traits_impl {
  static_assert(Always_false<GeometryTraits_1>::value,
                "No topology traits matches the side categories of the geometry traits; each side must be open or "
                "closed, or both sides must be identified.");
};

template <typename GeometryTraits_1, typename VertexData, typename EdgeData, bool UseVector, typename Allocator>
struct Default_topology_traits_impl<GeometryTraits_1, VertexData, EdgeData, UseVector, Allocator,
                                    Open_side_tag, Open_side_tag> {
  using type = Open_topology_traits<GeometryTraits_1, VertexData, EdgeData, UseVector, Allocator>;
};

template <typename GeometryTraits_1, typename VertexData, typename EdgeData, bool UseVector, typename Allocator>
struct Default_topology_traits_impl<GeometryTraits_1, VertexData, EdgeData, UseVector, Allocator,
                                    Closed_side_tag, Closed_side_tag> {
  using type = Closed_topology_traits<GeometryTraits_1, VertexData, EdgeData, UseVector, Allocator>;
};

template <typename GeometryTraits_1, typename VertexData, typename EdgeData, bool UseVector, typename Allocator>
struct Default_topology_traits_impl<GeometryTraits_1, VertexData, EdgeData, UseVector, Allocator,
                                    Closed_side_tag, Open_side_tag> {
  using type = Left_closed_topology_traits<GeometryTraits_1, VertexData, EdgeData, UseVector, Allocator>;
};

template <typename GeometryTraits_1, typename VertexData, typename EdgeData, bool UseVector, typename Allocator>
struct Default_topology_traits_impl<GeometryTraits_1, VertexData, EdgeData, UseVector, Allocator,
                                    Open_side_tag, Closed_side_tag> {
  using type = Right_closed_topology_traits<GeometryTraits_1, VertexData, EdgeData, UseVector, Allocator>;
};

template <typename GeometryTraits_1, typename VertexData, typename EdgeData, bool UseVector, typename Allocator>
struct Default_topology_traits_impl<GeometryTraits_1, VertexData, EdgeData, UseVector, Allocator,
                                    Identified_side_tag, Identified_side_tag> {
  using type = Identified_topology_traits<GeometryTraits_1, VertexData, EdgeData, UseVector, Allocator>;
};

} // namespace internal

/*! The nested type `type` is the topology traits that matches the side
 * categories of `GeometryTraits_1`:
 *
 * - `Open_side_tag` on both sides: `Open_topology_traits`
 * - `Closed_side_tag` on both sides: `Closed_topology_traits`
 * - `Closed_side_tag` on the left and `Open_side_tag` on the right: `Left_closed_topology_traits`
 * - `Open_side_tag` on the left and `Closed_side_tag` on the right: `Right_closed_topology_traits`
 * - `Identified_side_tag` on both sides: `Identified_topology_traits`
 *
 * Any other combination is rejected at compile time.
 */
template <typename GeometryTraits_1,
          typename VertexData = void,
          typename EdgeData = void,
          bool UseVector = false,
          typename Allocator = std::allocator<char>>
struct Default_topology_traits {
  using type =
    typename internal::Default_topology_traits_impl<GeometryTraits_1, VertexData, EdgeData, UseVector, Allocator,
                                                    typename GeometryTraits_1::Left_side_category,
                                                    typename GeometryTraits_1::Right_side_category>::type;
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
