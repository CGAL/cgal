// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_SIDE_TAGS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_SIDE_TAGS_H

/*! \file
 * Boundary-category tags for the 1D Arrangement on Curves package.
 *
 * Open_side_tag---the parameter space is open on that side; the endpoint of
 * the curve on that side, if it exists, is not part of the curve (e.g., either
 * side of a line, or a side of an open segment).  The edge of the arrangement
 * on that side has no vertex on that side.
 *
 * Closed_side_tag---the parameter space is closed on that side; the endpoint
 * of the curve on that side is part of the curve (e.g., a side of a closed
 * segment).  Once inserted, the endpoint is the boundary vertex on that side.
 *
 * Identified_side_tag---the left and right sides are identified (glued
 * together); the curve is cyclic, such as a circle or the direction space
 * \f$\mathbb{S}^1\f$.  The identification point is part of the curve only
 * once it is inserted.
 *
 * The categories are declared in the geometry traits.  This is because the categories govern both:
 * (a) which combinatorial structure the topology traits builds, and
 * (b) which geometric predicates the geometry traits must provide.
 *
 * A geometry traits class must export two nested type aliases:
 *
 *   using Left_side_category = <one of the three tags above>;
 *   using Right_side_category = <one of the three tags above>;
 *
 * The topology traits class reads these from the geometry traits.
 * `Arrangement_on_curve_1` and all free functions (e.g., insert, locate, overlay)
 * dispatch on them at compile time using if constexpr.
 *
 * Additional geometry-traits predicates required per category
 * -----------------------------------------------------------
 * Identified topology (Left_side_category == Identified_side_tag AND
 *                      Right_side_category == Identified_side_tag):
 *
 *   - `Is_on_x_identification_1`
 *       Functor; operator()(const Point_1& p) -> bool.
 *       Returns true iff `p` is the identification point, i.e. the point where
 *       the two identified sides meet.  Calling `Compare_x_1` with such a
 *       point as an argument is forbidden; use `Is_on_x_identification_1` to
 *       guard first.
 *       Obtained via: traits.is_on_x_identification_1_object().
 *
 * Open topology (side is Open_side_tag):
 *   No additional predicate.
 *
 * Closed topology (side is Closed_side_tag):
 *   No additional predicate for now; the boundary point is a well-defined
 *   finite `Point_1` stored as the boundary vertex of the arrangement.
 *   (A Compare_x_on_boundary_1 functor may be added in the future when
 *   rays that end exactly at the boundary need to be compared.)
 */

#include <type_traits>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// ---------------------------------------------------------------------------
// Tag types
// ---------------------------------------------------------------------------

//! The parameter space is open on this side: the curve excludes its endpoint on this side, if any.
struct Open_side_tag {};

//! The parameter space is closed on this side: the curve includes its endpoint on this side.
struct Closed_side_tag {};

//! The left and right sides are identified: the curve is cyclic.
//! Meaningful only when *both* sides carry this tag simultaneously.
struct Identified_side_tag {};

// ---------------------------------------------------------------------------
// Compile-time trait helpers
// ---------------------------------------------------------------------------

template <typename Tag>
struct Side_is_open : std::false_type {};
template <>
struct Side_is_open<Open_side_tag> : std::true_type {};

template <typename Tag>
struct Side_is_closed : std::false_type {};
template <>
struct Side_is_closed<Closed_side_tag> : std::true_type {};

template <typename Tag>
struct Side_is_identified : std::false_type {};
template <>
struct Side_is_identified<Identified_side_tag> : std::true_type {};

// Helper: true when both sides are identified (the curve is cyclic).
// GeometryTraits must expose Left_side_category and Right_side_category.
template <typename GeometryTraits>
inline constexpr bool Is_identified_v =
  Side_is_identified<typename GeometryTraits::Left_side_category>::value &&
  Side_is_identified<typename GeometryTraits::Right_side_category>::value;

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
