// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_CIRCLE_TRAITS_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_CIRCLE_TRAITS_1_H

/*! Geometry traits for a 1D arrangement of points on the unit circle
 * \f$\mathbb{S}^1\f$.
 *
 * This is the 1D analogue of `Arr_geodesic_arc_on_sphere_traits_2`, which
 * provides a traits for points (and geodesic arcs) on the unit sphere
 * \f$\mathbb{S}^2\f$.  Here, "points" are directions on the plane (elements of
 * \f$\mathbb{S}^1\f$) represented as `Direction_2` objects from the supplied
 * Kernel.  Using `Direction_2` rather than a floating-point angle gives exact
 * arithmetic when the kernel is exact (e.g.,
 * `CGAL::Exact_predicates_exact_constructions_kernel`).
 *
 * Parameter space and boundary categories
 * ---------------------------------------
 * The parameter space of \f$\mathbb{S}^1\f$ is \f$[\theta_{\mathrm{ref}}, \theta_{\mathrm{ref}} + 2\pi)\f$ with the
 * two endpoints identified, where \f$\theta_{\mathrm{ref}}\f$ is the angle of the reference
 * direction `d_ref` (see below).  Therefore:
 *
 * - `Left_side_category`  = `Identified_side_tag`
 * - `Right_side_category` = `Identified_side_tag`
 *
 * A "reference direction" `d_ref` (default: positive `x`-axis) defines the
 * position of the identification point.  The total order on directions is
 * the counter-clockwise order starting strictly after `d_ref` and wrapping
 * around back to `d_ref`:
 *
 *   \f$d_{\mathrm{ref}} < d_1 < d_2 < \ldots < d_n < d_{\mathrm{ref}}\f$
 *
 * The identification point is `d_ref`: a direction equal to `d_ref` lies
 * exactly on the identification point.  `Is_on_x_identification_1` returns
 * `true` for such a direction.  `Compare_x_1` must not be called with `d_ref`
 * as an argument; guard first with `Is_on_x_identification_1`.
 *
 * Angular total order
 * -------------------
 * Given a reference direction `d_ref` with angle \f$\theta_{\mathrm{ref}}\f$, define the
 * open angular interval \f$(\theta_{\mathrm{ref}}, \theta_{\mathrm{ref}} + 2\pi)\f$.  A
 * direction \f$d_1 < d_2\f$ in this order iff \f$d_1\f$ appears before \f$d_2\f$ when sweeping counter-clockwise
 * starting strictly after `d_ref`.
 *
 * `Direction_2` representation
 * ----------------------------
 * A CGAL `Direction_2` wraps a non-zero vector and considers two vectors equal
 * iff they are positive multiples of each other.  Opposite directions are
 * distinct.  The comparison below operates on the directions directly through
 * the kernel predicate `Counterclockwise_in_between_2`.
 *
 * Required Kernel concepts
 * ------------------------
 * The Kernel must provide:
 * - `Kernel::Direction_2`
 * - `Kernel::Counterclockwise_in_between_2`
 */

#include <CGAL/enum.h>
#include <CGAL/Arrangement_on_curve_1/Side_tags.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// ============================================================================
// Circle_traits_1
// ============================================================================
template <typename Kernel_>
class Circle_traits_1 : public Kernel_ {
public:
  using Kernel = Kernel_;

  // --------------------------------------------------------------------------
  // Type aliases imported from the kernel
  // --------------------------------------------------------------------------
  using FT = typename Kernel::FT;
  using Direction_2 = typename Kernel::Direction_2;

  // In this traits, a "point" on \f$\mathbb{S}^1\f$ is a `Direction_2`.
  using Point_1 = Direction_2;

  // --------------------------------------------------------------------------
  // Boundary-category tags (in the geometry traits, following the convention
  // of `Arrangement_on_surface_2` where the traits owns the side categories).
  //
  // Both sides of the circle's parameter space are identified: the curve is
  // the cyclic curve \f$\mathbb{S}^1\f$.
  // --------------------------------------------------------------------------
  using Left_side_category = Identified_side_tag;
  using Right_side_category = Identified_side_tag;

  /*! constructs from a reference direction defining the identification point.
   *
   * An optional reference direction sets the position of the identification
   * point.  The reference direction is the unique point where the two
   * identified sides meet.  The default identification point is the positive
   * `x`-axis direction \f$(1, 0)\f$.
   *
   * \param reference the identification-point.
   */
  explicit Circle_traits_1(const Direction_2& reference = Direction_2(1, 0)) : m_ref(reference) {}

  const Direction_2& reference_direction() const { return m_ref; }

  /*! Is_on_x_identification_1
   *
   * Returns `true` iff the given direction is the identification point, i.e.,
   * it equals the reference direction.
   *
   * Calling `Compare_x_1` with such a direction as an argument is forbidden.
   * Guard all `Compare_x_1` calls with Is_on_x_identification_1 first.
   */
  class Is_on_x_identification_1 {
  public:
    explicit Is_on_x_identification_1(const Circle_traits_1& tr) : m_tr(tr) {}

    bool operator()(const Point_1& d) const { return d == m_tr.m_ref; }

  private:
    const Circle_traits_1& m_tr;
  };

  Is_on_x_identification_1 is_on_x_identification_1_object() const { return Is_on_x_identification_1(*this); }

  /*! Compare_x_1
   *
   * Strict total order on the directions whose angles lie in
   * \f$(\theta_{\mathrm{ref}}, \theta_{\mathrm{ref}} + 2\pi)\f$:
   *   SMALLER if `d1` comes before `d2` in the CCW sweep starting after `d_ref`,
   *   EQUAL   if `d1 == d2`,
   *   LARGER  if `d1` comes after `d2`.
   *
   * PRECONDITION: neither `d1` nor `d2` is the identification point.
   */
  class Compare_x_1 {
  public:
    explicit Compare_x_1(const Circle_traits_1& tr) : m_tr(tr) {}

    Comparison_result operator()(const Point_1& d1, const Point_1& d2) const {
      // Precondition: neither argument is the identification point.
      CGAL_precondition(! m_tr.is_on_x_identification_1_object()(d1));
      CGAL_precondition(! m_tr.is_on_x_identification_1_object()(d2));

      // SMALLER: da comes before db in the CCW order starting after m_ref
      auto ccwib = m_tr.counterclockwise_in_between_2_object();
      if (ccwib(d1, m_tr.m_ref, d2)) return SMALLER;
      if (ccwib(d2, m_tr.m_ref, d1)) return LARGER;
      return EQUAL;
    }

  private:
    const Circle_traits_1& m_tr;
  };

  Compare_x_1 compare_x_1_object() const { return Compare_x_1(*this); }

private:
  Direction_2 m_ref;  // the identification-point direction
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
