// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/** \ingroup PkgArrangementOnCurve1TraitsClasses
 *
 * \brief `Circle_traits_1` is a geometry traits model of the
 * `AocIdentifiedBoundaryTraits_1` concept for 1D arrangements of points on the
 * unit circle \f$\mathbb{S}^1\f$.
 *
 * Points are represented as `Kernel::Direction_2` objects.  Using
 * `Direction_2` rather than a floating-point angle gives exact arithmetic when
 * the kernel is exact (e.g., `CGAL::Exact_predicates_exact_constructions_kernel`).
 *
 * ### Parameter space and boundary categories
 *
 * The parameter space of \f$\mathbb{S}^1\f$ is
 * \f$[\theta_{\mathrm{ref}}, \theta_{\mathrm{ref}} + 2\pi)\f$ with the two endpoints
 * identified, where \f$\theta_{\mathrm{ref}}\f$ is the angle of the reference
 * direction \f$d_{\mathrm{ref}}\f$ described below.  By default,
 * \f$\theta_{\mathrm{ref}} = 0\f$; the reference direction can be overridden through
 * the constructor.  Both sides of the parameter space therefore carry
 * `Identified_side_tag`:
 * \code
 *   Left_side_category = Identified_side_tag
 *   Right_side_category = Identified_side_tag
 * \endcode
 * `Circle_traits_1` must be paired with `Identified_topology_traits`.
 *
 * ### Reference direction and total order
 *
 * An optional reference direction \f$d_{\mathrm{ref}}\f$ (default: positive
 * `x`-axis \f$(1,0)\f$) defines the identification point.  The total order on
 * directions is the counter-clockwise (CCW) order starting strictly after
 * \f$d_{\mathrm{ref}}\f$: \f$d_{\mathrm{ref}} < d_1 < d_2 < \ldots < d_n <
 * d_{\mathrm{ref}}\f$ (one full turn).
 *
 * A direction equal to \f$d_{\mathrm{ref}}\f$ is the identification point:
 * `AocIdentifiedBoundaryTraits_1::Is_on_x_identification_1` returns `true` for
 * it, and `AocTraits_1::Compare_x_1` must not be called with it as either
 * argument.
 *
 * ### Exact angular comparison
 *
 * The comparison is implemented using only the kernel predicate
 * `Kernel::CompareAngleWithXAxis_2` (no trigonometric computation), which
 * compares the angles in \f$[0, 2\pi)\f$ between the positive \f$x\f$-axis and
 * two directions.  If \f$d_{\mathrm{ref}}\f$ is the positive \f$x\f$-axis, the
 * order of two directions \f$d_1\f$ and \f$d_2\f$ is the order of their
 * angles.  Otherwise, the directions whose angles are larger than the angle of
 * \f$d_{\mathrm{ref}}\f$ precede the others, and within each of these two
 * groups, the order is the order of the angles.
 *
 * \tparam Kernel must be a standard \cgal `Kernel` model (e.g.,
 *   `Exact_predicates_exact_constructions_kernel`) providing
 *   `Direction_2` and the `Compare_angle_with_x_axis_2` predicate.
 *
 * \cgalModels{AocIdentifiedBoundaryTraits_1}
 */
template <typename Kernel>
class Circle_traits_1 {
public:
  //! \name Types
  //! @{

  //! The 1D arrangement point type is a direction in the plane.
  typedef typename Kernel::Direction_2 Point_1;

  //! The left side of the parameter space is identified.
  typedef Identified_side_tag Left_side_category;

  //! The right side of the parameter space is identified.
  typedef Identified_side_tag Right_side_category;

  //! @}

  //! \name Creation
  //! @{

  /*! constructs from a reference direction defining the identification point.
   *
   * An optional reference direction sets the position of the identification
   * point.  The reference direction is the unique point where the two
   * identified sides meet.  The default identification point is the positive
   * `x`-axis direction \f$(1, 0)\f$.
   *
   * \param reference the identification-point.
   */
  explicit Circle_traits_1(const Point_1& reference = Point_1(1, 0));

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL
