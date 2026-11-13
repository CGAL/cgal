// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1TraitsClasses
 *
 * \brief `Left_ray_traits_1` is a model of the `AocClosedBoundaryTraits_1` concept.
 *
 * It parameterizes a 1D arrangement using a coordinate field type (\c Kernel::FT)
 * representing a continuous ray \f$[a, +\infty)\f$.
 * The left side of the parameter space is closed and the right side is open.
 *
 * \tparam Kernel must be a standard \cgal `Kernel` model (e.g., `Exact_predicates_exact_constructions_kernel`).
 *
 * \cgalModels{AocClosedBoundaryTraits_1}
 */
template <typename Kernel>
class Left_ray_traits_1 {
public:
  //! \name Creation
  //! @{

  /* construct a traits object from the boundary.
   * \param left_bound the left bound.
   */
  Left_ray_traits_1(const Point_1& left_bound);

  //! @}

  //! \name Types
  //! @{

  //! The left side of the interval parameter space is closed.
  typedef Closed_side_tag Left_side_category;

  //! The right side of the interval parameter space is open.
  typedef Open_side_tag Right_side_category;

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL
