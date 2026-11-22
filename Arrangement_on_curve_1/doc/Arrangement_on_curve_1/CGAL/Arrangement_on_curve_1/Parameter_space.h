// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1Enums
 *
 * The scoped enumeration `Interval_parameter_space` indicates where the
 * pre-image \f$t\f$ of a point lies in the interval parameter space \f$X\f$ of
 * the master curve.  It is the return type of the functor
 * `AocTraits::ParameterSpaceInX_1`, which every model of the concept
 * `AocClosedBoundaryTraits_1` provides.
 *
 * \sa `AocTraits::ParameterSpaceInX_1`
 * \sa `AocClosedBoundaryTraits_1`
 */
enum class Interval_parameter_space {
  LEFT_BOUNDARY = 0, //!< \f$t\f$ lies on the left boundary of \f$X\f$.
  RIGHT_BOUNDARY,    //!< \f$t\f$ lies on the right boundary of \f$X\f$.
  INTERIOR           //!< \f$t\f$ lies in the interior of \f$X\f$.
};

/*! \ingroup PkgArrangementOnCurve1Enums
 *
 * An alias of the enumeration `Interval_parameter_space`.
 */
typedef Interval_parameter_space Parameter_space;

} // namespace Arrangement_on_curve_1
} // namespace CGAL
