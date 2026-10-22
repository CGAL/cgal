// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_SIDE_TAGS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_SIDE_TAGS_H

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1Tags
 *
 * This type tag is used to indicate that a side of the parameter space, either
 * left or right is open. More precisely, when this condition occurs at the left
 * or right side, then the corresponding categories, namely,
 * `Left_side_category` and `Right_side_category`, nested in every geometry
 * traits class, must be convertible to the type `Open_side_tag`.
 *
 * \sa `Closed_side_tag`
 * \sa `Identified_side_tag`
 * \sa `AocTraits_1`
 */
struct Open_side_tag {};

/*! \ingroup PkgArrangementOnCurve1Tags
 *
 * This type tag is used to indicate that a side of the parameter space, either
 * left or right is closed. More precisely, when this condition occurs at the left
 * or right side, then the corresponding categories, namely,
 * `Left_side_category` and `Right_side_category`, nested in every geometry
 * traits class, must be convertible to the type `Closed_side_tag`.
 *
 * \sa `Open_side_tag`
 * \sa `Identified_side_tag`
 * \sa `AocTraits_1`
 */
struct Closed_side_tag {};

/*! \ingroup PkgArrangementOnCurve1Tags
 *
 * This type tag is used to indicate that a side of the parameter space, either
 * left or right is identified. More precisely, when this condition occurs at
 * the left or right side, then the corresponding categories, namely,
 * `Left_side_category` and `Right_side_category`, nested in every geometry
 * traits class, must be convertible to the type `Identified_side_tag`.
 *
 * \sa `Open_side_tag`
 * \sa `Closed_side_tag`
 * \sa `AocTraits_1`
 */
struct Identified_side_tag {};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
