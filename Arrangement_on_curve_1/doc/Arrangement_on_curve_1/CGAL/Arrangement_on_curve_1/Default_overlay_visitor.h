// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1Ref
 *
 * \brief The class template `Default_overlay_visitor` is an overlay visitor
 * whose callback functions do nothing; the cells of the resulting arrangement
 * are left with default-constructed data, if any.
 *
 * It is the visitor used by the overload of `overlay()` that does not accept a
 * visitor.  It can also serve as a base class for a visitor that overrides
 * only some of the callback functions.
 *
 * \tparam ArrangementA the type of the first input arrangement, an instance of
 *   `Arrangement_on_curve_1`.
 * \tparam ArrangementB the type of the second input arrangement, an instance
 *   of `Arrangement_on_curve_1`.
 * \tparam ArrangementR the type of the resulting arrangement, an instance of
 *   `Arrangement_on_curve_1`.
 *
 * \cgalModels{OverlayVisitor}
 *
 * \sa `Copy_overlay_visitor`
 * \sa `overlay()`
 */
template <typename ArrangementA, typename ArrangementB, typename ArrangementR>
class Default_overlay_visitor {
public:
  //! \name Creation
  //! @{

  /*! constructs a visitor for overlaying `arr_a` and `arr_b` into `arr_r`.
   */
  Default_overlay_visitor(const ArrangementA& arr_a, const ArrangementB& arr_b, ArrangementR& arr_r);

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL
