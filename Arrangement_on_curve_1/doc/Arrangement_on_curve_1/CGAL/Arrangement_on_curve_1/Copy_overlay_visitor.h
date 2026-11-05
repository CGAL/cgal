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
 * \brief The class template `Copy_overlay_visitor` is an overlay visitor that
 * copies the user data of the input cells to the cells of the resulting
 * arrangement.
 *
 * The data of a resulting cell is copied from the input cell it originates
 * from as follows:
 *
 * | Callback             | The resulting cell is induced by                  | The data is copied from |
 * |----------------------|---------------------------------------------------|-------------------------|
 * | `create_vertex_vv()` | a vertex of each input arrangement                | the vertex of the first |
 * | `create_vertex_ve()` | a vertex of the first and an edge of the second   | the vertex              |
 * | `create_vertex_ev()` | an edge of the first and a vertex of the second   | the vertex              |
 * | `create_edge_ee()`   | an edge of each input arrangement                 | the edge of the first   |
 *
 * Data is copied only if the corresponding data type of the resulting
 * arrangement is not `void` and the data type of the source cell is
 * convertible to it; otherwise, the callback function does nothing.
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
 * \sa `Default_overlay_visitor`
 * \sa `overlay()`
 */
template <typename ArrangementA, typename ArrangementB, typename ArrangementR>
class Copy_overlay_visitor {
public:
  //! \name Creation
  //! @{

  /*! constructs a visitor for overlaying `arr_a` and `arr_b` into `arr_r`.
   * The visitor reads the data of `arr_a` and `arr_b` and writes the data of
   * `arr_r` through their data property maps.
   */
  Copy_overlay_visitor(const ArrangementA& arr_a, const ArrangementB& arr_b, ArrangementR& arr_r);

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL
