// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1Funcs
 * computes the overlay of two input 1D arrangement objects `arr_a` and `arr_b`,
 * and sets the output arrangement `arr_r` to represent the overlaid
 * arrangement. The three arrangements can have different geometry traits and
 * different topology traits, as long as the side categories of their geometry
 * traits are the same. The geometry traits of the resulting arrangement is used
 * to construct the resulting arrangement. The point types
 * `GeometryTraitsA::Point_1` and `GeometryTraitsB::Point_1` of the input
 * arrangements must be convertible to the point type
 * `GeometryTraitsRes::Point_1` of the resulting arrangement.
 *
 * \pre `arr_r` is empty.  In particular, `arr_r` does not refer to either
 * `arr_a` or `arr_b` (that is, overlay in place is not supported).
 *
 * \tparam OverlayVisitor a model of the concept `OverlayVisitor`; its callback
 * functions update the vertices and edges of the resulting arrangement based
 * on the vertices and edges of the input arrangements that induce them.
 *
 * \sa `OverlayVisitor`
 */
template <typename GeometryTraitsA, typename GeometryTraitsB, typename GeometryTraitsRes,
          typename TopologyTraitsA, typename TopologyTraitsB, typename TopologyTraitsRes,
          bool BinarySearchA, bool BinarySearchB, bool BinarySearchRes,
          typename OverlayVisitor>
void overlay(const Arrangement_on_curve_1<GeometryTraitsA, TopologyTraitsA, BinarySearchA>& arr_a,
             const Arrangement_on_curve_1<GeometryTraitsB, TopologyTraitsB, BinarySearchB>& arr_b,
             Arrangement_on_curve_1<GeometryTraitsRes, TopologyTraitsRes, BinarySearchRes>& arr_r,
             OverlayVisitor& visitor);

/*! \ingroup PkgArrangementOnCurve1Funcs
 * computes the overlay of two input 1D arrangement objects `arr_a` and `arr_b`,
 * and sets the output arrangement `arr_r` to represent the overlaid
 * arrangement. The three arrangements can have different geometry traits and
 * different topology traits, as long as the side categories of their geometry
 * traits are the same. The geometry traits of the resulting arrangement is used
 * to construct the resulting arrangement. The point types
 * `GeometryTraitsA::Point_1` and `GeometryTraitsB::Point_1` of the input
 * arrangements must be convertible to the point type
 * `GeometryTraitsRes::Point_1` of the resulting arrangement.
 *
 * \pre `arr_r` is empty.  In particular, `arr_r` does not refer to either
 * `arr_a` or `arr_b` (that is, overlay in place is not supported).
 *
 * The function uses `Default_overlay_visitor`, so the data of the cells of the
 * resulting arrangement, if any, is not set.
 */
template <typename GeometryTraitsA, typename GeometryTraitsB, typename GeometryTraitsRes,
          typename TopologyTraitsA, typename TopologyTraitsB, typename TopologyTraitsRes,
          bool BinarySearchA, bool BinarySearchB, bool BinarySearchRes>
void overlay(const Arrangement_on_curve_1<GeometryTraitsA, TopologyTraitsA, BinarySearchA>& arr_a,
             const Arrangement_on_curve_1<GeometryTraitsB, TopologyTraitsB, BinarySearchB>& arr_b,
             Arrangement_on_curve_1<GeometryTraitsRes, TopologyTraitsRes, BinarySearchRes>& arr_r);

} // namespace Arrangement_on_curve_1
} // namespace CGAL
