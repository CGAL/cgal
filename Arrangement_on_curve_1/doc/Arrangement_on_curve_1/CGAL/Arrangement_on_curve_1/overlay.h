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
 * arrangement. All overlay template functions can be instantiated with
 * different geometric traits instances and different topology traits
 * instances. The geometry traits of the resulting arrangement is used to
 * construct the resulting arrangement. The type `GeometryTraitsRes::Point_2` of
 * both input arrangements must be convertible to the point type in the
 * resulting arrangement.
 *
 * \pre `arr_r` does not refer to either `arr_a` or `arr_b` (that is,
 * overlay in place is not supported).
 *
 * \pre The overlay-visitor object `visitor` must model the `OverlayVisitor`
 * concept; this object consists of callback functions that update the vertices
 * and edges of the resulting arrangement based of the vertices and edges of the
 * input arrangements that indice them.
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
 * arrangement. All overlay template functions can be instantiated with
 * different geometric traits instances and different topology traits
 * instances. The geometry traits of the resulting arrangement is used to
 * construct the resulting arrangement. The type `GeometryTraitsRes::Point_2` of
 * both input arrangements must be convertible to the point type in the
 * resulting arrangement.
 *
 * \pre `arr_r` does not refer to either `arr_a` or `arr_b` (that is,
 * "self overlay" (or overlay in place) is not supported).
 */
template <typename GeometryTraitsA, typename GeometryTraitsB, typename GeometryTraitsRes,
          typename TopologyTraitsA, typename TopologyTraitsB, typename TopologyTraitsRes,
          bool BinarySearchA, bool BinarySearchB, bool BinarySearchRes>
void overlay(const Arrangement_on_curve_1<GeometryTraitsA, TopologyTraitsA, BinarySearchA>& arr_a,
             const Arrangement_on_curve_1<GeometryTraitsB, TopologyTraitsB, BinarySearchB>& arr_b,
             Arrangement_on_curve_1<GeometryTraitsRes, TopologyTraitsRes, BinarySearchRes>& arr_r);

} // namespace Arrangement_on_curve_1
} // namespace CGAL
