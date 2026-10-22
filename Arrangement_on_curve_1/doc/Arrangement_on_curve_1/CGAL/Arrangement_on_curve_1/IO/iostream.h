// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1op_left_shift
 *
 * inserts the arrangement object `arr` into the output stream `os`. Only the
 * basic geometric and topological features of the arrangement are
 * inserted. Auxiliary data that may be attached to the vertices and edges is
 * ignored.
 *
 * The output consists of the number of vertices followed by their points in
 * increasing order along the curve, each on its own line.  If the parameter
 * space is identified, the identification point, if present, is last.  The
 * points are written with the output operator of `GeometryTraits_1::Point_1`.
 */
template <typename GeometryTraits_1, typename TopologyTraits, bool BinarySearch>
std::ostream& operator<<(std::ostream& os,
                         const Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>& arr);

/*! \ingroup PkgArrangementOnCurve1op_right_shift
 *
 * extracts an arrangement from a given input stream `is` into `arr`,
 * replacing its content. Only the basic geometric and topological features of
 * the arrangement are read and no auxiliary data is attached to the vertices
 * and edges.
 *
 * The input must have the format written by `operator<<()`; the points are
 * read with the input operator of `GeometryTraits_1::Point_1`.  If the input
 * is malformed, that is, the points are not in strictly increasing order, a
 * point on a closed boundary side is not first (left side) or last (right
 * side), or the identification point is not last, the `failbit` of `is` is set
 * and `arr` is left empty.
 */
template <typename GeometryTraits_1, typename TopologyTraits, bool BinarySearch>
std::istream& operator>>(std::istream& is,
                         Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>& arr);

} // namespace Arrangement_on_curve_1
} // namespace CGAL
