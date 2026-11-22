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
 * locates a query point in a given arrangement.
 *
 * \param arr the arrangement
 * \param q the query point
 *
 * \return a discriminated union of type
 * `ArrangementOnCurve_1::Location_result` (an instance of the `std::variant`
 * template) that identifies a cell, that is, a vertex or an edge.  In
 * particular, the returned object is a `ArrangementOnCurve_1::Vertex_descriptor`
 * if the point associated with the identified vertex coincides with the query
 * point `q`, or otherwise an `ArrangementOnCurve_1::Edge_descriptor` of the edge
 * that contains `q` in its interior.
 *
 * If `BinarySearch` is `true`, the search is a binary search over the vertices,
 * which takes \f$O(\log n)\f$ time.  Otherwise, it is a walk from left to right
 * that starts at
 * \link ArrangementOnCurve_1::left_boundary_edge() `left_boundary_edge()`\endlink,
 * which takes \f$O(n)\f$ time.
 */
template <typename GeometryTraits, typename TopologyTraits, bool BinarySearch>
typename Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>::Location_result
locate(Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>& arr,
       const typename GeometryTraits::Point_1& q);

/*! \ingroup PkgArrangementOnCurve1Funcs
 * locates a query point in a given arrangement (const version).
 *
 * \param arr the arrangement
 * \param q the query point
 *
 * \return a discriminated union of type
 * `ArrangementOnCurve_1::Const_location_result` (an instance of the
 * `std::variant` template) that identifies a cell, that is, a vertex or an
 * edge.  In particular, the returned object is a
 * `ArrangementOnCurve_1::Vertex_const_descriptor` if the point associated with
 * the identified vertex coincides with the query point `q`, or otherwise an
 * `ArrangementOnCurve_1::Edge_const_descriptor` of the edge that contains `q`
 * in its interior.
 */
template <typename GeometryTraits, typename TopologyTraits, bool BinarySearch>
typename Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>::Const_location_result
locate(const Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>& arr,
       const typename GeometryTraits::Point_1& q);

} // namespace Arrangement_on_curve_1
} // namespace CGAL
