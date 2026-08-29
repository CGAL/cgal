// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
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
 * \return a discriminated union container of type Const_`location_result` (an
 * instance of `std::variant` template) that identifies a cell (i.e., a vertex
 * or an edge). In particular, the returned object is a
 * `Const_vertex_descriptor` if the geometric embedding of the identified vertex
 * coincides exactly with the query point `q`, or otherwise an
 * `Const_edge_desriptor`, such that the geometric embedding of the identified
 * edge contains `q` in its interior. The search is implemented via a linear
 * topological walk from left to right starting at `unbounded_left_edge()`,
 * using the geometry traits functor to locate the precise cell containing the
 * query point parameter.
 */
template <typename GeometryTraits, typename TopologyTraits, bool BinarySearch>
typename Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>::Const_location_result
locate(Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>& arr,
       const typename GeometryTraits::Point_1& q);

/*! \ingroup PkgArrangementOnCurve1Funcs
 * locates a query point in a given arrangement.
 *
 * \param arr the arrangement
 * \param q the query point
 *
 * \return a discriminated union container of type `Location_result` (an
 * instance of `std::variant` template) that identifies a cell (i.e., a vertex
 * or an edge). In particular, the returned object is a `Vertex_descriptor` if
 * the geometric embedding of the identified vertex coincides exactly with the
 * query point `q`, or otherwise an `Edge_desriptor`, such that the geometric
 * embedding of the identified edge contains `q` in its interior. The search
 * is implemented via a linear topological walk from left to right starting at
 * `unbounded_left_edge()`, using the geometry traits functor to locate the
 * precise cell containing the query point parameter.
 */
template <typename GeometryTraits, typename TopologyTraits, bool BinarySearch>
typename Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>::Location_result
locate(Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>& arr,
       const typename GeometryTraits::Point_1& q);

} // namespace Arrangement_on_curve_1
} // namespace CGAL
