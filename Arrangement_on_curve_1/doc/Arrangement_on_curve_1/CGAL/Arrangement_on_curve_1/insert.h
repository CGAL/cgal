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
 *
 * inserts a point into an arrangement. The function first calls `locate(arr, p)`
 * internally to find where the point \f$p\f$ is located. If \f$p\f$ already
 * matches an existing vertex, the descriptor of that vertex is
 * returned. Otherwise, the function creates a new vertex and returns it.
 * Depending on where \f$p\f$ lies, it invokes one of the following:
 *
 * - \link ArrangementOnCurve_1::insert_empty() `insert_empty()`\endlink, if the
 *   arrangement is empty;
 * - \link CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1::insert_left_boundary_vertex()
 *   `insert_left_boundary_vertex()`\endlink
 *   or
 *   \link CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1::insert_right_boundary_vertex()
 *   `insert_right_boundary_vertex()`\endlink,
 *   if \f$p\f$ lies on a closed side;
 * - \link CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1::insert_identification_vertex()
 *   `insert_identification_vertex()`\endlink,
 *   if \f$p\f$ is the identification point;
 * - \link ArrangementOnCurve_1::insert_after() `insert_after()`\endlink, if the
 *   parameter space is identified and \f$p\f$ lies after the last vertex;
 * - \link ArrangementOnCurve_1::split_edge() `split_edge()`\endlink, otherwise.
 *
 * \param arr the arrangement
 * \param p the point to insert
 * \return the vertex descriptor that identifies the vertex associated with the `p`.
 */
template <typename GeometryTraits, typename TopologyTraits, bool BinarySearch>
typename Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>::Vertex_descriptor
insert(Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>& arr,
       const typename GeometryTraits::Point_1& p);

} // namespace Arrangement_on_curve_1
} // namespace CGAL
