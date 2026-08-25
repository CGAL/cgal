// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1Funcs
 *
 * inserts a point into an arrangement. The function first calls `locate(arr, p)`
 * internally to find where the point \f$p\f$ is located. If \f$p\f$ already
 * matches an existing vertex, the descriptor of that vertex is
 * returned. Otherwise, the function invokes `insert_empty()`,
 * `insert_before()`, `insert_after()`, or `split_edge()` as needed to safely
 * update the topology, and returns the newly created vertex.
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
