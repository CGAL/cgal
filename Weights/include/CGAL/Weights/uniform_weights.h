// Copyright (c) 2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_UNIFORM_WEIGHTS_H
#define CGAL_UNIFORM_WEIGHTS_H

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

#include <boost/graph/graph_traits.hpp>

namespace CGAL {
namespace Weights {

// 2D ==============================================================================================

/*!
  \ingroup PkgWeightsRefUniformWeights
  \brief returns `1`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`
*/
template<typename GeomTraits>
typename GeomTraits::FT uniform_weight(const typename GeomTraits::Point_2&,
                                       const typename GeomTraits::Point_2&,
                                       const typename GeomTraits::Point_2&,
                                       const typename GeomTraits::Point_2&,
                                       const GeomTraits&)
{
  return {1};
}

/*!
  \ingroup PkgWeightsRefUniformWeights
  \brief returns `1`.
  \tparam Kernel a model of `Kernel`
*/
template<typename GeomTraits>
typename GeomTraits::FT uniform_weight(const CGAL::Point_2<GeomTraits>& p0,
                                       const CGAL::Point_2<GeomTraits>& p1,
                                       const CGAL::Point_2<GeomTraits>& p2,
                                       const CGAL::Point_2<GeomTraits>& q)
{
  const GeomTraits traits;
  return uniform_weight(p0, p1, p2, q, traits);
}

// 3D ==============================================================================================

/*!
  \ingroup PkgWeightsRefUniformWeights
  \brief returns `1`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`
*/
template<typename GeomTraits>
typename GeomTraits::FT uniform_weight(const typename GeomTraits::Point_3&,
                                       const typename GeomTraits::Point_3&,
                                       const typename GeomTraits::Point_3&,
                                       const typename GeomTraits::Point_3&,
                                       const GeomTraits&)
{
  return {1};
}

/*!
  \ingroup PkgWeightsRefUniformWeights
  \brief returns `1`.
  \tparam Kernel a model of `Kernel`
*/
template<typename GeomTraits>
typename GeomTraits::FT uniform_weight(const CGAL::Point_3<GeomTraits>& p0,
                                       const CGAL::Point_3<GeomTraits>& p1,
                                       const CGAL::Point_3<GeomTraits>& p2,
                                       const CGAL::Point_3<GeomTraits>& q)
{
  const GeomTraits traits;
  return uniform_weight(p0, p1, p2, q, traits);
}

/// \cond SKIP_IN_MANUAL

// Undocumented uniform weight class taking as input a polygon mesh.
//
// It is currently used in:
// Polygon_mesh_processing -> triangulate_hole_Polyhedron_3_test.cpp
// Polygon_mesh_processing -> triangulate_hole_Polyhedron_3_no_delaunay_test.cpp
// CGAL Lab -> Fairing_plugin.cpp
// CGAL Lab -> Hole_filling_plugin.cpp
template<class PolygonMesh>
class Uniform_weight
{
public:
  using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  double w_i(vertex_descriptor) { return 1.; }
  double w_ij(halfedge_descriptor) { return 1.; }
};

/// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_UNIFORM_WEIGHTS_H
