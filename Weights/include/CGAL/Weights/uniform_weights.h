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

// Internal includes.
#include <CGAL/Weights/internal/utils.h>

namespace CGAL {
namespace Weights {

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightsRefUniformWeights

    \brief this function always returns 1, given four points in 2D and a traits class
    with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT uniform_weight(
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const GeomTraits&) { }

  /*!
    \ingroup PkgWeightsRefUniformWeights

    \brief this function always returns 1, given four points in 3D and a traits class
    with geometric objects, predicates, and constructions.
  */
  template<typename GeomTraits>
  typename GeomTraits::FT uniform_weight(
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const GeomTraits&) { }

  /*!
    \ingroup PkgWeightsRefUniformWeights

    \brief this function always returns 1, given four points in 2D which are
    parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT uniform_weight(
    const CGAL::Point_2<K>&,
    const CGAL::Point_2<K>&,
    const CGAL::Point_2<K>&,
    const CGAL::Point_2<K>&) { }

  /*!
    \ingroup PkgWeightsRefUniformWeights

    \brief this function always returns 1, given four points in 3D which are
    parameterized by a `Kernel` K.
  */
  template<typename K>
  typename K::FT uniform_weight(
    const CGAL::Point_3<K>&,
    const CGAL::Point_3<K>&,
    const CGAL::Point_3<K>&,
    const CGAL::Point_3<K>&) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  typename GeomTraits::FT uniform_weight(
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const GeomTraits&) {

    using FT = typename GeomTraits::FT;
    return FT(1);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT uniform_weight(
    const CGAL::Point_2<GeomTraits>& q,
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p) {

    const GeomTraits traits;
    return uniform_weight(q, t, r, p, traits);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT uniform_weight(
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const GeomTraits&) {

    using FT = typename GeomTraits::FT;
    return FT(1);
  }

  template<typename GeomTraits>
  typename GeomTraits::FT uniform_weight(
    const CGAL::Point_3<GeomTraits>& q,
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p) {

    const GeomTraits traits;
    return uniform_weight(q, t, r, p, traits);
  }

  // Undocumented uniform weight class taking as input a polygon mesh.
  // It is currently used in:
  // Polygon_mesh_processing -> triangulate_hole_Polyhedron_3_test.cpp
  // Polygon_mesh_processing -> triangulate_hole_Polyhedron_3_no_delaunay_test.cpp
  // Polyhedron demo -> Fairing_plugin.cpp
  // Polyhedron demo -> Hole_filling_plugin.cpp
  template<class PolygonMesh>
  class Uniform_weight {

  public:
    using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
    double w_i(vertex_descriptor) { return 1; }
    double w_ij(halfedge_descriptor) { return 1; }
  };

  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_UNIFORM_WEIGHTS_H
