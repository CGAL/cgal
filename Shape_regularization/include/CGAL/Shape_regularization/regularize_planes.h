// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot, Florent Lafarge
//

#ifndef CGAL_SHAPE_REGULARIZATION_REGULARIZE_PLANES_H
#define CGAL_SHAPE_REGULARIZATION_REGULARIZE_PLANES_H

/// \cond SKIP_IN_MANUAL
#include <CGAL/license/Shape_regularization.h>
/// \endcond

/**
* \ingroup PkgShapeRegularizationRef
* \file CGAL/Shape_regularization/regularize_planes.h
* This header includes all classes for regularizing planes.
* It also includes the corresponding free functions.
*/

// Boost includes.
/// \cond SKIP_IN_MANUAL
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
/// \endcond

// Shape detection includes.
/// \cond SKIP_IN_MANUAL
#include <CGAL/Shape_detection/Efficient_RANSAC.h>
/// \endcond

// Internal includes.
/// \cond SKIP_IN_MANUAL
#include <CGAL/Shape_regularization/internal/utils.h>
/// \endcond

namespace CGAL {
namespace Shape_regularization {
namespace Planes {

  /// \cond SKIP_IN_MANUAL
  template<
  typename PointRange,
  typename PointMap,
  typename PlaneRange,
  typename PlaneMap,
  typename IndexMap,
  typename GeomTraits>
  void regularize_planes(
    const PointRange& points,
    const PointMap point_map,
    PlaneRange& planes,
    const PlaneMap plane_map,
    const IndexMap index_map,
    const GeomTraits&,
    const bool regularize_parallelism,
    const bool regularize_orthogonality,
    const bool regularize_coplanarity,
    const bool regularize_axis_symmetry,
    typename GeomTraits::FT tolerance_angle
      = typename GeomTraits::FT(25),
    const typename GeomTraits::FT tolerance_coplanarity
      = typename GeomTraits::FT(1) / typename GeomTraits::FT(100),
    const typename GeomTraits::Vector_3 symmetry_direction
      = typename GeomTraits::Vector_3(
        typename GeomTraits::FT(0),
        typename GeomTraits::FT(0),
        typename GeomTraits::FT(1))) {

    using Kernel = GeomTraits;
    using FT = typename Kernel::FT;
    using Point = typename Kernel::Point_3;
    using Vector = typename Kernel::Vector_3;
    using Plane = typename Kernel::Plane_3;

    using Plane_cluster = typename internal::Plane_cluster<Kernel>;

    // Compute centroids and areas.
    std::vector<Point> centroids;
    std::vector<FT> areas;
    internal::compute_centroids_and_areas<Kernel>(
      points, point_map, planes.size(), index_map, centroids, areas);

    tolerance_angle *= static_cast<FT>(CGAL_PI) / FT(180);
    const FT tolerance_cosangle = FT(FT(1) - std::cos(tolerance_angle));
    const FT tolerance_cosangle_ortho =
      FT(std::cos((FT(1) / FT(2)) * static_cast<FT>(CGAL_PI) - FT(tolerance_angle)));

    // Cluster the parallel primitives and store them in clusters,
    // compute the normal, size and cos angle to the symmetry direction of each cluster.
    std::vector<Plane_cluster> clusters;
    internal::compute_parallel_clusters<Kernel>(
      planes, plane_map, clusters, areas,
      (regularize_parallelism ? tolerance_cosangle : FT(0)),
      (regularize_axis_symmetry ? symmetry_direction : CGAL::NULL_VECTOR));

    if (regularize_orthogonality) {
      // Discover orthogonal relationships between clusters.
      for (std::size_t i = 0; i < clusters.size(); ++i) {
        for (std::size_t j = i + 1; j < clusters.size(); ++j) {
          if (CGAL::abs(clusters[i].normal * clusters[j].normal) < tolerance_cosangle_ortho) {

            clusters[i].orthogonal_clusters.push_back(j);
            clusters[j].orthogonal_clusters.push_back(i);
          }
        }
      }
    }

    if (regularize_axis_symmetry) {
      // Cluster the symmetry cos angle and store their centroids in
      // cosangle_centroids and the centroid index of each cluster in
      // list_cluster_index.
      internal::cluster_symmetric_cosangles<Kernel>(
        clusters, tolerance_cosangle, tolerance_cosangle_ortho);
    }

    // Find subgraphs of mutually orthogonal clusters (store indices of
    // clusters in subgraph_clusters), and select the cluster of the largest area.
    if (regularize_orthogonality || regularize_axis_symmetry) {
      internal::subgraph_mutually_orthogonal_clusters<Kernel>(
        clusters, (regularize_axis_symmetry ? symmetry_direction : CGAL::NULL_VECTOR));
    }

    // Recompute optimal plane for each primitive after the normal regularization.
    for (std::size_t i = 0; i < clusters.size(); ++i) {
      Vector vec_reg = clusters[i].normal;
      for (std::size_t j = 0; j < clusters[i].planes.size(); ++j) {
        const std::size_t index_prim = clusters[i].planes[j];
        const Plane& plane = get(plane_map, *(planes.begin() + index_prim));

        const Point pt_reg = plane.projection(centroids[index_prim]);
        if (plane.orthogonal_vector() * vec_reg < FT(0)) {
          vec_reg = -vec_reg;
        }
        const Plane plane_reg(pt_reg, vec_reg);

        if (CGAL::abs(plane.orthogonal_vector() * vec_reg) > FT(1) - tolerance_cosangle) {
          put(plane_map, *(planes.begin() + index_prim), plane_reg);
        }
      }
    }

    if (regularize_coplanarity) {
      // Detect coplanarity and use list_coplanar_prim to store the results.
      for (std::size_t i = 0; i < clusters.size(); ++i) {
        Vector vec_reg = clusters[i].normal;
        for (std::size_t ip = 0; ip < clusters[i].planes.size(); ++ip) {
          clusters[i].coplanar_group.push_back(static_cast<std::size_t>(-1));
        }

        std::size_t cop_index = 0;
        for (std::size_t j = 0; j < clusters[i].planes.size(); ++j) {
          const std::size_t index_prim = clusters[i].planes[j];

          if (clusters[i].coplanar_group[j] == static_cast<std::size_t>(-1)) {
            const Plane& plane = get(plane_map, *(planes.begin() + index_prim));
            clusters[i].coplanar_group[j] = cop_index;

            const Point pt_reg = plane.projection(centroids[index_prim]);
            const Plane plan_reg(pt_reg, vec_reg);

            for (std::size_t k = j + 1; k < clusters[i].planes.size(); ++k) {
              if (clusters[i].coplanar_group[k] == static_cast<std::size_t>(-1)) {
                const std::size_t index_prim_next = clusters[i].planes[k];
                const Plane& plane_next = get(plane_map, *(planes.begin() + index_prim_next));
                const Point pt_reg_next = plane_next.projection(centroids[index_prim_next]);
                const Point pt_proj = plan_reg.projection(pt_reg_next);
                const FT distance = CGAL::sqrt(CGAL::squared_distance(pt_reg_next, pt_proj));

                if (distance < tolerance_coplanarity) {
                  clusters[i].coplanar_group[k] = cop_index;
                }
              }
            }
            ++cop_index;
          }
        }

        // Regularize primitive positions by computing barycenter of the coplanar planes.
        std::vector<Point> pt_bary(cop_index, Point(FT(0), FT(0), FT(0)));
        std::vector<FT> area(cop_index, FT(0));
        for (std::size_t j = 0; j < clusters[i].planes.size (); ++j) {
          const std::size_t index_prim = clusters[i].planes[j];
          const std::size_t group = clusters[i].coplanar_group[j];

          const Point pt_reg = get(plane_map,
            *(planes.begin() + index_prim)).projection(centroids[index_prim]);

          pt_bary[group] = CGAL::barycenter(
            pt_bary[group], area[group], pt_reg, areas[index_prim]);
          area[group] += areas[index_prim];
        }

        for (std::size_t j = 0; j < clusters[i].planes.size (); ++j) {
          const std::size_t index_prim = clusters[i].planes[j];
          const std::size_t group = clusters[i].coplanar_group[j];
          const Plane plane_reg(pt_bary[group], vec_reg);

          if (get(plane_map,
          *(planes.begin() + index_prim)).orthogonal_vector()
          * plane_reg.orthogonal_vector() < 0) {
            put(plane_map, *(planes.begin() + index_prim), plane_reg.opposite());
          } else {
            put(plane_map, *(planes.begin() + index_prim), plane_reg);
          }
        }
      }
    }
  }

  // Workaround for the bug reported here:
  // https://developercommunity.visualstudio.com/content/problem/340310/unaccepted-typename-that-other-compilers-require.html
  #if _MSC_VER == 1915
  #define CGAL_TYPENAME_FOR_MSC
  #else
  #define CGAL_TYPENAME_FOR_MSC typename
  #endif

  // This variant deduces the kernel from the point property map.
  template<
  typename PointRange,
  typename PointMap,
  typename PlaneRange,
  typename PlaneMap,
  typename IndexMap>
  void regularize_planes(
    const PointRange& points,
    const PointMap point_map,
    PlaneRange& planes,
    const PlaneMap plane_map,
    const IndexMap index_map,
    const bool regularize_parallelism,
    const bool regularize_orthogonality,
    const bool regularize_coplanarity,
    const bool regularize_axis_symmetry,
    const typename Kernel_traits<
      typename boost::property_traits<PointMap>::value_type>::Kernel::FT tolerance_angle =
      CGAL_TYPENAME_FOR_MSC Kernel_traits<
        typename boost::property_traits<PointMap>::value_type>::Kernel::FT(25),
    const typename Kernel_traits<
      typename boost::property_traits<PointMap>::value_type>::Kernel::FT tolerance_coplanarity =
      CGAL_TYPENAME_FOR_MSC Kernel_traits<
        typename boost::property_traits<PointMap>::value_type>::Kernel::FT(1) /
      CGAL_TYPENAME_FOR_MSC Kernel_traits<
        typename boost::property_traits<PointMap>::value_type>::Kernel::FT(100),
    const typename Kernel_traits<
      typename boost::property_traits<PointMap>::value_type>::Kernel::Vector_3 symmetry_direction
        = CGAL_TYPENAME_FOR_MSC Kernel_traits<
          typename boost::property_traits<PointMap>::value_type>::Kernel::Vector_3(
            CGAL_TYPENAME_FOR_MSC Kernel_traits<
              typename boost::property_traits<PointMap>::value_type>::Kernel::FT(0),
            CGAL_TYPENAME_FOR_MSC Kernel_traits<
              typename boost::property_traits<PointMap>::value_type>::Kernel::FT(0),
            CGAL_TYPENAME_FOR_MSC Kernel_traits<
              typename boost::property_traits<PointMap>::value_type>::Kernel::FT(1))) {

    typedef typename boost::property_traits<PointMap>::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;

    Kernel kernel;
    regularize_planes(
      points, point_map, planes, plane_map, index_map, kernel,
      regularize_parallelism, regularize_orthogonality,
      regularize_coplanarity, regularize_axis_symmetry,
      tolerance_angle, tolerance_coplanarity, symmetry_direction);
  }

  #ifdef CGAL_TYPENAME_FOR_MSC
  #undef CGAL_TYPENAME_FOR_MSC
  #endif

  template<
  typename PlaneRange,
  typename PlaneMap,
  typename PointRange,
  typename PointMap,
  typename NamedParameters = parameters::Default_named_parameters>
  void regularize_planes(
    PlaneRange& planes,
    const PlaneMap plane_map,
    const PointRange& points,
    const PointMap point_map,
    const NamedParameters& np = parameters::default_values()) {

    using parameters::get_parameter;
    using parameters::choose_parameter;
    using parameters::is_default_parameter;

    using PlaneIndexMap = typename CGAL::Point_set_processing_3::GetPlaneIndexMap<NamedParameters>::type;

    static_assert(!is_default_parameter<NamedParameters, internal_np::plane_index_t>::value,
                              "Error: no plane index map");

    const PlaneIndexMap index_map =
      choose_parameter(get_parameter(np, internal_np::plane_index_map), PlaneIndexMap());

    using Kernel = typename Kernel_traits<
      typename boost::property_traits<PointMap>::value_type>::type;

    using FT = typename Kernel::FT;
    using Vector_3 = typename Kernel::Vector_3;

    const bool reg_prll = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::regularize_parallelism),
      true);
    const bool reg_orth = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::regularize_orthogonality),
      true);
    const bool reg_copl = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::regularize_coplanarity),
      true);
    const bool reg_symm = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::regularize_axis_symmetry),
      true);

    const FT tol_angle = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::maximum_angle),
      FT(25));
    const FT tol_copln = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::maximum_offset),
      FT(1) / FT(100));
    const Vector_3 sym_dir = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::symmetry_direction),
      Vector_3(FT(0), FT(0), FT(1)));

    Kernel kernel;
    regularize_planes(
      points, point_map, planes, plane_map, index_map, kernel,
      reg_prll, reg_orth, reg_copl, reg_symm,
      tol_angle, tol_copln, sym_dir);
  }
  /// \endcond

  /*!
    \ingroup PkgShapeRegularizationRefPlanes

    \brief Hierarchical plane regularization.

    Given a set of detected planes with their corresponding inlier sets,
    this function enables to reinforce four types of regularities among these planes:
    - *Parallelism*: planes, which are detected as near parallel, are made exactly parallel.
    - *Orthogonality*: planes, which are detected as near orthogonal, are made exactly orthogonal.
    - *Coplanarity*: parallel planes, which are detected as near coplanar, are made exactly coplanar.
    - *Axis-Symmetry*: planes, which are detected as near symmetrical with respect to a user-specified axis,
    are made exactly symmetrical.

    %Planes are directly modified. Points are left unaltered, as well as their
    relationship to the planes (no transfer of a point from a primitive plane to another).

    This function infers a traits class `GeomTraits` from the `PointRange` value type.

    The implementation follows \cgalCite{cgal:vla-lod-15}.

    \tparam PlaneRange
    a model of `Range` whose iterator type is `RandomAccessIterator`

    \tparam PointRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam NamedParameters
    a sequence of \ref bgl_namedparameters "Named Parameters"

    \param planes
    a range of planes to be regularized

    \param points
    a const range of points assigned to planes, see the `plane_index_map` below
    for more details

    \param np
    a sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:

    \cgalNamedParamsBegin
      \cgalParamNBegin{plane_map}
        \cgalParamDescription{a property map that maps a plane from `planes` range to `CGAL::Plane_3<GeomTraits>`}
        \cgalParamType{a model of `WritablePropertyMap` with the value type `CGAL::Plane_3<GeomTraits>`}
        \cgalParamDefault{`PlaneMap()`}
      \cgalParamNEnd
      \cgalParamNBegin{point_map}
        \cgalParamDescription{a property map that maps a point from `points` range to `CGAL::Point_3<GeomTraits>`}
        \cgalParamType{a model of `ReadablePropertyMap` with the value type `CGAL::Point_3<GeomTraits>`}
        \cgalParamDefault{`PointMap()`}
      \cgalParamNEnd
      \cgalParamNBegin{plane_index_map}
        \cgalParamDescription{a property map that associates the index of a point
          in the `points` range to the index of a plane in the `planes` range (-1 if
          point is not assigned to a plane)}
        \cgalParamType{a model of `ReadablePropertyMap` with `std::size_t` as key type and `int` as value type}
        \cgalParamDefault{There is no default, this parameters is mandatory}
      \cgalParamNEnd
      \cgalParamNBegin{maximum_angle}
        \cgalParamDescription{maximum allowed angle in degrees between plane normals used
          for parallelism, orthogonality, and axis symmetry}
        \cgalParamType{`GeomTraits::FT`}
        \cgalParamDefault{25 degrees}
      \cgalParamNEnd
      \cgalParamNBegin{maximum_offset}
        \cgalParamDescription{maximum allowed orthogonal distance between two parallel planes
          such that they are considered to be coplanar}
        \cgalParamType{`GeomTraits::FT`}
        \cgalParamDefault{0.01 unit length}
      \cgalParamNEnd
      \cgalParamNBegin{regularize_parallelism}
        \cgalParamDescription{indicates whether parallelism should be regularized or not}
        \cgalParamType{boolean}
        \cgalParamDefault{true}
      \cgalParamNEnd
      \cgalParamNBegin{regularize_orthogonality}
        \cgalParamDescription{indicates whether orthogonality should be regularized or not}
        \cgalParamType{boolean}
        \cgalParamDefault{true}
      \cgalParamNEnd
      \cgalParamNBegin{regularize_coplanarity}
        \cgalParamDescription{indicates whether coplanarity should be regularized or not}
        \cgalParamType{boolean}
        \cgalParamDefault{true}
      \cgalParamNEnd
      \cgalParamNBegin{regularize_axis_symmetry}
        \cgalParamDescription{indicates whether axis symmetry should be regularized or not}
        \cgalParamType{boolean}
        \cgalParamDefault{true}
      \cgalParamNEnd
      \cgalParamNBegin{symmetry_direction}
        \cgalParamDescription{an axis for symmetry regularization}
        \cgalParamType{`GeomTraits::Vector_3`}
        \cgalParamDefault{Z axis that is `GeomTraits::Vector_3(0, 0, 1)`}
      \cgalParamNEnd
    \cgalNamedParamsEnd
  */
  template<
  typename PlaneRange,
  typename PointRange,
  typename NamedParameters>
  void regularize_planes(
    PlaneRange& planes,
    const PointRange& points,
    const NamedParameters& np) {

    using parameters::get_parameter;
    using parameters::choose_parameter;

    using PlaneMap = typename CGAL::Point_set_processing_3::
      GetPlaneMap<PlaneRange, NamedParameters>::type;
    const PlaneMap plane_map =
      choose_parameter(get_parameter(np, internal_np::plane_map), PlaneMap());

    using PointMap = typename CGAL::
      GetPointMap<PointRange, NamedParameters>::type;
    const PointMap point_map =
      choose_parameter(get_parameter(np, internal_np::point_map), PointMap());

    regularize_planes(planes, plane_map, points, point_map, np);
  }

} // namespace Planes
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_REGULARIZE_PLANES_H
