// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Jocelyn Meyron and Quentin Mérigot
//

#ifndef CGAL_VCM_ESTIMATE_NORMALS_H
#define CGAL_VCM_ESTIMATE_NORMALS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Point_set_processing_3/internal/Voronoi_covariance_3/voronoi_covariance_3.h>

#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Fuzzy_sphere.h>

#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Default_diagonalize_traits.h>

#include <iterator>
#include <vector>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace internal {

/// @cond SKIP_IN_MANUAL
/// Computes the VCM for each point in the property map.
/// The matrix is computed by intersecting the Voronoi cell
/// of a point and a sphere whose radius is `offset_radius` and discretized
/// by `N` planes.
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointMap is a model of `ReadablePropertyMap` with a value_type = `Kernel::Point_3`.
/// @tparam K Geometric traits class.
/// @tparam Covariance Covariance matrix type. It is similar to an array with a length of 6.
template < typename ForwardIterator,
           typename PointMap,
           class K,
           class Covariance
>
void
vcm_offset (ForwardIterator first, ///< iterator over the first input point.
            ForwardIterator beyond, ///< past-the-end iterator over the input points.
            PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3.
            std::vector<Covariance> &cov, ///< vector of covariance matrices.
            double offset_radius, ///< radius of the sphere.
            std::size_t N, ///< number of planes used to discretize the sphere.
            const K & /*kernel*/) ///< geometric traits.
{
    // Sphere discretization
    typename CGAL::Voronoi_covariance_3::Sphere_discretization<K> sphere(offset_radius, N);

    // Compute the Delaunay Triangulation
    std::vector<typename K::Point_3> points;
    points.reserve(std::distance(first, beyond));
    for (ForwardIterator it = first; it != beyond; ++it)
      points.push_back(get(point_map, *it));

    typedef Delaunay_triangulation_3<K> DT;
    DT dt(points.begin(), points.end());

    cov.clear();
    cov.reserve(points.size());
    // Compute the VCM
    for (typename std::vector<typename K::Point_3>::iterator
          it = points.begin(); it != points.end(); ++it)
    {
        typename DT::Vertex_handle vh = dt.nearest_vertex(*it);
        cov.push_back(
          Voronoi_covariance_3::voronoi_covariance_3(dt, vh, sphere)
        );
    }
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// Convolve using a radius.
template < class ForwardIterator,
           class PointMap,
           class K,
           class Covariance
>
void
vcm_convolve (ForwardIterator first,
              ForwardIterator beyond,
              PointMap point_map,
              const std::vector<Covariance> &cov,
              std::vector<Covariance> &ncov,
              double convolution_radius,
              const K &)
{
    typedef std::pair<typename K::Point_3, std::size_t>              Tree_point;
    typedef First_of_pair_property_map< Tree_point >                  Tree_map;
    typedef Search_traits_3<K>                                      Traits_base;
    typedef Search_traits_adapter<Tree_point, Tree_map, Traits_base>    Traits;
    typedef Kd_tree<Traits>                                                Tree;
    typedef Fuzzy_sphere<Traits>                                   Fuzzy_sphere;

    // Kd tree
    Tree tree;
    tree.reserve(cov.size());
    std::size_t i=0;
    for (ForwardIterator it = first; it != beyond; ++it, ++i)
        tree.insert( Tree_point(get(point_map, *it), i) );

    // Convolving
    ncov.clear();
    ncov.reserve(cov.size());
    for (ForwardIterator it = first; it != beyond; ++it) {
        std::vector<Tree_point> nn;
        tree.search(std::back_inserter(nn),
                    Fuzzy_sphere (get(point_map, *it), convolution_radius));

        Covariance m;
        std::fill(m.begin(), m.end(), typename K::FT(0));
        for (std::size_t k = 0; k < nn.size(); ++k)
        {
          std::size_t index = nn[k].second;
          for (int i=0; i<6; ++i)
            m[i] += cov[index][i];
        }
        ncov.push_back(m);
    }
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// Convolve using neighbors.
template < class ForwardIterator,
           class PointMap,
           class K,
           class Covariance
>
void
vcm_convolve (ForwardIterator first,
              ForwardIterator beyond,
              PointMap point_map,
              const std::vector<Covariance> &cov,
              std::vector<Covariance> &ncov,
              unsigned int nb_neighbors_convolve,
              const K &)
{
    typedef std::pair<typename K::Point_3, std::size_t>              Tree_point;
    typedef First_of_pair_property_map< Tree_point >                  Tree_map;
    typedef Search_traits_3<K>                                      Traits_base;
    typedef Search_traits_adapter<Tree_point, Tree_map, Traits_base>    Traits;
    typedef Orthogonal_k_neighbor_search<Traits>                Neighbor_search;
    typedef typename Neighbor_search::Tree                                 Tree;

    // Search tree
    Tree tree;
    tree.reserve(cov.size());
    std::size_t i=0;
    for (ForwardIterator it = first; it != beyond; ++it, ++i)
        tree.insert( Tree_point(get(point_map, *it), i) );

    // Convolving
    ncov.clear();
    ncov.reserve(cov.size());
    for (ForwardIterator it = first; it != beyond; ++it) {
        Neighbor_search search(tree, get(point_map, *it), nb_neighbors_convolve);

        Covariance m;
        for (typename Neighbor_search::iterator nit = search.begin();
             nit != search.end();
             ++nit)
        {
          std::size_t index = nit->first.second;
          for (int i=0; i<6; ++i)
            m[i] += cov[index][i];
        }

        ncov.push_back(m);
    }
}
/// @endcond

} // namespace internal

// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/**
   \ingroup PkgPointSetProcessingAlgorithms
   computes the Voronoi Covariance Measure (VCM) of a point cloud,
   a construction that can be used for normal estimation and sharp feature detection.

   The VCM associates to each point the covariance matrix of its Voronoi
   cell intersected with the ball of radius `offset_radius`.
   In addition, if the second radius `convolution_radius` is positive, the covariance matrices are smoothed
   via a convolution process. More specifically, each covariance matrix is replaced by
   the average of the matrices of the points located at a distance at most `convolution_radius`.
   The choice for parameter `offset_radius` should refer to the geometry of the underlying surface while
   the choice for parameter `convolution_radius` should refer to the noise level in the point cloud.
   For example, if the point cloud is a uniform and noise-free sampling of a smooth surface,
   `offset_radius` should be set to the minimum local feature size of the surface, while `convolution_radius` can be set to zero.

   The Voronoi covariance matrix of each vertex is stored in an array `a` of length 6 and is as follow:

   <center>
   \f$ \begin{bmatrix}
   a[0] & a[1] & a[2] \\
   a[1] & a[3] & a[4] \\
   a[2] & a[4] & a[5] \\
   \end{bmatrix}\f$
   </center>

   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param ccov output range of covariance matrices.
   \param offset_radius offset_radius.
   \param convolution_radius convolution_radius.
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`.
     If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd

   \sa `CGAL::vcm_is_on_feature_edge()`
   \sa `CGAL::vcm_estimate_normals()`

*/
template <typename PointRange,
          typename NamedParameters>
void
compute_vcm (const PointRange& points,
             std::vector< cpp11::array<double, 6> > &ccov,
             double offset_radius,
             double convolution_radius,
             const NamedParameters& np)
{
    using boost::choose_param;
    // basic geometric types
    typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
    typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;

    PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
    Kernel kernel;
    
    // First, compute the VCM for each point
    std::vector< cpp11::array<double, 6> > cov;
    std::size_t N = 20;
    internal::vcm_offset (points.begin(), points.end(),
                          point_map,
                          cov,
                          offset_radius,
                          N,
                          kernel);
    // Then, convolve it (only when convolution_radius != 0)
    if (convolution_radius == 0) {
        ccov.reserve(cov.size());
        std::copy(cov.begin(), cov.end(), std::back_inserter(ccov));
    } else {
      internal::vcm_convolve(points.begin(), points.end(),
                               point_map,
                               cov,
                               ccov,
                               convolution_radius,
                               kernel);
    }
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename PointRange>
void
compute_vcm (const PointRange& points,
             std::vector< cpp11::array<double, 6> > &ccov,
             double offset_radius,
             double convolution_radius)
{
  compute_vcm (points, ccov, offset_radius, convolution_radius,
               CGAL::Point_set_processing_3::parameters::all_default (points));
}

#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API
template < class ForwardIterator,
           class PointMap,
           class Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::compute_vcm(), please update your code")
void
compute_vcm (ForwardIterator first,
             ForwardIterator beyond,
             PointMap point_map,
             std::vector< cpp11::array<double, 6> > &ccov,
             double offset_radius,
             double convolution_radius,
             const Kernel & kernel)
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  compute_vcm (points, ccov, offset_radius, convolution_radius,
               CGAL::parameters::point_map (point_map).
               geom_traits (kernel));
}
#endif // CGAL_NO_DEPRECATED_CODE
/// \endcond
  
/// \cond SKIP_IN_MANUAL
template <typename PointRange,
          typename NamedParameters
>
void
vcm_estimate_normals_internal (PointRange& points,
                               double offset_radius, ///< offset radius.
                               double convolution_radius, ///< convolution radius.
                               const NamedParameters& np,
                               int nb_neighbors_convolve = -1 ///< number of neighbors used during the convolution.
)
{
    using boost::choose_param;
    // basic geometric types
    typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
    typedef typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::type NormalMap;
    typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;
    typedef typename GetDiagonalizeTraits<NamedParameters, double, 3>::type DiagonalizeTraits;

    CGAL_static_assertion_msg(!(boost::is_same<NormalMap,
                                typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::NoMap>::value),
                              "Error: no normal map");

    PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
    NormalMap normal_map = choose_param(get_param(np, internal_np::normal_map), NormalMap());
    
    typedef cpp11::array<double, 6> Covariance;
    
    // Compute the VCM and convolve it
    std::vector<Covariance> cov;
    if (nb_neighbors_convolve == -1) {
        compute_vcm(points,
                    cov,
                    offset_radius,
                    convolution_radius,
                    np);
    } else {
        internal::vcm_offset(points.begin(), points.end(),
                             point_map,
                             cov,
                             offset_radius,
                             20,
                             Kernel());

        if (nb_neighbors_convolve > 0)
        {
          std::vector<Covariance> ccov;
          ccov.reserve(cov.size());
          internal::vcm_convolve(points.begin(), points.end(),
                                 point_map,
                                 cov,
                                 ccov,
                                 (unsigned int) nb_neighbors_convolve,
                                 Kernel());

          cov.clear();
          std::copy(ccov.begin(), ccov.end(), std::back_inserter(cov));
        }
    }

    // And finally, compute the normals
    int i = 0;
    for (typename PointRange::iterator it = points.begin(); it != points.end(); ++it) {
        cpp11::array<double, 3> enormal = {{ 0,0,0 }};
        DiagonalizeTraits::extract_largest_eigenvector_of_covariance_matrix
          (cov[i], enormal);

        typename Kernel::Vector_3 normal(enormal[0],
                                         enormal[1],
                                         enormal[2]);
        put(normal_map, *it, normal);
        i++;
    }
}
/// @endcond


/**  
   \ingroup PkgPointSetProcessingAlgorithms
   Estimates normal directions of the range of `points`
   using the Voronoi Covariance Measure with a radius for the convolution.
   The output normals are randomly oriented.

   See `compute_vcm()` for a detailed description of the parameters `offset_radius` and `convolution_radius`
   and of the Voronoi Covariance Measure.

   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param offset_radius offset_radius.
   \param convolution_radius convolution_radius.
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`.
     If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{normal_map} a model of `WritablePropertyMap` with value type
     `geom_traits::Vector_3`.\cgalParamEnd
       \cgalParamBegin{diagonalize_traits} a model of `DiagonalizeTraits`. It can be omitted:
       if Eigen 3 (or greater) is available and `CGAL_EIGEN3_ENABLED` is defined then an overload
       using `Eigen_diagonalize_traits` is provided. Otherwise, the internal implementation
       `CGAL::Diagonalize_traits` is used.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd

*/
template <typename PointRange,
          typename NamedParameters
>
void
vcm_estimate_normals (PointRange& points,
                      double offset_radius,
                      double convolution_radius,
                      const NamedParameters& np
)
{
  vcm_estimate_normals_internal(points, offset_radius, convolution_radius, np);
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename PointRange>
void
vcm_estimate_normals (PointRange& points,
                      double offset_radius, ///< offset radius.
                      double convolution_radius) ///< convolution radius.
{
  return vcm_estimate_normals
    (points, offset_radius, convolution_radius,
     CGAL::Point_set_processing_3::parameters::all_default(points));
}

#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API
template < typename ForwardIterator,
           typename PointMap,
           typename NormalMap,
           typename VCMTraits
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::vcm_estimate_normals(), please update your code")
void
vcm_estimate_normals (ForwardIterator first, ///< iterator over the first input point.
                      ForwardIterator beyond, ///< past-the-end iterator over the input points.
                      PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3.
                      NormalMap normal_map, ///< property map: value_type of ForwardIterator -> Vector_3.
                      double offset_radius, ///< offset radius.
                      double convolution_radius, ///< convolution radius.
                      VCMTraits
)
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  vcm_estimate_normals
    (points,
     offset_radius, convolution_radius,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map).
     diagonalize_traits (VCMTraits()));
}
#endif // CGAL_NO_DEPRECATED_CODE
/// \endcond


/**
   \ingroup PkgPointSetProcessingAlgorithms
   Estimates normal directions of the range of `points`
   using the Voronoi Covariance Measure with a number of neighbors for the convolution.
   The output normals are randomly oriented.

   See `compute_vcm()` for a detailed description of the parameter `offset_radius`
   and of the Voronoi Covariance Measure.

   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param offset_radius offset_radius.
   \param k number of neighbor points used for convolution.
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`.
     If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{normal_map} a model of `WritablePropertyMap` with value type
     `geom_traits::Vector_3`.\cgalParamEnd
       \cgalParamBegin{diagonalize_traits} a model of `DiagonalizeTraits`. It can be omitted:
       if Eigen 3 (or greater) is available and `CGAL_EIGEN3_ENABLED` is defined then an overload
       using `Eigen_diagonalize_traits` is provided. Otherwise, the internal implementation
       `CGAL::Diagonalize_traits` is used.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd
*/
template < typename PointRange,
           typename NamedParameters
>
void
vcm_estimate_normals (PointRange& points,
                      double offset_radius,
                      unsigned int k,
                      const NamedParameters& np
)
{
  vcm_estimate_normals_internal(points, offset_radius, 0, np, k);
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename PointRange>
void
vcm_estimate_normals (PointRange& points,
                      double offset_radius, ///< offset radius.
                      unsigned int k)
{
  return vcm_estimate_normals
    (points, offset_radius, k,
     CGAL::Point_set_processing_3::parameters::all_default(points));
}

#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API
template < typename ForwardIterator,
           typename PointMap,
           typename NormalMap,
           typename VCMTraits
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::vcm_estimate_normals(), please update your code")
void
vcm_estimate_normals (ForwardIterator first, ///< iterator over the first input point.
                      ForwardIterator beyond, ///< past-the-end iterator over the input points.
                      PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3.
                      NormalMap normal_map, ///< property map: value_type of ForwardIterator -> Vector_3.
                      double offset_radius, ///< offset radius.
                      unsigned int k, ///< number of neighbor points used for the convolution.
                      VCMTraits
)
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  vcm_estimate_normals
    (points,
     offset_radius, k,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map).
     diagonalize_traits (VCMTraits()));
}

// deprecated API  
template < typename ForwardIterator,
           typename PointMap,
           typename NormalMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::vcm_estimate_normals(), please update your code")
void
vcm_estimate_normals (ForwardIterator first,
                      ForwardIterator beyond,
                      PointMap point_map,
                      NormalMap normal_map,
                      double offset_radius,
                      double convolution_radius)
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  vcm_estimate_normals
    (points,
     offset_radius, convolution_radius,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map));
}

// deprecated API  
template < typename ForwardIterator,
           typename PointMap,
           typename NormalMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::vcm_estimate_normals(), please update your code")
void
vcm_estimate_normals (ForwardIterator first,
                      ForwardIterator beyond,
                      PointMap point_map,
                      NormalMap normal_map,
                      double offset_radius,
                      unsigned int nb_neighbors_convolve)
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  vcm_estimate_normals
    (points,
     offset_radius, nb_neighbors_convolve,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map));
}


// deprecated API  
template < typename ForwardIterator,
           typename NormalMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::vcm_estimate_normals(), please update your code")
void
vcm_estimate_normals (ForwardIterator first,
                      ForwardIterator beyond,
                      NormalMap normal_map,
                      double offset_radius,
                      double convolution_radius) {
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  vcm_estimate_normals
    (points,
     offset_radius, convolution_radius,
     CGAL::parameters::normal_map (normal_map));
}

// deprecated API  
template < typename ForwardIterator,
           typename NormalMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::vcm_estimate_normals(), please update your code")
void
vcm_estimate_normals (ForwardIterator first,
                      ForwardIterator beyond,
                      NormalMap normal_map,
                      double offset_radius,
                      unsigned int nb_neighbors_convolve) {
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  vcm_estimate_normals
    (points,
     offset_radius, nb_neighbors_convolve,
     CGAL::parameters::normal_map (normal_map));
}
#endif // CGAL_NO_DEPRECATED_CODE
/// \endcond

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_VCM_ESTIMATE_NORMALS_H
