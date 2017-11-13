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
//
// Author(s) : Jocelyn Meyron and Quentin Mérigot
//

#ifndef CGAL_VCM_ESTIMATE_NORMALS_H
#define CGAL_VCM_ESTIMATE_NORMALS_H

#include <CGAL/license/Point_set_processing_3.h>


#include <CGAL/internal/Voronoi_covariance_3/voronoi_covariance_3.h>

#include <CGAL/property_map.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Fuzzy_sphere.h>

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
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = `Kernel::Point_3`.
/// @tparam K Geometric traits class.
/// @tparam Covariance Covariance matrix type. It is similar to an array with a length of 6.
template < typename ForwardIterator,
           typename PointPMap,
           class K,
           class Covariance
>
void
vcm_offset (ForwardIterator first, ///< iterator over the first input point.
            ForwardIterator beyond, ///< past-the-end iterator over the input points.
            PointPMap point_pmap, ///< property map: value_type of ForwardIterator -> Point_3.
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
      points.push_back(get(point_pmap, *it));

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
           class PointPMap,
           class K,
           class Covariance
>
void
vcm_convolve (ForwardIterator first,
              ForwardIterator beyond,
              PointPMap point_pmap,
              const std::vector<Covariance> &cov,
              std::vector<Covariance> &ncov,
              double convolution_radius,
              const K &)
{
    typedef std::pair<typename K::Point_3, std::size_t>              Tree_point;
    typedef First_of_pair_property_map< Tree_point >                  Tree_pmap;
    typedef Search_traits_3<K>                                      Traits_base;
    typedef Search_traits_adapter<Tree_point, Tree_pmap, Traits_base>    Traits;
    typedef Kd_tree<Traits>                                                Tree;
    typedef Fuzzy_sphere<Traits>                                   Fuzzy_sphere;

    // Kd tree
    Tree tree;
    tree.reserve(cov.size());
    std::size_t i=0;
    for (ForwardIterator it = first; it != beyond; ++it, ++i)
        tree.insert( Tree_point(get(point_pmap, *it), i) );

    // Convolving
    ncov.clear();
    ncov.reserve(cov.size());
    for (ForwardIterator it = first; it != beyond; ++it) {
        std::vector<Tree_point> nn;
        tree.search(std::back_inserter(nn),
                    Fuzzy_sphere (get(point_pmap, *it), convolution_radius));

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
           class PointPMap,
           class K,
           class Covariance
>
void
vcm_convolve (ForwardIterator first,
              ForwardIterator beyond,
              PointPMap point_pmap,
              const std::vector<Covariance> &cov,
              std::vector<Covariance> &ncov,
              unsigned int nb_neighbors_convolve,
              const K &)
{
    typedef std::pair<typename K::Point_3, std::size_t>              Tree_point;
    typedef First_of_pair_property_map< Tree_point >                  Tree_pmap;
    typedef Search_traits_3<K>                                      Traits_base;
    typedef Search_traits_adapter<Tree_point, Tree_pmap, Traits_base>    Traits;
    typedef Orthogonal_k_neighbor_search<Traits>                Neighbor_search;
    typedef typename Neighbor_search::Tree                                 Tree;

    // Search tree
    Tree tree;
    tree.reserve(cov.size());
    std::size_t i=0;
    for (ForwardIterator it = first; it != beyond; ++it, ++i)
        tree.insert( Tree_point(get(point_pmap, *it), i) );

    // Convolving
    ncov.clear();
    ncov.reserve(cov.size());
    for (ForwardIterator it = first; it != beyond; ++it) {
        Neighbor_search search(tree, get(point_pmap, *it), nb_neighbors_convolve);

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

/// \ingroup PkgPointSetProcessingAlgorithms
/// computes the Voronoi Covariance Measure (VCM) of a point cloud,
/// a construction that can be used for normal estimation and sharp feature detection.
///
/// The VCM associates to each point the covariance matrix of its Voronoi
/// cell intersected with the ball of radius `offset_radius`.
/// In addition, if the second radius `convolution_radius` is positive, the covariance matrices are smoothed
/// via a convolution process. More specifically, each covariance matrix is replaced by
/// the average of the matrices of the points located at a distance at most `convolution_radius`.
/// The choice for parameter `offset_radius` should refer to the geometry of the underlying surface while
/// the choice for parameter `convolution_radius` should refer to the noise level in the point cloud.
/// For example, if the point cloud is a uniform and noise-free sampling of a smooth surface,
/// `offset_radius` should be set to the minimum local feature size of the surface, while `convolution_radius` can be set to zero.
///
/// The Voronoi covariance matrix of each vertex is stored in an array `a` of length 6 and is as follow:
/*!
<center>
\f$ \begin{bmatrix}
  a[0] & a[1] & a[2] \\
  a[1] & a[3] & a[4] \\
  a[2] & a[4] & a[5] \\
 \end{bmatrix}\f$
</center>
*/
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = `Kernel::Point_3`.
/// @tparam CovariancePMap is a model of `ReadWritePropertyMap` with a value_type = `cpp11::array<Kernel::FT,6>`.
/// @tparam Kernel Geometric traits class.
///
/// \sa `CGAL::vcm_is_on_feature_edge()`
/// \sa `CGAL::vcm_estimate_normals()`
///
/// \todo replace ccov by a property map.
template < class ForwardIterator,
           class PointPMap,
           class Kernel
>
void
compute_vcm (ForwardIterator first,
             ForwardIterator beyond,
             PointPMap point_pmap,
             std::vector< cpp11::array<typename Kernel::FT, 6> > &ccov,
             double offset_radius,
             double convolution_radius,
             const Kernel & kernel)
{
    // First, compute the VCM for each point
    std::vector< cpp11::array<typename Kernel::FT, 6> > cov;
    std::size_t N = 20;
    internal::vcm_offset (first, beyond,
                          point_pmap,
                          cov,
                          offset_radius,
                          N,
                          kernel);
    // Then, convolve it (only when convolution_radius != 0)
    if (convolution_radius == 0) {
        ccov.reserve(cov.size());
        std::copy(cov.begin(), cov.end(), std::back_inserter(ccov));
    } else {
        internal::vcm_convolve(first, beyond,
                               point_pmap,
                               cov,
                               ccov,
                               convolution_radius,
                               kernel);
    }
}

/// @cond SKIP_IN_MANUAL
/// Estimates normal directions of the `[first, beyond)` range of points
/// using the Voronoi Covariance Measure (see `compute_vcm` for more details on the VCM).
/// The output normals are randomly oriented.
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = `Kernel::Point_3`.
/// @tparam NormalPMap is a model of `WritablePropertyMap` with a value_type = `Kernel::Vector_3`.
/// @tparam Kernel Geometric traits class.
/// @tparam Covariance Covariance matrix type. It is similar to an array with a length of 6.
/// @pre If `nb_neighbors_convolve` is equal to -1, then the convolution is made using a radius.
/// On the contrary, if `nb_neighbors_convolve` is different from -1, the convolution is made using
/// this number of neighbors.

// This variant requires all of the parameters.
template < typename VCMTraits,
           typename ForwardIterator,
           typename PointPMap,
           typename NormalPMap,
           typename Kernel
>
void
vcm_estimate_normals (ForwardIterator first, ///< iterator over the first input point.
                      ForwardIterator beyond, ///< past-the-end iterator over the input points.
                      PointPMap point_pmap, ///< property map: value_type of ForwardIterator -> Point_3.
                      NormalPMap normal_pmap, ///< property map: value_type of ForwardIterator -> Vector_3.
                      double offset_radius, ///< offset radius.
                      double convolution_radius, ///< convolution radius.
                      const Kernel & kernel, ///< geometric traits.
                      int nb_neighbors_convolve = -1 ///< number of neighbors used during the convolution.
)
{
    typedef cpp11::array<double, 6> Covariance;
    // Compute the VCM and convolve it
    std::vector<Covariance> cov;
    if (nb_neighbors_convolve == -1) {
        compute_vcm(first, beyond,
                    point_pmap,
                    cov,
                    offset_radius,
                    convolution_radius,
                    kernel);
    } else {
        internal::vcm_offset(first, beyond,
                             point_pmap,
                             cov,
                             offset_radius,
                             20,
                             kernel);

        if (nb_neighbors_convolve > 0)
        {
          std::vector<Covariance> ccov;
          ccov.reserve(cov.size());
          internal::vcm_convolve(first, beyond,
                                 point_pmap,
                                 cov,
                                 ccov,
                                 (unsigned int) nb_neighbors_convolve,
                                 kernel);

          cov.clear();
          std::copy(ccov.begin(), ccov.end(), std::back_inserter(cov));
        }
    }

    // And finally, compute the normals
    int i = 0;
    for (ForwardIterator it = first; it != beyond; ++it) {
        cpp11::array<double, 3> enormal = {{ 0,0,0 }};
        VCMTraits::extract_largest_eigenvector_of_covariance_matrix
          (cov[i], enormal);

        typename Kernel::Vector_3 normal(enormal[0],
                                         enormal[1],
                                         enormal[2]);
        put(normal_pmap, *it, normal);
        i++;
    }
}
/// @endcond

/// \ingroup PkgPointSetProcessingAlgorithms
/// Estimates normal directions of the points in the range `[first, beyond)`
/// using the Voronoi Covariance Measure with a radius for the convolution.
/// The output normals are randomly oriented.
///
/// See `compute_vcm()` for a detailed description of the parameters `offset_radius` and `convolution_radius`
/// and of the Voronoi Covariance Measure.
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = `Kernel::Point_3`.
/// @tparam NormalPMap is a model of `WritablePropertyMap` with a value_type = `Kernel::Vector_3`.
/// \tparam VCMTraits is a model of `DiagonalizeTraits`. It can be
/// omitted: if Eigen 3 (or greater) is available and
/// `CGAL_EIGEN3_ENABLED` is defined then an overload using
/// `Eigen_diagonalize_traits` is provided. Otherwise, the internal
/// implementation `Diagonalize_traits` is used.
// This variant deduces the kernel from the point property map
// and uses a radius for the convolution.
template < typename ForwardIterator,
           typename PointPMap,
           typename NormalPMap,
           typename VCMTraits
>
void
vcm_estimate_normals (ForwardIterator first, ///< iterator over the first input point.
                      ForwardIterator beyond, ///< past-the-end iterator over the input points.
                      PointPMap point_pmap, ///< property map: value_type of ForwardIterator -> Point_3.
                      NormalPMap normal_pmap, ///< property map: value_type of ForwardIterator -> Vector_3.
                      double offset_radius, ///< offset radius.
                      double convolution_radius, ///< convolution radius.
                      VCMTraits
)
{
    typedef typename boost::property_traits<PointPMap>::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;

    vcm_estimate_normals<VCMTraits>(first, beyond,
				    point_pmap, normal_pmap,
				    offset_radius, convolution_radius,
				    Kernel());
}


/// \ingroup PkgPointSetProcessingAlgorithms
/// Estimates normal directions of the points in the range `[first, beyond)`
/// using the Voronoi Covariance Measure with a number of neighbors for the convolution.
/// The output normals are randomly oriented.
///
/// See `compute_vcm()` for a detailed description of the parameter `offset_radius`
/// and of the Voronoi Covariance Measure.
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = `Kernel::Point_3`.
/// @tparam NormalPMap is a model of `WritablePropertyMap` with a value_type = `Kernel::Vector_3`.
/// \tparam VCMTraits is a model of `DiagonalizeTraits`. It can be
/// omitted: if Eigen 3 (or greater) is available and
/// `CGAL_EIGEN3_ENABLED` is defined then an overload using
/// `Eigen_diagonalize_traits` is provided. Otherwise, the internal
/// implementation `Diagonalize_traits` is used.

// This variant deduces the kernel from the point property map
// and uses a number of neighbors for the convolution.
template < typename ForwardIterator,
           typename PointPMap,
           typename NormalPMap,
           typename VCMTraits
>
void
vcm_estimate_normals (ForwardIterator first, ///< iterator over the first input point.
                      ForwardIterator beyond, ///< past-the-end iterator over the input points.
                      PointPMap point_pmap, ///< property map: value_type of ForwardIterator -> Point_3.
                      NormalPMap normal_pmap, ///< property map: value_type of ForwardIterator -> Vector_3.
                      double offset_radius, ///< offset radius.
                      unsigned int k, ///< number of neighbor points used for the convolution.
                      VCMTraits
)
{
    typedef typename boost::property_traits<PointPMap>::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;

    vcm_estimate_normals<VCMTraits>(first, beyond,
				    point_pmap, normal_pmap,
				    offset_radius, 0,
				    Kernel(),
				    k);
}


template < typename ForwardIterator,
           typename PointPMap,
           typename NormalPMap
>
void
vcm_estimate_normals (ForwardIterator first,
                      ForwardIterator beyond,
                      PointPMap point_pmap,
                      NormalPMap normal_pmap,
                      double offset_radius,
                      double convolution_radius)
{
  vcm_estimate_normals(first, beyond, point_pmap, normal_pmap, offset_radius, convolution_radius,
		       CGAL::Default_diagonalize_traits<double, 3>());
}

template < typename ForwardIterator,
           typename PointPMap,
           typename NormalPMap
>
void
vcm_estimate_normals (ForwardIterator first,
                      ForwardIterator beyond,
                      PointPMap point_pmap,
                      NormalPMap normal_pmap,
                      double offset_radius,
                      unsigned int nb_neighbors_convolve)
{
  vcm_estimate_normals(first, beyond, point_pmap, normal_pmap, offset_radius, nb_neighbors_convolve,
		       CGAL::Default_diagonalize_traits<double, 3>());

}


/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map
// and use a radius for the convolution.
template < typename ForwardIterator,
           typename NormalPMap
>
void
vcm_estimate_normals (ForwardIterator first,
                      ForwardIterator beyond,
                      NormalPMap normal_pmap,
                      double offset_radius,
                      double convolution_radius) {
    vcm_estimate_normals(first, beyond,
                         make_identity_property_map(typename std::iterator_traits<ForwardIterator>::value_type()),
                         normal_pmap,
                         offset_radius, convolution_radius);
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map
// and use a number of neighbors for the convolution.
template < typename ForwardIterator,
           typename NormalPMap
>
void
vcm_estimate_normals (ForwardIterator first,
                      ForwardIterator beyond,
                      NormalPMap normal_pmap,
                      double offset_radius,
                      unsigned int nb_neighbors_convolve) {
    vcm_estimate_normals(first, beyond,
                         make_identity_property_map(typename std::iterator_traits<ForwardIterator>::value_type()),
                         normal_pmap,
                         offset_radius, nb_neighbors_convolve);
}
/// @endcond

} // namespace CGAL

#endif // CGAL_VCM_ESTIMATE_NORMALS_H
