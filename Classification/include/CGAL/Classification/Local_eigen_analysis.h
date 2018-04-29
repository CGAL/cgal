// Copyright (c) 2017 GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H
#define CGAL_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H

#include <CGAL/license/Classification.h>

#include <vector>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Default_diagonalize_traits.h>
#include <CGAL/centroid.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/array.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>
#include <tbb/mutex.h>
#endif // CGAL_LINKED_WITH_TBB

namespace CGAL {

namespace Classification {

  /*!
    \ingroup PkgClassificationDataStructures

    \brief Class that precomputes and stores the eigenvectors and
    eigenvalues of the covariance matrices of all points of a point
    set using a local neighborhood.

    This class can be used to compute eigen features (see
    \ref PkgClassificationFeatures) and to estimate local normal vectors
    and tangent planes.

  */
class Local_eigen_analysis
{
public:
  typedef CGAL::cpp11::array<float, 3> Eigenvalues; ///< Eigenvalues (sorted in ascending order)
  
private:

#ifdef CGAL_LINKED_WITH_TBB
  template <typename PointRange, typename PointMap, typename NeighborQuery, typename DiagonalizeTraits>
  class Compute_eigen_values
  {
    Local_eigen_analysis& m_eigen;
    const PointRange& m_input;
    PointMap m_point_map;
    const NeighborQuery& m_neighbor_query;
    float& m_mean_range;
    tbb::mutex& m_mutex;
    
  public:
    
    Compute_eigen_values (Local_eigen_analysis& eigen,
                          const PointRange& input,
                          PointMap point_map,
                          const NeighborQuery& neighbor_query,
                          float& mean_range,
                          tbb::mutex& mutex)
      : m_eigen (eigen), m_input (input), m_point_map (point_map),
        m_neighbor_query (neighbor_query), m_mean_range (mean_range), m_mutex (mutex)
    { }
    
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      std::vector<std::size_t> neighbors;
      for (std::size_t i = r.begin(); i != r.end(); ++ i)
      {
        neighbors.clear();
        m_neighbor_query (get(m_point_map, *(m_input.begin()+i)), std::back_inserter (neighbors));

        std::vector<typename PointMap::value_type> neighbor_points;
        neighbor_points.reserve(neighbors.size());
        for (std::size_t j = 0; j < neighbors.size(); ++ j)
          neighbor_points.push_back (get(m_point_map, *(m_input.begin()+neighbors[j])));

        m_mutex.lock();
        m_mean_range += float(CGAL::sqrt
                              (CGAL::squared_distance (get(m_point_map, *(m_input.begin() + i)),
                                                       get(m_point_map, *(m_input.begin() + neighbors.back())))));
        m_mutex.unlock();
          
        m_eigen.compute<typename PointMap::value_type,
                        DiagonalizeTraits> (i, get(m_point_map, *(m_input.begin()+i)), neighbor_points);
      }
    }

  };
#endif

  template <typename ClusterRange, typename PointMap, typename DiagonalizeTraits>
  class Compute_clusters_eigen_values
  {
    Local_eigen_analysis& m_eigen;
    const ClusterRange& m_input;
    PointMap m_point_map;
    
  public:
    
    Compute_clusters_eigen_values (Local_eigen_analysis& eigen,
                                   const ClusterRange& input,
                                   PointMap point_map)
      : m_eigen (eigen), m_input (input), m_point_map (point_map)
    { }

#ifdef CGAL_LINKED_WITH_TBB
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t i = r.begin(); i != r.end(); ++ i)
        apply (i);
    }
#endif

    inline void apply (std::size_t i) const
    {
      const typename ClusterRange::value_type& cluster = m_input[i];

      std::vector<typename PointMap::value_type> points;
      for (std::size_t j = 0; j < cluster.size(); ++ j)
        points.push_back (get(m_point_map, cluster[j]));

      m_eigen.compute<typename PointMap::value_type,
                      DiagonalizeTraits> (i, typename PointMap::value_type(0.,0.,0.), points);
    }

  };

  
  typedef CGAL::cpp11::array<float, 3> float3;
  std::vector<float3> m_eigenvalues;
  std::vector<float> m_sum_eigenvalues;
  std::vector<float3> m_centroids;
  std::vector<float3> m_smallest_eigenvectors;
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
  std::vector<float3> m_middle_eigenvectors;
  std::vector<float3> m_largest_eigenvectors;
#endif
  float m_mean_range;

  
public:

  /// \cond SKIP_IN_MANUAL
  Local_eigen_analysis () { }
  /// \endcond

  /*! 
    \brief Computes the local eigen analysis of an input range based
    on a local neighborhood.

    \tparam PointRange model of `ConstRange`. Its iterator type is
    `RandomAccessIterator` and its value type is the key type of
    `PointMap`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `CGAL::Point_3`.
    \tparam NeighborQuery model of `NeighborQuery`
    \tparam ConcurrencyTag enables sequential versus parallel
    algorithm. Possible values are `Parallel_tag` (default value is %CGAL
    is linked with TBB) or `Sequential_tag` (default value otherwise).
    \tparam DiagonalizeTraits model of `DiagonalizeTraits` used for
    matrix diagonalization. It can be omitted: if Eigen 3 (or greater)
    is available and `CGAL_EIGEN3_ENABLED` is defined then an overload
    using `Eigen_diagonalize_traits` is provided. Otherwise, the
    internal implementation `Diagonalize_traits` is used.

    \param input point range.
    \param point_map property map to access the input points.
    \param neighbor_query object used to access neighborhoods of points.
  */
  template <typename PointRange,
            typename PointMap,
            typename NeighborQuery,
#if defined(DOXYGEN_RUNNING)
            typename ConcurrencyTag,
#elif defined(CGAL_LINKED_WITH_TBB)
            typename ConcurrencyTag = CGAL::Parallel_tag,
#else
            typename ConcurrencyTag = CGAL::Sequential_tag,
#endif
#if defined(DOXYGEN_RUNNING)
            typename DiagonalizeTraits>
#else
            typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<float, 3> >
#endif
    Local_eigen_analysis (const PointRange& input,
                          PointMap point_map,
                          const NeighborQuery& neighbor_query,
                          const ConcurrencyTag& = ConcurrencyTag(),
                          const DiagonalizeTraits& = DiagonalizeTraits())
  {
    m_eigenvalues.resize (input.size());
    m_sum_eigenvalues.resize (input.size());
    m_centroids.resize (input.size());
    m_smallest_eigenvectors.resize (input.size());
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
    m_middle_eigenvectors.resize (input.size());
    m_largest_eigenvectors.resize (input.size());
#endif
    
    m_mean_range = 0.;
      
#ifndef CGAL_LINKED_WITH_TBB
    CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                               "Parallel_tag is enabled but TBB is unavailable.");
#else
    if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
    {
      tbb::mutex mutex;
      Compute_eigen_values<PointRange, PointMap, NeighborQuery, DiagonalizeTraits>
        f(*this, input, point_map, neighbor_query, m_mean_range, mutex);
      tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size ()), f);
    }
    else
#endif
    {
      for (std::size_t i = 0; i < input.size(); i++)
      {
        std::vector<std::size_t> neighbors;
        neighbor_query (get(point_map, *(input.begin()+i)), std::back_inserter (neighbors));

        std::vector<typename PointMap::value_type> neighbor_points;
        for (std::size_t j = 0; j < neighbors.size(); ++ j)
          neighbor_points.push_back (get(point_map, *(input.begin()+neighbors[j])));

        m_mean_range += float(CGAL::sqrt (CGAL::squared_distance
                                          (get(point_map, *(input.begin() + i)),
                                           get(point_map, *(input.begin() + neighbors.back())))));
        
        compute<typename PointMap::value_type, DiagonalizeTraits>
          (i, get(point_map, *(input.begin()+i)), neighbor_points);
      }
    }
    m_mean_range /= input.size();
  }

  // Experimental feature, not used officially
  /// \cond SKIP_IN_MANUAL
  struct Input_is_clusters { };
  
  template <typename ClusterRange,
            typename PointMap,
#if defined(DOXYGEN_RUNNING)
            typename ConcurrencyTag,
#elif defined(CGAL_LINKED_WITH_TBB)
            typename ConcurrencyTag = CGAL::Parallel_tag,
#else
            typename ConcurrencyTag = CGAL::Sequential_tag,
#endif
#if defined(DOXYGEN_RUNNING)
            typename DiagonalizeTraits>
#else
            typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<float, 3> >
#endif
  Local_eigen_analysis (const Input_is_clusters&,
                        const ClusterRange& input,
                        PointMap point_map,
                        const ConcurrencyTag& = ConcurrencyTag(),
                        const DiagonalizeTraits& = DiagonalizeTraits())
  {
    m_eigenvalues.resize (input.size());
    m_sum_eigenvalues.resize (input.size());
    m_centroids.resize (input.size());
    m_smallest_eigenvectors.resize (input.size());
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
    m_middle_eigenvectors.resize (input.size());
    m_largest_eigenvectors.resize (input.size());
#endif
    
    m_mean_range = 0.;

    Compute_clusters_eigen_values<ClusterRange, PointMap, DiagonalizeTraits>
    f(*this, input, point_map);

    
#ifndef CGAL_LINKED_WITH_TBB
    CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                               "Parallel_tag is enabled but TBB is unavailable.");
#else
    if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size ()), f);
    }
    else
#endif
    {
      for (std::size_t i = 0; i < input.size(); ++ i)
        f.apply (i);
    }
  }
  /// \endcond

  /*!
    \brief Returns the estimated unoriented normal vector of the point at position `index`.
    \tparam GeomTraits model of \cgal Kernel.
  */
  template <typename GeomTraits>
  typename GeomTraits::Vector_3 normal_vector (std::size_t index) const
  {
    return typename GeomTraits::Vector_3(double(m_smallest_eigenvectors[index][0]),
                                          double(m_smallest_eigenvectors[index][1]),
                                          double(m_smallest_eigenvectors[index][2]));
  }

  /*!
    \brief Returns the estimated local tangent plane of the point at position `index`.
    \tparam GeomTraits model of \cgal Kernel.
  */
  template <typename GeomTraits>
  typename GeomTraits::Plane_3 plane (std::size_t index) const
  {
    return typename GeomTraits::Plane_3
      (typename GeomTraits::Point_3 (double(m_centroids[index][0]),
                                      double(m_centroids[index][1]),
                                      double(m_centroids[index][2])),
       typename GeomTraits::Vector_3 (double(m_smallest_eigenvectors[index][0]),
                                       double(m_smallest_eigenvectors[index][1]),
                                       double(m_smallest_eigenvectors[index][2])));
  }

  /*!
    \brief Returns the normalized eigenvalues of the point at position `index`.
  */
  const Eigenvalues& eigenvalue (std::size_t index) const { return m_eigenvalues[index]; }

  /*!
    \brief Returns the sum of eigenvalues of the point at position `index`.
  */
  float sum_of_eigenvalues (std::size_t index) const { return m_sum_eigenvalues[index]; }

  /// \cond SKIP_IN_MANUAL
  float mean_range() const { return m_mean_range; }
  /// \endcond

private:

  template <typename Point, typename DiagonalizeTraits>
  void compute (std::size_t index, const Point& query, std::vector<Point>& neighbor_points)
  {
    typedef typename Kernel_traits<Point>::Kernel::Vector_3 Vector;
      
    if (neighbor_points.size() == 0)
    {
      Eigenvalues v = make_array( 0.f, 0.f, 0.f );
      m_eigenvalues[index] = v;
      m_centroids[index] = make_array(float(query.x()), float(query.y()), float(query.z()) );
      m_smallest_eigenvectors[index] = make_array( 0.f, 0.f, 1.f );
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
      m_middle_eigenvectors[index] = make_array( 0.f, 1.f, 0.f );
      m_largest_eigenvectors[index] = make_array( 1.f, 0.f, 0.f );
#endif
      return;
    }

    Point centroid = CGAL::centroid (neighbor_points.begin(), neighbor_points.end());
    m_centroids[index] = make_array( float(centroid.x()), float(centroid.y()), float(centroid.z()) );
    
    CGAL::cpp11::array<float, 6> covariance = make_array( 0.f, 0.f, 0.f, 0.f, 0.f, 0.f );
      
    for (std::size_t i = 0; i < neighbor_points.size(); ++ i)
    {
      Vector d = neighbor_points[i] - centroid;
      covariance[0] += float(d.x () * d.x ());
      covariance[1] += float(d.x () * d.y ());
      covariance[2] += float(d.x () * d.z ());
      covariance[3] += float(d.y () * d.y ());
      covariance[4] += float(d.y () * d.z ());
      covariance[5] += float(d.z () * d.z ());
    }

    CGAL::cpp11::array<float, 3> evalues = make_array( 0.f, 0.f, 0.f );
    CGAL::cpp11::array<float, 9> evectors = make_array( 0.f, 0.f, 0.f,
                                               0.f, 0.f, 0.f,
                                               0.f, 0.f, 0.f );

    DiagonalizeTraits::diagonalize_selfadjoint_covariance_matrix
      (covariance, evalues, evectors);

    // Normalize
    float sum = evalues[0] + evalues[1] + evalues[2];
    if (sum > 0.f)
      for (std::size_t i = 0; i < 3; ++ i)
        evalues[i] = evalues[i] / sum;
    m_sum_eigenvalues[index] = float(sum);
    m_eigenvalues[index] = make_array( float(evalues[0]), float(evalues[1]), float(evalues[2]) );
    m_smallest_eigenvectors[index] = make_array( float(evectors[0]), float(evectors[1]), float(evectors[2]) );
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
    m_middle_eigenvectors[index] = make_array( float(evectors[3]), float(evectors[4]), float(evectors[5]) );
    m_largest_eigenvectors[index] = make_array( float(evectors[6]), float(evectors[7]), float(evectors[8]) );
#endif
  }

};
  

}
  
}


#endif // CGAL_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H
