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
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H
#define CGAL_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H

#include <vector>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Default_diagonalize_traits.h>
#include <CGAL/centroid.h>

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

    \brief Class that precomputes and stored the eigenvectors and
    eigenvalues of the covariance matrices of all points of a point
    set using a local neighborhood.

    \tparam Geom_traits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type is
    `RandomAccessIterator`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `Geom_traits::Point_3`.
    \tparam DiagonalizeTraits model of `DiagonalizeTraits` used
    for matrix diagonalization.
  */
template <typename Geom_traits, typename PointRange, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Local_eigen_analysis
{
public:
  typedef typename Geom_traits::Point_3 Point;
  typedef typename Geom_traits::Vector_3 Vector;
  typedef typename Geom_traits::Plane_3 Plane;
  
  typedef CGAL::cpp11::array<double, 3> Eigenvalues; ///< Eigenvalues (sorted in ascending order)


private:

#ifdef CGAL_LINKED_WITH_TBB
  template <typename NeighborQuery>
  class Compute_eigen_values
  {
    Local_eigen_analysis& m_eigen;
    const PointRange& m_input;
    PointMap m_point_map;
    const NeighborQuery& m_neighbor_query;
    double& m_mean_range;
    tbb::mutex& m_mutex;
    
  public:
    
    Compute_eigen_values (Local_eigen_analysis& eigen,
                          const PointRange& input,
                          PointMap point_map,
                          const NeighborQuery& neighbor_query,
                          double& mean_range,
                          tbb::mutex& mutex)
      : m_eigen (eigen), m_input (input), m_point_map (point_map),
        m_neighbor_query (neighbor_query), m_mean_range (mean_range), m_mutex (mutex)
    { }
    
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t i = r.begin(); i != r.end(); ++ i)
        {
          std::vector<std::size_t> neighbors;
          m_neighbor_query (get(m_point_map, *(m_input.begin()+i)), std::back_inserter (neighbors));

          std::vector<Point> neighbor_points;
          for (std::size_t j = 0; j < neighbors.size(); ++ j)
            neighbor_points.push_back (get(m_point_map, *(m_input.begin()+neighbors[j])));

          m_mutex.lock();
          m_mean_range += CGAL::sqrt (CGAL::squared_distance (get(m_point_map, *(m_input.begin() + i)),
                                                              get(m_point_map, *(m_input.begin() + neighbors.back()))));
          m_mutex.unlock();
          
          m_eigen.compute (i, get(m_point_map, *(m_input.begin()+i)), neighbor_points);
        }
    }

  };
#endif
  
  std::vector<Eigenvalues> m_eigenvalues;
  std::vector<double> m_sum_eigenvalues;
  std::vector<Point> m_centroids;
  std::vector<Vector> m_smallest_eigenvectors;
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
  std::vector<Vector> m_middle_eigenvectors;
  std::vector<Vector> m_largest_eigenvectors;
#endif
  double m_mean_range;

  
public:

  /// \cond SKIP_IN_MANUAL
  Local_eigen_analysis () { }
  /// \endcond

  /*! 

    \brief Computes the local eigen analysis of an input range based
    on a local neighborhood.

    \tparam NeighborQuery model of `NeighborQuery`
    \tparam ConcurrencyTag enables sequential versus parallel
    algorithm. Possible values are `Parallel_tag` (default value is %CGAL
    is linked with TBB) or `Sequential_tag` (default value otherwise).
    \param input input range.
    \param point_map property map to access the input points
    \param neighbor_query object used to access neighborhoods of points.
  */
  template <typename NeighborQuery,
#if defined(DOXYGEN_RUNNING)
            typename ConcurrencyTag>
#elif defined(CGAL_LINKED_WITH_TBB)
            typename ConcurrencyTag = CGAL::Parallel_tag>
#else
            typename ConcurrencyTag = CGAL::Sequential_tag>
#endif
  Local_eigen_analysis (const PointRange& input,
                        PointMap point_map,
                        const NeighborQuery& neighbor_query)
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
        Compute_eigen_values<NeighborQuery> f(*this, input, point_map, neighbor_query, m_mean_range, mutex);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size ()), f);
      }
    else
#endif
      {
        for (std::size_t i = 0; i < input.size(); i++)
          {
            std::vector<std::size_t> neighbors;
            neighbor_query (get(point_map, *(input.begin()+i)), std::back_inserter (neighbors));

            std::vector<Point> neighbor_points;
            for (std::size_t j = 0; j < neighbors.size(); ++ j)
              neighbor_points.push_back (get(point_map, *(input.begin()+neighbors[j])));

            m_mean_range += CGAL::sqrt (CGAL::squared_distance (get(point_map, *(input.begin() + i)),
                                                                get(point_map, *(input.begin() + neighbors.back()))));
        
            compute (i, get(point_map, *(input.begin()+i)), neighbor_points);
          }
      }
    m_mean_range /= input.size();
  }

  /*!
    \brief Returns the estimated unoriented normal vector of the point at position `index`.
  */
  const Vector& normal_vector (std::size_t index) const { return m_smallest_eigenvectors[index]; }

  /*!
    \brief Returns the estimated local tangent plane of the point at position `index`.
  */
  Plane plane (std::size_t index) const { return Plane (m_centroids[index], m_smallest_eigenvectors[index]); }

  /*!
    \brief Returns the normalized eigenvalues of the point at position `index`.
  */
  const Eigenvalues& eigenvalue (std::size_t index) const { return m_eigenvalues[index]; }

  /*!
    \brief Returns the sum of eigenvalues of the point at position `index`.
  */
  const double& sum_of_eigenvalues (std::size_t index) const { return m_sum_eigenvalues[index]; }

  /// \cond SKIP_IN_MANUAL
  double mean_range() const { return m_mean_range; }
  /// \endcond

private:

  void compute (std::size_t index, const Point& query, std::vector<Point>& neighbor_points)
  {
    if (neighbor_points.size() == 0)
      {
        Eigenvalues v = {{ 0., 0., 0. }};
        m_eigenvalues[index] = v;
        m_centroids[index] = query;
        m_smallest_eigenvectors[index] = Vector (0., 0., 1.);
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
        m_middle_eigenvectors[index] = Vector (0., 1., 0.);
        m_largest_eigenvectors[index] = Vector (1., 0., 0.);
#endif
        return;
      }

    Point centroid = CGAL::centroid (neighbor_points.begin(), neighbor_points.end());
    m_centroids[index] = centroid;
    
    CGAL::cpp11::array<double, 6> covariance = {{ 0., 0., 0., 0., 0., 0. }};
      
    for (std::size_t i = 0; i < neighbor_points.size(); ++ i)
      {
        Vector d = neighbor_points[i] - centroid;
        covariance[0] += d.x () * d.x ();
        covariance[1] += d.x () * d.y ();
        covariance[2] += d.x () * d.z ();
        covariance[3] += d.y () * d.y ();
        covariance[4] += d.y () * d.z ();
        covariance[5] += d.z () * d.z ();
      }

    Eigenvalues evalues = {{ 0., 0., 0. }};
    CGAL::cpp11::array<double, 9> evectors = {{ 0., 0., 0.,
                                                0., 0., 0.,
                                                0., 0., 0. }};

    DiagonalizeTraits::diagonalize_selfadjoint_covariance_matrix
      (covariance, evalues, evectors);

    // Normalize
    double sum = evalues[0] + evalues[1] + evalues[2];
    if (sum > 0.)
      for (std::size_t i = 0; i < 3; ++ i)
        evalues[i] = evalues[i] / sum;
    m_sum_eigenvalues[index] = sum;
    m_eigenvalues[index] = evalues;
    m_smallest_eigenvectors[index] = Vector (evectors[0], evectors[1], evectors[2]);
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
    m_middle_eigenvectors[index] = Vector (evectors[3], evectors[4], evectors[5]);
    m_largest_eigenvectors[index] = Vector (evectors[6], evectors[7], evectors[8]);
#endif
  }

};
  

}
  
}


#endif // CGAL_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H
