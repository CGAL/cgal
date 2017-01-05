#ifndef CGAL_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H
#define CGAL_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H

#include <vector>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Default_diagonalize_traits.h>
#include <CGAL/centroid.h>

namespace CGAL {

namespace Classification {

  /*!
    \ingroup PkgClassificationDataStructures

    \brief Class that precomputes and stored the eigenvectors and
    eigenvalues of the covariance matrices of all points of a point
    set using a local neighborhood.

    \tparam Kernel model of \cgal Kernel.
    \tparam Range range of items, model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `Range` and value type
    is `Point_3<Kernel>`.
    \tparam DiagonalizeTraits model of `DiagonalizeTraits` used
    for matrix diagonalization.
  */
template <typename Kernel, typename Range, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Local_eigen_analysis
{
public:
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Plane_3 Plane;
  
  typedef CGAL::cpp11::array<double, 3> Eigenvalues; ///< Eigenvalues (sorted in ascending order)


private:
  std::vector<Eigenvalues> m_eigenvalues;
  std::vector<double> m_sum_eigenvalues;
  std::vector<Point> m_centroids;
  std::vector<Vector> m_smallest_eigenvectors;
  std::vector<Vector> m_middle_eigenvectors;
  std::vector<Vector> m_largest_eigenvectors;
  double m_mean_range;

  
public:

  /// \cond SKIP_IN_MANUAL
  Local_eigen_analysis () { }
  /// \endcond

  /*! 

    \brief Computes the local eigen analysis of an input range based
    on a local neighborhood.

    \tparam NeighborQuery model of `NeighborQuery`
    \param input Input range.
    \param point_map property map to access the input points
    \param neighbor_query object used to access neighborhoods of points
  */
  template <typename NeighborQuery>
  Local_eigen_analysis (const Range& input,
                        PointMap point_map,
                        const NeighborQuery& neighbor_query)
  {
    m_eigenvalues.reserve (input.size());
    m_centroids.reserve (input.size());
    m_smallest_eigenvectors.reserve (input.size());
    m_middle_eigenvectors.reserve (input.size());
    m_largest_eigenvectors.reserve (input.size());
    
    m_mean_range = 0.;
      
    for (std::size_t i = 0; i < input.size(); i++)
      {
        std::vector<std::size_t> neighbors;
        neighbor_query (get(point_map, input[i]), std::back_inserter (neighbors));

        std::vector<Point> neighbor_points;
        for (std::size_t j = 0; j < neighbors.size(); ++ j)
          neighbor_points.push_back (get(point_map, input[neighbors[j]]));

        m_mean_range += CGAL::sqrt (CGAL::squared_distance (get(point_map, input[i]),
                                                          get(point_map, input[neighbors.back()])));
        
        compute (get(point_map, input[i]), neighbor_points);
      }

    m_mean_range /= input.size();
  }

  /*!
    \brief Returns the estimated unoriented normal vector of the indexed point.
  */
  const Vector& normal_vector (std::size_t index) const { return m_smallest_eigenvectors[index]; }

  /*!
    \brief Returns the estimated local tangent plane of the index point.
  */
  Plane plane (std::size_t index) const { return Plane (m_centroids[index], m_smallest_eigenvectors[index]); }

  /*!
    \brief Returns the normalized eigenvalues of the index point.
  */
  const Eigenvalues& eigenvalue (std::size_t index) const { return m_eigenvalues[index]; }

  /*!
    \brief Returns the sum of eigenvalues of the index point.
  */
  const double& sum_of_eigenvalues (std::size_t index) const { return m_sum_eigenvalues[index]; }

  /// \cond SKIP_IN_MANUAL
  double mean_range() const { return m_mean_range; }
  /// \endcond

private:

  void compute (const Point& query, std::vector<Point>& neighbor_points)
  {
    if (neighbor_points.size() == 0)
      {
        Eigenvalues v = {{ 0., 0., 0. }};
        m_eigenvalues.push_back (v);
        m_centroids.push_back (query);
        m_smallest_eigenvectors.push_back (Vector (0., 0., 1.));
        m_middle_eigenvectors.push_back (Vector (0., 1., 0.));
        m_largest_eigenvectors.push_back (Vector (1., 0., 0.));
        return;
      }

    Point centroid = CGAL::centroid (neighbor_points.begin(), neighbor_points.end());
    m_centroids.push_back (centroid);
    
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
    m_sum_eigenvalues.push_back(sum);    
    m_eigenvalues.push_back (evalues);
    m_smallest_eigenvectors.push_back (Vector (evectors[0], evectors[1], evectors[2]));
    m_middle_eigenvectors.push_back (Vector (evectors[3], evectors[4], evectors[5]));
    m_largest_eigenvectors.push_back (Vector (evectors[6], evectors[7], evectors[8]));

  }

};
  

}
  
}


#endif // CGAL_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H
