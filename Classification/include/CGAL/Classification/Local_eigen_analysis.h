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
    \ingroup PkgClassification

    \brief Class that precomputes and stored the eigenvectors and
    eigenvalues of the covariance matrices of all points of a point
    set, using a local neighborhood.

    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
    \tparam DiagonalizeTraits Solver used for matrix diagonalization.
  */
template <typename Kernel, typename RandomAccessIterator, typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Local_eigen_analysis
{
public:
  /// \cond SKIP_IN_MANUAL
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Plane_3 Plane;
  /// \endcond
  
  typedef CGAL::cpp11::array<double, 3> Eigenvalues; ///< Eigenvalues (sorted in ascending order)


private:
  std::vector<Eigenvalues> m_eigenvalues;
  std::vector<double> m_sum_eigenvalues;
  std::vector<Point> m_centroids;
  std::vector<Vector> m_smallest_eigenvectors;
  std::vector<Vector> m_middle_eigenvectors;
  std::vector<Vector> m_largest_eigenvectors;


  
public:

  /// \cond SKIP_IN_MANUAL
  Local_eigen_analysis () { }
  /// \endcond

  /*! 

    \brief Computes the local eigen analysis of an input range based
    on a local neighborhood.

    \tparam NeighborQuery is a model of `NeighborQuery`
    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_map Property map to access the input points
    \param neighbor_query Object used to access neighborhoods of points

    \param mean_range The mean value of the range corresponding to the
    `knn` number of neighbors is returned by the constructor through
    this reference.
  */
  template <typename NeighborQuery>
  Local_eigen_analysis (RandomAccessIterator begin,
                        RandomAccessIterator end,
                        PointMap point_map,
                        const NeighborQuery& neighbor_query,
                        double& mean_range)
  {
    std::size_t size = end - begin;
    m_eigenvalues.reserve (size);
    m_centroids.reserve (size);
    m_smallest_eigenvectors.reserve (size);
    m_middle_eigenvectors.reserve (size);
    m_largest_eigenvectors.reserve (size);
    
    mean_range = 0.;
      
    for (std::size_t i = 0; i < size; i++)
      {
        std::vector<std::size_t> neighbors;
        neighbor_query (get(point_map, begin[i]), std::back_inserter (neighbors));

        std::vector<Point> neighbor_points;
        for (std::size_t j = 0; j < neighbors.size(); ++ j)
          neighbor_points.push_back (get(point_map, begin[neighbors[j]]));

        mean_range += CGAL::sqrt (CGAL::squared_distance (get(point_map, begin[i]),
                                                          get(point_map, begin[neighbors.back()])));
        
        compute (get(point_map, begin[i]), neighbor_points);
      }

    mean_range /= size;
  }

  /// \cond SKIP_IN_MANUAL
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
  /// \endcond

  /*!
    \brief Returns the estimated normal vector of the indexed point.
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
  const double& sum_eigenvalues (std::size_t index) const { return m_sum_eigenvalues[index]; }

};
  

}
  
}


#endif // CGAL_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H
