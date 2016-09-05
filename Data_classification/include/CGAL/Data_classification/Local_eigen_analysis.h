#ifndef CGAL_DATA_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H
#define CGAL_DATA_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H

#include <vector>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Default_diagonalize_traits.h>
#include <CGAL/centroid.h>

#include <CGAL/Data_classification/Neighborhood.h>

namespace CGAL {

namespace Data_classification {

  /*!
    \ingroup PkgDataClassification

    \brief 

    \tparam Kernel The geometric kernel used.

  */
template <typename Kernel, typename RandomAccessIterator, typename PointPMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Local_eigen_analysis
{
public:
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Plane_3 Plane;

  typedef CGAL::cpp11::array<double, 3> Eigenvalues;

  typedef Data_classification::Neighborhood<Kernel, RandomAccessIterator, PointPMap> Neighborhood;

private:
  std::vector<Eigenvalues> eigenvalues;
  std::vector<Point> centroids;
  std::vector<Vector> smallest_eigenvectors;
  std::vector<Vector> middle_eigenvectors;
  std::vector<Vector> largest_eigenvectors;

  
public:

  Local_eigen_analysis () { }
  
  Local_eigen_analysis (RandomAccessIterator begin,
                        RandomAccessIterator end,
                        PointPMap point_pmap,
                        Neighborhood& neighborhood,
                        double radius_neighbors)
  {
    std::size_t size = end - begin;
    eigenvalues.reserve (size);
    centroids.reserve (size);
    smallest_eigenvectors.reserve (size);
    middle_eigenvectors.reserve (size);
    largest_eigenvectors.reserve (size);
    
    for (std::size_t i = 0; i < size; i++)
      {
        std::vector<std::size_t> neighbors;
        neighborhood.range_neighbors (i, radius_neighbors, std::back_inserter (neighbors));

        std::vector<Point> neighbor_points;
        for (std::size_t j = 0; j < neighbors.size(); ++ j)
          neighbor_points.push_back (get(point_pmap, begin[neighbors[j]]));
        
        compute (get(point_pmap, begin[i]), neighbor_points);
      }
  }

  void compute (const Point& query, std::vector<Point>& neighbor_points)
  {
    if (neighbor_points.size() == 0)
      {
        Eigenvalues v = {{ 0., 0., 0. }};
        eigenvalues.push_back (v);
        centroids.push_back (query);
        smallest_eigenvectors.push_back (Vector (0., 0., 1.));
        middle_eigenvectors.push_back (Vector (0., 1., 0.));
        largest_eigenvectors.push_back (Vector (1., 0., 0.));
        return;
      }

    Point centroid = CGAL::centroid (neighbor_points.begin(), neighbor_points.end());
    centroids.push_back (centroid);
    
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
    
    eigenvalues.push_back (evalues);
    smallest_eigenvectors.push_back (Vector (evectors[0], evectors[1], evectors[2]));
    middle_eigenvectors.push_back (Vector (evectors[3], evectors[4], evectors[5]));
    largest_eigenvectors.push_back (Vector (evectors[6], evectors[7], evectors[8]));
  }

  const Vector& normal_vector (std::size_t index) const { return smallest_eigenvectors[index]; }
  Plane plane (std::size_t index) const { return Plane (centroids[index], smallest_eigenvectors[index]); }

  const Eigenvalues& eigenvalue (std::size_t index) const { return eigenvalues[index]; }
};
  

}
  
}


#endif // CGAL_DATA_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H
