#ifndef CGAL_DATA_CLASSIFICATION_ATTRIBUTE_DISTANCE_TO_PLANE_H
#define CGAL_DATA_CLASSIFICATION_ATTRIBUTE_DISTANCE_TO_PLANE_H

#include <vector>

#include <CGAL/Point_set_classification.h>

namespace CGAL {

namespace Data_classification {

  /*!
    \ingroup PkgDataClassification

    \brief Segmentation attribute based on local distance to a fitted plane.

    Characterizing a level of non-planarity can help identify noisy
    parts of the input such as vegetation. This attribute computes the
    distance of a point to a locally fitted plane.
    
    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_pmap Property map to access the input points
    \tparam DiagonalizeTraits Solver used for matrix diagonalization.
  */
template <typename Kernel, typename RandomAccessIterator, typename PointPMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Attribute_distance_to_plane : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointPMap, DiagonalizeTraits> Local_eigen_analysis;

  std::vector<double> distance_to_plane_attribute;
  
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_pmap Property map to access the input points
    \param eigen Class with precompute eigenvectors and eigenvalues
    \param weight Weight of the attribute
  */
  Attribute_distance_to_plane (RandomAccessIterator begin,
                               RandomAccessIterator end,
                               PointPMap point_pmap,
                               const Local_eigen_analysis& eigen,
                               double weight = 1.)
  {
    this->weight = weight;
    for(std::size_t i = 0; i < (std::size_t)(end - begin); i++)
      distance_to_plane_attribute.push_back
        (std::sqrt (CGAL::squared_distance (get(point_pmap, begin[i]), eigen.plane(i))));
    
    this->compute_mean_max (distance_to_plane_attribute, this->mean, this->max);
    //    max *= 2;
  }

  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return distance_to_plane_attribute[pt_index];
  }

  virtual std::string id() { return "distance_to_plane"; }
  /// \endcond
};


} // namespace Data_classification
  
} // namespace CGAL

#endif // CGAL_DATA_CLASSIFICATION_ATTRIBUTE_DISTANCE_TO_PLANE_H
