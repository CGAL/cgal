#ifndef CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_H
#define CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_H

#include <vector>

namespace CGAL {


/*!
  \ingroup PkgDataClassification

  \brief Abstract class describing a segmentation attribute.

  A segmentation attribute must have a unique ID and associate a
  scalar value to each point index of the classification input.
*/


class Segmentation_attribute
{
public:
  /// \cond SKIP_IN_MANUAL
  virtual ~Segmentation_attribute() { }
  /// \endcond

  /*!
    \brief Get the ID of the attribute.
    \return the ID of the attribute
  */
  virtual std::string id() { return "abstract_attribute"; }

  /*!
    \brief Value taken by the attribute on the given point.

    This method must be implemented by inherited classes.

    \param pt_index Index of the point in the classification input

    \return the value of the attribute on the point of index `pt_index`
  */
  virtual double value (std::size_t pt_index) = 0;

  /// \cond SKIP_IN_MANUAL
  virtual double favored (std::size_t pt_index) { return (1. - value (pt_index)); }
  virtual double penalized (std::size_t pt_index) { return value (pt_index); }
  //  virtual double ignored (std::size_t pt_index) { return std::min (favored(pt_index), penalized(pt_index)); }
  virtual double ignored (std::size_t) { return 0.5; }

  void compute_mean_max (std::vector<double>& vect, double& mean, double& max)
  {
    mean = 0.;
    max = 0.;
  
    for (std::size_t i = 0; i < vect.size(); ++ i)
      {
        mean += vect[i];
        if (vect[i] > max)
          max = vect[i];
      }
    mean /= vect.size();

  }
  /// \endcond
};

}

#endif // CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_H
