#ifndef CGAL_DATA_CLASSIFICATION_ATTRIBUTE_H
#define CGAL_DATA_CLASSIFICATION_ATTRIBUTE_H

#include <vector>

namespace CGAL {

namespace Data_classification {
  
/*!
  \ingroup PkgDataClassification

  \brief Abstract class describing a segmentation attribute.

  A segmentation attribute associates a scalar value to each point
  index of the classification input.
*/


class Attribute
{
public:
  /// \cond SKIP_IN_MANUAL
  double weight;
  double mean;
  double max;

  virtual ~Attribute() { }
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
  virtual double normalized (std::size_t pt_index)
  {
    return std::max (0., std::min (1., value(pt_index) / weight));
  }
  virtual double favored (std::size_t pt_index) { return (1. - normalized (pt_index)); }
  virtual double penalized (std::size_t pt_index) { return normalized (pt_index); }
  //  virtual double ignored (std::size_t pt_index) { return std::min (favored(pt_index), penalized(pt_index)); }
  virtual double ignored (std::size_t) { return 0.5; }

  void compute_mean_max (std::vector<double>& vect, double& mean, double& max)
  {
    mean = 0.;
    max = -std::numeric_limits<double>::max();
    double min = std::numeric_limits<double>::max();
    
    for (std::size_t i = 0; i < vect.size(); ++ i)
      {
        mean += vect[i];
        if (vect[i] > max)
          max = vect[i];
        if (vect[i] < min)
          min = vect[i];
      }
    //    std::cerr << id() << " Min/max = " << min << " / " << max << std::endl;
    mean /= vect.size();

  }
  /// \endcond
};


typedef boost::shared_ptr<Attribute> Attribute_handle;
  
} // namespace Data_classification

} // namespace CGAL

#endif // CGAL_DATA_CLASSIFICATION_ATTRIBUTE_H
