#ifndef CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_COLOR_H
#define CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_COLOR_H

#include <vector>

#include <CGAL/Point_set_classification.h>

namespace CGAL {

  /*!
    \ingroup PkgDataClassification

    \brief Segmentation attribute based on colorimetric information.

    If the input point cloud comes with colorimetric information, it
    can be useful for classification purposes. This attribute computes
    a distance between the color of a point and a user-defined color
    region in the HSV color-space, defined by a Gaussian
    distribution with user-defined mean and standard deviation values.

    \tparam Kernel The geometric kernel used.

  */
template <typename Kernel>
class Segmentation_attribute_color : public Segmentation_attribute
{
  
  typedef Point_set_classification<Kernel> PSC;
  typedef typename PSC::Image_float Image_float;
  typedef typename Data_classification::RGB_Color RGB_Color;
  typedef typename Data_classification::HSV_Color HSV_Color;
  
  std::vector<double> color_attribute;
  double mean_h, mean_s, mean_v, sd_h, sd_s, sd_v;
  
public:
  /// \cond SKIP_IN_MANUAL
  double weight;
  double mean;
  double max;
  /// \endcond

  /*!
    \brief Constructs an attribute based on the given color.

    \param M The point set classification object
    \param weight The relative weight of this attribute
    \param mean_h Mean hue of the selected color
    \param mean_s Mean saturation of the selected color
    \param mean_v Mean value of the selected color
    \param sd_h Standard deviation of the hue of the selected color
    \param sd_s Standard deviation of the saturation of the selected color
    \param sd_v Standard deviation of the value of the selected color

    \note The default values describe a gray color region
    corresponding to the color of a concrete road.
  */
  Segmentation_attribute_color (Point_set_classification<Kernel>& M,
                                double weight,
                                double mean_h = 156., double mean_s = 5., double mean_v = 76.,
                                double sd_h = 70., double sd_s = 12., double sd_v = 8.4)
    : mean_h (mean_h)
    , mean_s (mean_s)
    , mean_v (mean_v)
    , sd_h (sd_h)
    , sd_s (sd_s)
    , sd_v (sd_v)
    , weight (weight)
  {

    std::cerr << "Using colors " << mean_h << " " << mean_s << " " << mean_v << std::endl;
    for(std::size_t i=0; i < M.HPS.size();i++)
      {
        HSV_Color c = Data_classification::rgb_to_hsv (M.HPS[i].color);
        color_attribute.push_back (std::exp (-(c[0] - mean_h) * (c[0] - mean_h) / (2. * sd_h * sd_h))
                                   * std::exp (-(c[1] - mean_s) * (c[1] - mean_s) / (2. * sd_s * sd_s))
                                   * std::exp (-(c[2] - mean_v) * (c[2] - mean_v) / (2. * sd_v * sd_v)));
      }
    this->compute_mean_max (color_attribute, mean, max);
  }


  virtual double value (std::size_t pt_index)
  {
    return std::max (0., std::min (1., color_attribute[pt_index] / weight));
  }

  virtual std::string id() { return "color"; }

};


}

#endif // CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_COLOR_H
