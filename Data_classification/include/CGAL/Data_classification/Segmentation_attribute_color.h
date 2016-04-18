#ifndef CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_COLOR_H
#define CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_COLOR_H

#include <vector>

#include <CGAL/Point_set_classification.h>

namespace CGAL {

template <typename Kernel>
class Segmentation_attribute_color : public Segmentation_attribute
{
  typedef Point_set_classification<Kernel> PSC;
  typedef typename PSC::Image_float Image_float;
  typedef typename Data_classification::RGB_Color RGB_Color;
  typedef typename Data_classification::HSV_Color HSV_Color;
  
  std::vector<double> color_attribute;
  
public:
  double weight;
  double mean;
  double max;

  Segmentation_attribute_color (PSC& M,
                                double weight) : weight (weight)
  {

    std::cerr << "Using colors" << std::endl;
    for(int i=0;i<(int)M.HPS.size();i++)
      {
        HSV_Color c = Data_classification::rgb_to_hsv (M.HPS[i].color);
        color_attribute.push_back (std::exp (-(c[0] - 156.) * (c[0] - 156.) / (2. * 81. * 81.))
                                   * std::exp (-(c[1] - 5.) * (c[1] - 5.) / (2. * 4. * 4.))
                                   * std::exp (-(c[2] - 76.) * (c[2] - 76.) / (2. * 6.8 * 6.8)));
      }
    this->compute_mean_max (color_attribute, mean, max);
  }

  virtual double value (int pt_index)
  {
    return std::max (0., std::min (1., color_attribute[pt_index] / weight));
  }

  virtual std::string id() { return "color"; }
};


}

#endif // CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_COLOR_H
