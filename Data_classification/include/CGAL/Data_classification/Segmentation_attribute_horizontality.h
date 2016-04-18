#ifndef CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_HORIZONTALITY_H
#define CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_HORIZONTALITY_H

#include <vector>

#include <CGAL/Point_set_classification.h>

namespace CGAL {

template <typename Kernel>
class Segmentation_attribute_horizontality : public Segmentation_attribute
{
  typedef Point_set_classification<Kernel> PSC;

  std::vector<double> horizontality_attribute;
  
public:
  double weight;
  double mean;
  double max;
  
  Segmentation_attribute_horizontality (PSC& M, double weight) : weight (weight)
  {
    typename Kernel::Vector_3 vertical (0., 0., 1.);
    
    for(int i=0; i<(int)M.HPS.size(); i++)
      {
        typename Kernel::Vector_3 normal = M.planes[i].orthogonal_vector();
        normal = normal / std::sqrt (normal * normal);
        horizontality_attribute.push_back (1. - std::fabs(normal * vertical));
      }
    
    this->compute_mean_max (horizontality_attribute, mean, max);
    //    max *= 2;
  }

  virtual double value (int pt_index)
  {
    return std::max (0., std::min (1., horizontality_attribute[pt_index] / weight));
  }

  virtual std::string id() { return "horizontality"; }
};


}

#endif // CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_HORIZONTALITY_H
