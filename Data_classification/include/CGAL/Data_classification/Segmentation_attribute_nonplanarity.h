#ifndef CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_NONPLANARITY_H
#define CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_NONPLANARITY_H

#include <vector>

#include <CGAL/Point_set_classification.h>

namespace CGAL {

template <typename Kernel>
class Segmentation_attribute_nonplanarity : public Segmentation_attribute
{
  typedef Point_set_classification<Kernel> PSC;

  std::vector<double> distance_to_plane_attribute;
  
public:
  double weight;
  double mean;
  double max;
  
  Segmentation_attribute_nonplanarity (PSC& M, double weight) : weight (weight)
  {
    for(int i=0; i<(int)M.HPS.size(); i++)
      distance_to_plane_attribute.push_back (std::sqrt (CGAL::squared_distance (M.HPS[i].position, M.planes[i])));

    this->compute_mean_max (distance_to_plane_attribute, mean, max);
    //    max *= 2;
  }

  virtual double value (int pt_index)
  {
    return std::max (0., std::min (1., distance_to_plane_attribute[pt_index] / weight));
  }

  virtual std::string id() { return "nonplanarity"; }
};


}

#endif // CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_NONPLANARITY_H
