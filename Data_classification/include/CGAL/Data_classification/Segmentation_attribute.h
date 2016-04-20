#ifndef CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_H
#define CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_H

#include <vector>

namespace CGAL {

class Segmentation_attribute
{
public:

  virtual ~Segmentation_attribute() { }
  
  virtual double value (std::size_t pt_index) = 0;

  virtual double favored (std::size_t pt_index) { return (1. - value (pt_index)); }
  virtual double penalized (std::size_t pt_index) { return value (pt_index); }
  //  virtual double ignored (std::size_t pt_index) { return std::min (favored(pt_index), penalized(pt_index)); }
  virtual double ignored (std::size_t) { return 0.5; }

  virtual std::string id() { return "abstract_attribute"; }

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

};

}

#endif // CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_H
