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
  
  Segmentation_attribute_nonplanarity (PSC& M, double weight, bool on_groups = false) : weight (weight)
  {
    if (on_groups)
      {
        std::vector<double> sq_dist (M.groups.size(), 0.);
        std::vector<std::size_t> nb (M.groups.size(), 0);
        for (std::size_t i = 0; i < M.HPS.size(); ++ i)
          {
            std::size_t g = M.HPS[i].group;
            if (g == (std::size_t)(-1))
              continue;
            sq_dist[g] += CGAL::squared_distance (M.groups[g], M.HPS[i].position);
            nb[g] ++;
          }
        for (std::size_t i = 0; i < sq_dist.size(); ++ i)
          if (nb[i] > 0)
            //            sq_dist[i] = std::sqrt (sq_dist[i] / nb[i]);
            sq_dist[i] = 1. / nb[i];
        for (std::size_t i = 0; i < M.HPS.size(); ++ i)
          {
            std::size_t g = M.HPS[i].group;
            if (g == (std::size_t)(-1))
              distance_to_plane_attribute.push_back (0.);
            else
              distance_to_plane_attribute.push_back (sq_dist[g]);
          }
      }
    else
      {
        for(int i=0; i<(int)M.HPS.size(); i++)
          distance_to_plane_attribute.push_back (std::sqrt (CGAL::squared_distance (M.HPS[i].position, M.planes[i])));
      }
    
    this->compute_mean_max (distance_to_plane_attribute, mean, max);
    //    max *= 2;
  }

  virtual double value (std::size_t pt_index)
  {
    return std::max (0., std::min (1., distance_to_plane_attribute[pt_index] / weight));
  }

  virtual std::string id() { return "nonplanarity"; }
};


}

#endif // CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_NONPLANARITY_H
