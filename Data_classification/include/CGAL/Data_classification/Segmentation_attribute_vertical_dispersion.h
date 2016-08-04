#ifndef CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_VERTICAL_DISPERSION_H
#define CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_VERTICAL_DISPERSION_H

#include <vector>

#include <CGAL/Data_classification/Image.h>
#include <CGAL/Data_classification/Planimetric_grid.h>

namespace CGAL {

  /*!
    \ingroup PkgDataClassification

    \brief Segmentation attribute based on local vertical dispersion of points.

    Urban scenes can usually be described as a set of 2D regions with
    different heights. While these heights are usually piecewise
    constant or piecewise linear, on some specific parts of the scene
    such as vegetation, they can become extremely unstable. This
    attribute quantifies the vertical dispersion of the points on a local
    Z-cylinder around the points.

    \tparam Kernel The geometric kernel used.

  */
template <typename Kernel, typename RandomAccessIterator, typename PointPMap>
class Segmentation_attribute_vertical_dispersion : public Segmentation_attribute
{
  typedef Data_classification::Image<float> Image_float;
  typedef Data_classification::Planimetric_grid<Kernel, RandomAccessIterator, PointPMap> Grid;
  std::vector<double> vertical_dispersion;
  
public:
  /// \cond SKIP_IN_MANUAL
  double weight;
  double mean;
  double max;
  /// \endcond

  /*!
    \brief Constructs the attribute.
  */
  Segmentation_attribute_vertical_dispersion (RandomAccessIterator begin,
                                              RandomAccessIterator end,
                                              PointPMap point_pmap,
                                              Grid& grid,
                                              const double grid_resolution,
                                              double radius_neighbors = -1.,
                                              double weight = 1.) : weight (weight)
  {
    if (radius_neighbors < 0.)
      radius_neighbors = 5. * grid_resolution;
    
    Image_float Dispersion(grid.width(), grid.height());
    for (std::size_t j = 0; j < grid.height(); j++)	
      for (std::size_t i = 0; i < grid.width(); i++)
        Dispersion(i,j)=0;
    
    std::size_t square = (std::size_t)(0.5 * radius_neighbors / grid_resolution) + 1;
    typename Kernel::Vector_3 verti (0., 0., 1.);
    
    for (std::size_t j = 0; j < grid.height(); j++){	
      for (std::size_t i = 0; i < grid.width(); i++){
						
        if(!(grid.mask(i,j)))
          continue;
        std::vector<double> hori;
            
        std::size_t squareXmin = (i < square ? 0 : i - square);
        std::size_t squareXmax = std::min (grid.width()-1, i + square);
        std::size_t squareYmin = (j < square ? 0 : j - square);
        std::size_t squareYmax = std::min (grid.height()-1, j + square);

        for(std::size_t k = squareXmin; k <= squareXmax; k++)
          for(std::size_t l = squareYmin; l <= squareYmax; l++)
            if(sqrt(pow((double)k-i,2)+pow((double)l-j,2))
               <=(double)0.5*radius_neighbors/grid_resolution
               && (grid.indices(k,l).size()>0))
              for(int t=0; t<(int)grid.indices(k,l).size();t++)
                {
                  int ip = grid.indices(k,l)[t];
                  hori.push_back (get(point_pmap, begin[ip]).z());
                }
        if (hori.empty())
          continue;
              
        std::sort (hori.begin(), hori.end());

        std::size_t nb_layers = 1;

        std::vector<bool> occupy (1 + (std::size_t)((hori.back() - hori.front()) / grid_resolution), false);
              
        std::size_t last_index = 0;
        for (std::size_t k = 0; k < hori.size(); ++ k)
          {
            std::size_t index = (std::size_t)((hori[k] - hori.front()) / grid_resolution);
            occupy[index] = true;
            if (index > last_index + 1)
              ++ nb_layers;
            last_index = index;
          }

        std::size_t nb_occ = 0;
        for (std::size_t k = 0; k < occupy.size(); ++ k)
          if (occupy[k])
            ++ nb_occ;
					
        Dispersion(i,j)= 1.f - (nb_occ / (float)(occupy.size()));
			
      }
		
    }
    for (std::size_t i = 0; i < (std::size_t)(end - begin);i++)
      {
        int I= grid.x(i);
        int J= grid.y(i);
        vertical_dispersion.push_back((double)Dispersion(I,J));
      }

    this->compute_mean_max (vertical_dispersion, mean, max);
  }

  virtual double value (std::size_t pt_index)
  {
    return std::max (0., std::min (1., vertical_dispersion[pt_index] / weight));
  }

  virtual std::string id() { return "vertical_dispersion"; }
};

}

#endif // CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_VERTICAL_DISPERSION_H
