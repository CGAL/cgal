#ifndef CGAL_DATA_CLASSIFICATION_ATTRIBUTE_ECHO_SCATTER_H
#define CGAL_DATA_CLASSIFICATION_ATTRIBUTE_ECHO_SCATTER_H

#include <vector>


namespace CGAL {

namespace Data_classification {
  
  /*!
    \ingroup PkgDataClassification

    \brief Attribute based on echo scatter.

    The number of returns (echo number) is a useful information
    provided by most LIDAR sensor. It can help identifying trees.

    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
    \tparam EchoMap is a model of `ReadablePropertyMap` with value type `std::size_t`.
  */
template <typename Kernel, typename RandomAccessIterator, typename PointMap, typename EchoMap>
class Attribute_echo_scatter : public Attribute
{
  typedef Data_classification::Image<float> Image_float;
  typedef Data_classification::Planimetric_grid<Kernel, RandomAccessIterator, PointMap> Grid;

  std::vector<double> echo_scatter;
  
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param echo_map Property map to access the echo values of the input points
    \param grid Precomputed `Planimetric_grid`
    \param grid_resolution Resolution of the planimetric grid
    \param radius_neighbors Radius of local neighborhoods
  */
  Attribute_echo_scatter (RandomAccessIterator begin,
                          RandomAccessIterator end,
                          EchoMap echo_map,
                          Grid& grid,
                          const double grid_resolution,
                          double radius_neighbors = 1.)
  {
    this->weight = 1.;
    Image_float Scatter(grid.width(), grid.height());
    for (std::size_t j = 0; j < grid.height(); j++)
      for (std::size_t i = 0; i < grid.width(); i++)
        Scatter(i,j)=0;

    std::size_t square = (std::size_t)(0.5 * radius_neighbors / grid_resolution) + 1;

    for (std::size_t j = 0; j < grid.height(); j++){	
      for (std::size_t i = 0; i < grid.width(); i++){
						
        if(grid.mask(i,j)){

          std::size_t squareXmin = (i < square ? 0 : i - square);
          std::size_t squareXmax = (std::min) (grid.width()-1, i + square);
          std::size_t squareYmin = (j < square ? 0 : j - square);
          std::size_t squareYmax = (std::min) (grid.height()-1, j + square);
			
          std::size_t NB_echo_sup=0;
          std::size_t NB_echo_total=0;

          for(std::size_t k = squareXmin; k <= squareXmax; k++){
            for(std::size_t l = squareYmin; l <= squareYmax; l++){
									
              if(CGAL::sqrt(pow((double)k-i,2)+pow((double)l-j,2))<=(double)0.5*radius_neighbors/grid_resolution){
										
                if(grid.indices(k,l).size()>0){
									
                  for(std::size_t t=0; t<grid.indices(k,l).size();t++){
												
                    std::size_t ip = grid.indices(k,l)[t]; 
                    if(get(echo_map, begin[ip]) > 1)
                      NB_echo_sup++;
                  }
									
                  NB_echo_total=NB_echo_total+grid.indices(k,l).size();
									
                }
							
              }
						
            }
					
          }
					
          Scatter(i,j)=(float)NB_echo_sup/NB_echo_total;
				
        }
			
      }
		
    }
    for(std::size_t i = 0; i < (std::size_t)(end - begin); i++){
      std::size_t I= grid.x(i);
      std::size_t J= grid.y(i);
      echo_scatter.push_back((double)Scatter(I,J));
    }
    this->compute_mean_max (echo_scatter, this->mean, this->max);
  }

  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return echo_scatter[pt_index];
  }

  virtual std::string id() { return "echo_scatter"; }
  /// \endcond
};

} // namespace Data_classification
  
} // namespace CGAL

#endif // CGAL_DATA_CLASSIFICATION_ATTRIBUTE_ECHO_SCATTER_H
