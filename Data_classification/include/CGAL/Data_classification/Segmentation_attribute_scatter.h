#ifndef CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_SCATTER_H
#define CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_SCATTER_H

#include <vector>

#include <CGAL/Point_set_classification.h>

namespace CGAL {


template <typename Kernel>
class Segmentation_attribute_scatter : public Segmentation_attribute
{
  typedef Point_set_classification<Kernel> PSC;
  typedef typename PSC::Image_float Image_float;
  
  std::vector<double> vegetation_attribute;
  
public:
  double weight;
  double mean;
  double max;

  Segmentation_attribute_scatter (PSC& M,
                                  double weight) : weight (weight)
  {
    Image_float Vegetation(M.grid_HPS.width(), M.grid_HPS.height());

    for (int j=0;j<(int)M.DTM.height();j++)	
      for (int i=0;i<(int)M.DTM.width();i++)
        Vegetation(i,j)=0;
		
    std::size_t square = (std::size_t)(0.5 * M.m_radius_neighbors / M.m_grid_resolution) + 1;

    if(M.is_echo_given){
      std::cerr << "Echo given" << std::endl;				
      for (std::size_t j = 0; j < M.grid_HPS.height(); j++){	
        for (std::size_t i = 0; i < M.grid_HPS.width(); i++){
						
          if(M.Mask(i,j)){

            std::size_t squareXmin = (i < square ? 0 : i - square);
            std::size_t squareXmax = std::min (M.grid_HPS.width()-1, i + square);
            std::size_t squareYmin = (j < square ? 0 : j - square);
            std::size_t squareYmax = std::min (M.grid_HPS.height()-1, j + square);
			
            int NB_echo_sup=0;
            int NB_echo_total=0;

            for(std::size_t k = squareXmin; k <= squareXmax; k++){
              for(std::size_t l = squareYmin; l <= squareYmax; l++){
									
                if(sqrt(pow((double)k-i,2)+pow((double)l-j,2))<=(double)0.5*M.m_radius_neighbors/M.m_grid_resolution){
										
                  if(M.grid_HPS(k,l).size()>0){
									
                    for(int t=0; t<(int)M.grid_HPS(k,l).size();t++){
												
                      int ip = M.grid_HPS(k,l)[t]; 
                      if(M.HPS[ip].echo>1) NB_echo_sup++;
                    }
									
                    NB_echo_total=NB_echo_total+M.grid_HPS(k,l).size();
									
                  }
							
                }
						
              }
					
            }
					
            Vegetation(i,j)=(float)NB_echo_sup/NB_echo_total;
				
          }
			
        }
		
      }
      for(int i=0;i<(int)M.HPS.size();i++){
        int I= M.HPS[i].ind_x;
        int J= M.HPS[i].ind_y;
        vegetation_attribute.push_back((double)Vegetation(I,J));
      }
	
    }
    else 
      {
        typename Kernel::Vector_3 verti (0., 0., 1.);
        
        for (std::size_t j = 0; j < M.grid_HPS.height(); j++){	
          for (std::size_t i = 0; i < M.grid_HPS.width(); i++){
						
            if(!(M.Mask(i,j)))
              continue;
            std::vector<double> hori;
            
            std::size_t squareXmin = (i < square ? 0 : i - square);
            std::size_t squareXmax = std::min (M.grid_HPS.width()-1, i + square);
            std::size_t squareYmin = (j < square ? 0 : j - square);
            std::size_t squareYmax = std::min (M.grid_HPS.height()-1, j + square);

            for(std::size_t k = squareXmin; k <= squareXmax; k++)
              for(std::size_t l = squareYmin; l <= squareYmax; l++)
                if(sqrt(pow((double)k-i,2)+pow((double)l-j,2))
                   <=(double)0.5*M.m_radius_neighbors/M.m_grid_resolution
                   && (M.grid_HPS(k,l).size()>0))
                  for(int t=0; t<(int)M.grid_HPS(k,l).size();t++)
                    {
                      int ip = M.grid_HPS(k,l)[t];
                      hori.push_back (M.HPS[ip].position.z());
                    }
            if (hori.empty())
              continue;
              
            std::sort (hori.begin(), hori.end());

            std::size_t nb_layers = 1;

            std::vector<bool> occupy (1 + (std::size_t)((hori.back() - hori.front()) / M.m_grid_resolution), false);
              
            std::size_t last_index = 0;
            for (std::size_t k = 0; k < hori.size(); ++ k)
              {
                std::size_t index = (std::size_t)((hori[k] - hori.front()) / M.m_grid_resolution);
                occupy[index] = true;
                if (index > last_index + 1)
                  ++ nb_layers;
                last_index = index;
              }

            std::size_t nb_occ = 0;
            for (std::size_t k = 0; k < occupy.size(); ++ k)
              if (occupy[k])
                ++ nb_occ;
					
            Vegetation(i,j)= 1.f - (nb_occ / (float)(occupy.size()));
			
          }
		
        }
        for(int i=0;i<(int)M.HPS.size();i++){
          int I= M.HPS[i].ind_x;
          int J= M.HPS[i].ind_y;
          vegetation_attribute.push_back((double)Vegetation(I,J));
        }
        
    }

    this->compute_mean_max (vegetation_attribute, mean, max);
    //    max *= 2;
  }

  virtual double value (int pt_index)
  {
    return std::max (0., std::min (1., vegetation_attribute[pt_index] / weight));
  }

  virtual std::string id() { return "scatter"; }
};

}

#endif // CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_SCATTER_H
