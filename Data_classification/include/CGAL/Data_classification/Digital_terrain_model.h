#ifndef CGAL_DATA_CLASSIFICATION_DIGITAL_TERRAIN_MODEL_H
#define CGAL_DATA_CLASSIFICATION_DIGITAL_TERRAIN_MODEL_H

#include <vector>

#include <CGAL/Data_classification/Image.h>
#include <CGAL/Data_classification/Planimetric_grid.h>

namespace CGAL {

namespace Data_classification {

  /*!
    \ingroup PkgDataClassification

    \brief 

    \tparam Kernel The geometric kernel used.

  */
template <typename Kernel, typename RandomAccessIterator, typename PointPMap>
class Digital_terrain_model
{
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Iso_cuboid_3 Iso_cuboid_3;
  typedef Image<float> Image_float;
  typedef Planimetric_grid<Kernel, RandomAccessIterator, PointPMap> Grid;
  
  Image_float m_dtm;
  
public:

  Digital_terrain_model () { }
  
  Digital_terrain_model (const RandomAccessIterator& begin,
                         const RandomAccessIterator& end,
                         PointPMap point_pmap,
                         Grid& grid,
                         const FT grid_resolution,
                         FT radius_dtm = -1.)
  {
    if (radius_dtm < 0.)
      radius_dtm = 5 * grid_resolution;

    std::size_t step=15;
    std::size_t size = end - begin;
    
    std::size_t square = (std::size_t)(radius_dtm/ grid_resolution) + 1;

    m_dtm = Image_float (grid.width(), grid.height());
    Image_float im_Zfront(grid.width(), grid.height());

    //round 1
    for (std::size_t j = 0; j < m_dtm.height(); j++)
      for (std::size_t i = 0; i < m_dtm.width(); i++)
        {
          m_dtm(i,j)=0;
          im_Zfront(i,j)=0;
        }

    for (std::size_t j = step/2+1; j < m_dtm.height(); j = j+step){	
      for (std::size_t i = step+1/2; i < m_dtm.width(); i = i+step){
		
        //storage of the points in the disk
        std::vector<float> list_pointsZ;
        std::size_t squareXmin = (i < square ? 0 : i - square);
        std::size_t squareXmax = std::min(m_dtm.width() - 1, i + square);
        std::size_t squareYmin = (j < square ? 0 : j - square);
        std::size_t squareYmax = std::min(m_dtm.height() - 1, j + square);

        for(std::size_t k = squareXmin; k <= squareXmax; k++){
          for(std::size_t l = squareYmin; l <= squareYmax; l++){
			
            double distance=sqrt(pow((double)i-k,2)+pow((double)j-l,2))*grid_resolution;
            if(distance<=radius_dtm){
              for(int nb=0;nb<(int)grid.indices(k,l).size();nb++)
                list_pointsZ.push_back(get(point_pmap, begin[grid.indices(k,l)[nb]]).z());
            }
          }
        }

        //ordering of the points
        float G1=0; 
        float G2=0;

        if(list_pointsZ.size()>0){
          std::sort(list_pointsZ.begin(),list_pointsZ.end());
          int ind_k2= (int)floor((double)(list_pointsZ.size()*0.6)-0.5);
          int ind_k1= (int)floor((double)(list_pointsZ.size()*0.3)-0.5);
          G1=list_pointsZ[ind_k1];
          G2=list_pointsZ[ind_k2];
        }
		
        float Gfront=(G1+G2)/2;	
        for(int iter=0;iter<3;iter++){
			
          float G1_temp=0; 
          float G2_temp=0;
          int count1=0; 

          for (std::size_t j=0;j<list_pointsZ.size();j++){
            if(list_pointsZ[j]<=Gfront){
              G1_temp+=list_pointsZ[j]; 
              count1++;
            }
            else{G2_temp+=list_pointsZ[j];}
          }

          if(count1>0) G1=(float)G1_temp/count1;
          if(count1<(int)list_pointsZ.size()) G2=(float)G2_temp/(list_pointsZ.size()-count1);
          Gfront=(G1+G2)/2;

        }

        m_dtm(i,j)=G1;
        im_Zfront(i,j)=Gfront;

        //extension by duplication
        std::size_t IsquareXmin = (i < step/2 ? 0 : i-step/2);
        std::size_t IsquareXmax = std::min(m_dtm.width()-1,i+step/2);
        std::size_t JsquareYmin = (j < step/2 ? 0 : j-step/2);
        std::size_t JsquareYmax = std::min(m_dtm.height()-1,j+step/2);
		
        if(m_dtm.width()-1-IsquareXmax<step) IsquareXmax=m_dtm.width()-1;
        if(m_dtm.height()-1-JsquareYmax<step) JsquareYmax=m_dtm.height()-1;
        if(IsquareXmin<step) IsquareXmin=0;
        if(JsquareYmin<step) JsquareYmin=0;

        for(std::size_t k = IsquareXmin; k <= IsquareXmax; k++){
          for(std::size_t l = JsquareYmin; l <= JsquareYmax; l++){
            if(grid.mask(k,l)){
              m_dtm(k,l)=G1;
              im_Zfront(k,l)=Gfront;
            }
          }
        }
      }
    }

    //Gaussian smoothness
    for(std::size_t ik = 0; ik < 10; ik++) {
      for(std::size_t j = 0; j < m_dtm.height(); j++){ 
        for (std::size_t i = 0; i < m_dtm.width(); i++) {
		
          std::size_t IsquareXmin = (i < 1 ? 0 : i-1);
          std::size_t IsquareXmax = std::min(m_dtm.width()-1,i+1); 
          std::size_t JsquareYmin = (j < 1 ? 0 : j-1);
          std::size_t JsquareYmax = std::min(m_dtm.height()-1,j+1); 
          int count=0; 

          for(std::size_t k = IsquareXmin; k <= IsquareXmax; k++){ 
            for(std::size_t l = JsquareYmin; l <= JsquareYmax; l++){ 
              if(grid.mask(k,l)) count++;
            }
          }	
			
          if(count==9) im_Zfront(i,j)=(im_Zfront(i-1,j-1)+2*im_Zfront(i,j-1)+im_Zfront(i+1,j-1)+2*im_Zfront(i-1,j)+4*im_Zfront(i,j)+2*im_Zfront(i+1,j)+im_Zfront(i-1,j+1)+2*im_Zfront(i,j+1)+im_Zfront(i+1,j+1))/16;
        }
      }
    }

    //round 2 
    std::vector < bool > test_ground;
    for (std::size_t i = 0; i < size; i ++)
      {
        std::size_t I = grid.x(i);
        std::size_t J = grid.y(i);
        bool test=false;
        if(get(point_pmap, begin[i]).z()-im_Zfront(I,J)<=0) {test=true; test_ground.push_back(test);}
        else{test_ground.push_back(test);}
    }

    for (int j=0;j<(int)m_dtm.height();j++){	
      for (int i=0;i<(int)m_dtm.width();i++){
        m_dtm(i,j)=0;
        im_Zfront(i,j)=0;
      }
    }

    for (std::size_t j = step/2+1; j < m_dtm.height(); j = j+step){	
      for (std::size_t i = step/2+1; i < m_dtm.width(); i = i+step){
		
        //stockage des nouveaux points
        std::vector < float > list_pointsZ;
        std::size_t squareXmin = (i < square ? 0 : i-square);
        std::size_t squareXmax = std::min(m_dtm.width()-1,i+square);
        std::size_t squareYmin = (j < square ? 0 : j-square);
        std::size_t squareYmax = std::min(m_dtm.height()-1,j+square);

        for(std::size_t k = squareXmin; k <= squareXmax; k++){
          for(std::size_t l = squareYmin; l <= squareYmax; l++){
				
            double distance=sqrt(pow((double)i-k,2)+pow((double)j-l,2))*grid_resolution;
				
            if(distance<=radius_dtm){
					
              for(int nb=0;nb<(int)grid.indices(k,l).size();nb++){
                if(test_ground[grid.indices(k,l)[nb]])
                  list_pointsZ.push_back(get(point_pmap, begin[grid.indices(k,l)[nb]]).z());
              }
            }
          }
        }

        //ordering the new points
        float G1=0; 
        float G2=0;

        if(list_pointsZ.size()>0){
          std::sort(list_pointsZ.begin(),list_pointsZ.end());
          int ind_k2= (int)floor((double)(list_pointsZ.size()*0.6)-0.5);
          int ind_k1= (int)floor((double)(list_pointsZ.size()*0.3)-0.5);
          G1=list_pointsZ[ind_k1];
          G2=list_pointsZ[ind_k2];
        }

        float Gfront=(G1+G2)/2;	
        for(int iter=0;iter<3;iter++){
          float G1_temp=0; 
          float G2_temp=0;
          int count1=0; 

          for (int j=0;j<(int)list_pointsZ.size();j++){
            if(list_pointsZ[j]<=Gfront){G1_temp+=list_pointsZ[j]; count1++;}
            else{G2_temp+=list_pointsZ[j];}
          }
          if(count1>0) G1=(float)G1_temp/count1;
          if(count1<(int)list_pointsZ.size()) G2=(float)G2_temp/(list_pointsZ.size()-count1);
          Gfront=(G1+G2)/2;
        }
        m_dtm(i,j)=G1;
        im_Zfront(i,j)=Gfront;

        //extension by duplication
        std::size_t IsquareXmin = (i < step/2 ? 0 : i-step/2);
        std::size_t IsquareXmax = std::min(m_dtm.width()-1,i+step/2);
        std::size_t JsquareYmin = (j < step/2 ? 0 : j-step/2);
        std::size_t JsquareYmax = std::min(m_dtm.height()-1,j+step/2);
		
        if(m_dtm.width()-1-IsquareXmax<step) IsquareXmax=m_dtm.width()-1;
        if(m_dtm.height()-1-JsquareYmax<step) JsquareYmax=m_dtm.height()-1;
        if(IsquareXmin<step) IsquareXmin=0;
        if(JsquareYmin<step) JsquareYmin=0;

        for(std::size_t k = IsquareXmin; k <= IsquareXmax; k++){
          for(std::size_t l = JsquareYmin; l <= JsquareYmax; l++){
            if(grid.mask(k,l)){
              m_dtm(k,l)=G1; 
              im_Zfront(k,l)=Gfront;
            }
          }
        }
      }
    }

    //Gaussian smoothness
    for(std::size_t ik = 0; ik < 10; ik++) {
      for(std::size_t j = 0; j < m_dtm.height(); j++){ 
        for (std::size_t i = 0; i < m_dtm.width(); i++) {
		
          std::size_t IsquareXmin = (i < 1 ? 0 : i - 1);
          std::size_t IsquareXmax = std::min(m_dtm.width()-1,i+1); 
          std::size_t JsquareYmin = (j < 1 ? 0 : j - 1);
          std::size_t JsquareYmax = std::min(m_dtm.height()-1,j+1); 
          int count=0; 

          for(std::size_t k=IsquareXmin; k<=IsquareXmax; k++){ 
            for(std::size_t l=JsquareYmin; l<=JsquareYmax; l++){ 
              if(grid.mask(k,l)) count++;
            }
          }	
          if(count==9) m_dtm(i,j)=(m_dtm(i-1,j-1)+2*m_dtm(i,j-1)+m_dtm(i+1,j-1)+2*m_dtm(i-1,j)+4*m_dtm(i,j)+2*m_dtm(i+1,j)+m_dtm(i-1,j+1)+2*m_dtm(i,j+1)+m_dtm(i+1,j+1))/16;
        }
      }
    }

    
  }
  
};
  

}
  
}


#endif // CGAL_DATA_CLASSIFICATION_DIGITAL_TERRAIN_MODEL_H
