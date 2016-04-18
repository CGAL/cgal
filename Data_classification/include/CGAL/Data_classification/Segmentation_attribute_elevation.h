#ifndef CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_ELEVATION_H
#define CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_ELEVATION_H

#include <vector>

#include <CGAL/Point_set_classification.h>

namespace CGAL {

template <typename Kernel>
class Segmentation_attribute_elevation : public Segmentation_attribute
{
  typedef Point_set_classification<Kernel> PSC;
  typedef typename PSC::Image_float Image_float;
  
  std::vector<double> elevation_attribute;
  
public:
  double weight;
  double mean;
  double max;
  
  Segmentation_attribute_elevation (PSC& M, double weight) : weight (weight)
  {
    //DEM
    Image_float DEM(M.grid_HPS.width(),M.grid_HPS.height());
    Image_float DEMtemp(M.grid_HPS.width(),M.grid_HPS.height());

    for (int j=0;j<(int)DEM.height();j++){	
      for (int i=0;i<(int)DEM.width();i++){
		
        double mean_height=0;
        std::vector < double > list_Z;
		
        for(std::size_t k=0;k<M.grid_HPS(i,j).size();k++) list_Z.push_back(M.HPS[M.grid_HPS(i,j)[k]].position.z());
			
        if(list_Z.size()>0){
			
          std::sort(list_Z.begin(),list_Z.end());
          int ind_k= (int)floor((double)(list_Z.size()*0.9)-0.5);
          mean_height=list_Z[ind_k];
        }

        DEM(i,j)=(float)mean_height;
        DEMtemp(i,j)=DEM(i,j);
      }
    }

    std::size_t square = (std::size_t)(M.m_radius_neighbors / M.m_grid_resolution) + 1;

    for (std::size_t j = 0; j < DEM.height(); j++){	
      for (std::size_t i = 0; i < DEM.width(); i++){
		
        if((M.Mask(i,j))&&(DEM(i,j)==0)){
				
          double distance_tot=0;	
          double val=0;
          std::size_t squareXmin = (i < square ? 0 : i - square);
          std::size_t squareXmax = std::min(DEM.width() - 1, i + square);
          std::size_t squareYmin = (j < square ? 0 : j - square);
          std::size_t squareYmax = std::min(DEM.height() - 1, j + square);
			
          for(std::size_t k = squareXmin; k <= squareXmax; k++){
            for(std::size_t l = squareYmin; l <= squareYmax; l++){
					
              double distance=sqrt(pow((double)i-k,2)+pow((double)j-l,2))*M.m_grid_resolution;
					
              if((distance<=M.m_radius_neighbors)&&(DEM(k,l)>0)&&(distance!=0)){
                double dista=distance*distance;
                val=val+DEM(k,l)/dista;
                distance_tot=distance_tot+1/dista;
              }
            }
          }
			
          val=val/distance_tot;
          DEMtemp(i,j)=(float)val;			
        }
      }
    }

    float scale_byte_min_dem = M.BBox_scan.zmax();
    float scale_byte_max_dem = M.BBox_scan.zmin();

    for (std::size_t j = 0; j < DEM.height(); j++){	
      for (std::size_t i = 0; i < DEM.width(); i++){
		
        DEM(i,j)=DEMtemp(i,j);

        if(M.Mask(i,j)){
          if(DEM(i,j)>scale_byte_max_dem) scale_byte_max_dem=DEM(i,j);
          if(DEM(i,j)<scale_byte_min_dem) scale_byte_min_dem=DEM(i,j);
        }
      }
    }

    //TODO: smooth the DEM
    for (std::size_t j = 1; j < DEM.height()-1; j++){	
      for (std::size_t i = 1; i < DEM.width()-1; i++){
        if(M.Mask(i,j)){
          //DEM(i,j)=
        }
      }
    }

    //DTM computation
    std::size_t step=15;
    square = (std::size_t)(M.m_radius_dtm/M.m_grid_resolution)+1;
    Image_float toto(M.grid_HPS.width(),M.grid_HPS.height());
    Image_float im_Zfront(M.grid_HPS.width(),M.grid_HPS.height());
    M.DTM=toto;

    //round 1
    for (std::size_t j = 0; j < M.DTM.height(); j++){	
      for (std::size_t i = 0; i < M.DTM.width(); i++){
        M.DTM(i,j)=0;
        im_Zfront(i,j)=0;
      }
    }


    for (std::size_t j = step/2+1; j < M.DTM.height(); j = j+step){	
      for (std::size_t i = step+1/2; i < M.DTM.width(); i = i+step){
		
        //storage of the points in the disk
        std::vector<float> list_pointsZ;
        std::size_t squareXmin = (i < square ? 0 : i - square);
        std::size_t squareXmax = std::min(M.DTM.width() - 1, i + square);
        std::size_t squareYmin = (j < square ? 0 : j - square);
        std::size_t squareYmax = std::min(M.DTM.height() - 1, j + square);

        for(std::size_t k = squareXmin; k <= squareXmax; k++){
          for(std::size_t l = squareYmin; l <= squareYmax; l++){
			
            double distance=sqrt(pow((double)i-k,2)+pow((double)j-l,2))*M.m_grid_resolution;
            if(distance<=M.m_radius_dtm){
              for(int nb=0;nb<(int)M.grid_HPS(k,l).size();nb++) list_pointsZ.push_back(M.HPS[M.grid_HPS(k,l)[nb]].position.z());
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

        M.DTM(i,j)=G1;
        im_Zfront(i,j)=Gfront;

        //extension by duplication
        std::size_t IsquareXmin = (i < step/2 ? 0 : i-step/2);
        std::size_t IsquareXmax = std::min(M.DTM.width()-1,i+step/2);
        std::size_t JsquareYmin = (j < step/2 ? 0 : j-step/2);
        std::size_t JsquareYmax = std::min(M.DTM.height()-1,j+step/2);
		
        if(M.DTM.width()-1-IsquareXmax<step) IsquareXmax=M.DTM.width()-1;
        if(M.DTM.height()-1-JsquareYmax<step) JsquareYmax=M.DTM.height()-1;
        if(IsquareXmin<step) IsquareXmin=0;
        if(JsquareYmin<step) JsquareYmin=0;

        for(std::size_t k = IsquareXmin; k <= IsquareXmax; k++){
          for(std::size_t l = JsquareYmin; l <= JsquareYmax; l++){
            if(M.Mask(k,l)){
              M.DTM(k,l)=G1;
              im_Zfront(k,l)=Gfront;
            }
          }
        }
      }
    }

    //Gaussian smoothness
    for(std::size_t ik = 0; ik < 10; ik++) {
      for(std::size_t j = 0; j < M.DTM.height(); j++){ 
        for (std::size_t i = 0; i < M.DTM.width(); i++) {
		
          std::size_t IsquareXmin = (i < 1 ? 0 : i-1);
          std::size_t IsquareXmax = std::min(M.DTM.width()-1,i+1); 
          std::size_t JsquareYmin = (j < 1 ? 0 : j-1);
          std::size_t JsquareYmax = std::min(M.DTM.height()-1,j+1); 
          int count=0; 

          for(std::size_t k = IsquareXmin; k <= IsquareXmax; k++){ 
            for(std::size_t l = JsquareYmin; l <= JsquareYmax; l++){ 
              if(M.Mask(k,l)) count++;
            }
          }	
			
          if(count==9) im_Zfront(i,j)=(im_Zfront(i-1,j-1)+2*im_Zfront(i,j-1)+im_Zfront(i+1,j-1)+2*im_Zfront(i-1,j)+4*im_Zfront(i,j)+2*im_Zfront(i+1,j)+im_Zfront(i-1,j+1)+2*im_Zfront(i,j+1)+im_Zfront(i+1,j+1))/16;
        }
      }
    }

    //round 2 
    std::vector < bool > test_ground;
    for(int i=0;i<(int)M.HPS.size();i++){
      int I=M.HPS[i].ind_x;
      int J=M.HPS[i].ind_y;
      bool test=false;
      if(M.HPS[i].position.z()-im_Zfront(I,J)<=0) {test=true; test_ground.push_back(test);}
      else{test_ground.push_back(test);}
    }

    for (int j=0;j<(int)M.DTM.height();j++){	
      for (int i=0;i<(int)M.DTM.width();i++){
        M.DTM(i,j)=0;
        im_Zfront(i,j)=0;
      }
    }

    for (std::size_t j = step/2+1; j < M.DTM.height(); j = j+step){	
      for (std::size_t i = step/2+1; i < M.DTM.width(); i = i+step){
		
        //stockage des nouveaux points
        std::vector < float > list_pointsZ;
        std::size_t squareXmin = (i < square ? 0 : i-square);
        std::size_t squareXmax = std::min(M.DTM.width()-1,i+square);
        std::size_t squareYmin = (j < square ? 0 : j-square);
        std::size_t squareYmax = std::min(M.DTM.height()-1,j+square);

        for(std::size_t k = squareXmin; k <= squareXmax; k++){
          for(std::size_t l = squareYmin; l <= squareYmax; l++){
				
            double distance=sqrt(pow((double)i-k,2)+pow((double)j-l,2))*M.m_grid_resolution;
				
            if(distance<=M.m_radius_dtm){
					
              for(int nb=0;nb<(int)M.grid_HPS(k,l).size();nb++){
                if(test_ground[M.grid_HPS(k,l)[nb]]) list_pointsZ.push_back(M.HPS[M.grid_HPS(k,l)[nb]].position.z());
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
        M.DTM(i,j)=G1;
        im_Zfront(i,j)=Gfront;

        //extension by duplication
        std::size_t IsquareXmin = (i < step/2 ? 0 : i-step/2);
        std::size_t IsquareXmax = std::min(M.DTM.width()-1,i+step/2);
        std::size_t JsquareYmin = (j < step/2 ? 0 : j-step/2);
        std::size_t JsquareYmax = std::min(M.DTM.height()-1,j+step/2);
		
        if(M.DTM.width()-1-IsquareXmax<step) IsquareXmax=M.DTM.width()-1;
        if(M.DTM.height()-1-JsquareYmax<step) JsquareYmax=M.DTM.height()-1;
        if(IsquareXmin<step) IsquareXmin=0;
        if(JsquareYmin<step) JsquareYmin=0;

        for(std::size_t k = IsquareXmin; k <= IsquareXmax; k++){
          for(std::size_t l = JsquareYmin; l <= JsquareYmax; l++){
            if(M.Mask(k,l)){
              M.DTM(k,l)=G1; 
              im_Zfront(k,l)=Gfront;
            }
          }
        }
      }
    }

    //Gaussian smoothness
    for(std::size_t ik = 0; ik < 10; ik++) {
      for(std::size_t j = 0; j < M.DTM.height(); j++){ 
        for (std::size_t i = 0; i < M.DTM.width(); i++) {
		
          std::size_t IsquareXmin = (i < 1 ? 0 : i - 1);
          std::size_t IsquareXmax = std::min(M.DTM.width()-1,i+1); 
          std::size_t JsquareYmin = (j < 1 ? 0 : j - 1);
          std::size_t JsquareYmax = std::min(M.DTM.height()-1,j+1); 
          int count=0; 

          for(std::size_t k=IsquareXmin; k<=IsquareXmax; k++){ 
            for(std::size_t l=JsquareYmin; l<=JsquareYmax; l++){ 
              if(M.Mask(k,l)) count++;
            }
          }	
          if(count==9) M.DTM(i,j)=(M.DTM(i-1,j-1)+2*M.DTM(i,j-1)+M.DTM(i+1,j-1)+2*M.DTM(i-1,j)+4*M.DTM(i,j)+2*M.DTM(i+1,j)+M.DTM(i-1,j+1)+2*M.DTM(i,j+1)+M.DTM(i+1,j+1))/16;
        }
      }
    }

    //ranger les valeurs dans scans lidars
    for(int i=0;i<(int)M.HPS.size();i++){
      int I=M.HPS[i].ind_x;
      int J=M.HPS[i].ind_y;
      elevation_attribute.push_back ((double)(M.HPS[i].position.z()-M.DTM(I,J)));
    }

    this->compute_mean_max (elevation_attribute, mean, max);
    //    max *= 5;
  }

  virtual double value (int pt_index)
  {
    return std::max (0., std::min (1., elevation_attribute[pt_index] / weight));
  }
  
  virtual std::string id() { return "elevation"; }

};
  


}

#endif // CGAL_DATA_CLASSIFICATION_SEGMENTATION_ATTRIBUTE_ELEVATION_H
