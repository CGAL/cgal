#ifndef CGAL_DATA_CLASSIFICATION_ATTRIBUTE_ELEVATION_H
#define CGAL_DATA_CLASSIFICATION_ATTRIBUTE_ELEVATION_H

#include <vector>

#include <CGAL/Data_classification/Attribute.h>
#include <CGAL/Data_classification/Image.h>
#include <CGAL/Data_classification/Planimetric_grid.h>

namespace CGAL {

namespace Data_classification {

  /*!
    \ingroup PkgDataClassification

    \brief Segmentation attribute based on local elevation.

    The local position of the ground can be computed for urban
    scenes. This attribute computes the distance to the local
    estimation of the ground.

    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointPMap Property map to access the input points.

  */
template <typename Kernel, typename RandomAccessIterator, typename PointPMap>
class Attribute_elevation : public Attribute
{
  typedef typename Kernel::Iso_cuboid_3 Iso_cuboid_3;

  typedef Image<float> Image_float;
  typedef Planimetric_grid<Kernel, RandomAccessIterator, PointPMap> Grid;
   
  std::vector<double> elevation_attribute;
  
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_pmap Property map to access the input points
    \param bbox Bounding box of the input range
    \param grid Precomputed `Planimetric_grid`
    \param grid_resolution Resolution of the planimetric grid
    \param radius_neighbors Radius of local neighborhoods
    \param radius_dtm Radius for digital terrain modeling (must be bigger than the size of a building)
    \param weight Weight of the attribute

  */
  Attribute_elevation (RandomAccessIterator begin,
                       RandomAccessIterator end,
                       PointPMap point_pmap,
                       const Iso_cuboid_3& bbox,
                       const Grid& grid,
                       const double grid_resolution,
                       double radius_neighbors = -1.,
                       double radius_dtm = -1.,
                       double weight = 1.)
  {
    this->weight = weight;
    if (radius_neighbors < 0.)
      radius_neighbors = 5. * grid_resolution;
    if (radius_dtm < 0.)
      radius_dtm = 5. * radius_neighbors;

    //DEM
    Image_float dem(grid.width(),grid.height());
    Image_float demtemp(grid.width(),grid.height());

    for (int j=0;j<(int)dem.height();j++){	
      for (int i=0;i<(int)dem.width();i++){
		
        double mean_height=0;
        std::vector < double > list_Z;
		
        for(std::size_t k=0;k<grid.indices(i,j).size();k++) list_Z.push_back(get(point_pmap, begin[grid.indices(i,j)[k]]).z());
			
        if(list_Z.size()>0){
			
          std::sort(list_Z.begin(),list_Z.end());
          int ind_k= (int)floor((double)(list_Z.size()*0.9)-0.5);
          mean_height=list_Z[ind_k];
        }

        dem(i,j)=(float)mean_height;
        demtemp(i,j)=dem(i,j);
      }
    }

    std::size_t square = (std::size_t)(radius_neighbors / grid_resolution) + 1;

    for (std::size_t j = 0; j < dem.height(); j++){	
      for (std::size_t i = 0; i < dem.width(); i++){
		
        if((grid.mask(i,j))&&(dem(i,j)==0)){
				
          double distance_tot=0;	
          double val=0;
          std::size_t squareXmin = (i < square ? 0 : i - square);
          std::size_t squareXmax = std::min(dem.width() - 1, i + square);
          std::size_t squareYmin = (j < square ? 0 : j - square);
          std::size_t squareYmax = std::min(dem.height() - 1, j + square);
			
          for(std::size_t k = squareXmin; k <= squareXmax; k++){
            for(std::size_t l = squareYmin; l <= squareYmax; l++){
					
              double distance=sqrt(pow((double)i-k,2)+pow((double)j-l,2))*grid_resolution;
					
              if((distance<=radius_neighbors)&&(dem(k,l)>0)&&(distance!=0)){
                double dista=distance*distance;
                val=val+dem(k,l)/dista;
                distance_tot=distance_tot+1/dista;
              }
            }
          }
			
          val=val/distance_tot;
          demtemp(i,j)=(float)val;			
        }
      }
    }

    float scale_byte_min_dem = bbox.zmax();
    float scale_byte_max_dem = bbox.zmin();

    for (std::size_t j = 0; j < dem.height(); j++){	
      for (std::size_t i = 0; i < dem.width(); i++){
		
        dem(i,j)=demtemp(i,j);

        if(grid.mask(i,j)){
          if(dem(i,j)>scale_byte_max_dem) scale_byte_max_dem=dem(i,j);
          if(dem(i,j)<scale_byte_min_dem) scale_byte_min_dem=dem(i,j);
        }
      }
    }

    //TODO: smooth the dem
    for (std::size_t j = 1; j < dem.height()-1; j++){	
      for (std::size_t i = 1; i < dem.width()-1; i++){
        if(grid.mask(i,j)){
          //dem(i,j)=
        }
      }
    }

    //DTM computation
    std::size_t step=15;
    square = (std::size_t)(radius_dtm/grid_resolution)+1;
    
    Image_float dtm (grid.width(), grid.height());
    Image_float im_Zfront(grid.width(),grid.height());


    //round 1
    for (std::size_t j = 0; j < dtm.height(); j++){	
      for (std::size_t i = 0; i < dtm.width(); i++){
        dtm(i,j)=0;
        im_Zfront(i,j)=0;
      }
    }


    for (std::size_t j = step/2+1; j < dtm.height(); j = j+step){	
      for (std::size_t i = step+1/2; i < dtm.width(); i = i+step){
		
        //storage of the points in the disk
        std::vector<float> list_pointsZ;
        std::size_t squareXmin = (i < square ? 0 : i - square);
        std::size_t squareXmax = std::min(dtm.width() - 1, i + square);
        std::size_t squareYmin = (j < square ? 0 : j - square);
        std::size_t squareYmax = std::min(dtm.height() - 1, j + square);

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

        dtm(i,j)=G1;
        im_Zfront(i,j)=Gfront;

        //extension by duplication
        std::size_t IsquareXmin = (i < step/2 ? 0 : i-step/2);
        std::size_t IsquareXmax = std::min(dtm.width()-1,i+step/2);
        std::size_t JsquareYmin = (j < step/2 ? 0 : j-step/2);
        std::size_t JsquareYmax = std::min(dtm.height()-1,j+step/2);
		
        if(dtm.width()-1-IsquareXmax<step) IsquareXmax=dtm.width()-1;
        if(dtm.height()-1-JsquareYmax<step) JsquareYmax=dtm.height()-1;
        if(IsquareXmin<step) IsquareXmin=0;
        if(JsquareYmin<step) JsquareYmin=0;

        for(std::size_t k = IsquareXmin; k <= IsquareXmax; k++){
          for(std::size_t l = JsquareYmin; l <= JsquareYmax; l++){
            if(grid.mask(k,l)){
              dtm(k,l)=G1;
              im_Zfront(k,l)=Gfront;
            }
          }
        }
      }
    }

    //Gaussian smoothness
    for(std::size_t ik = 0; ik < 10; ik++) {
      for(std::size_t j = 0; j < dtm.height(); j++){ 
        for (std::size_t i = 0; i < dtm.width(); i++) {
		
          std::size_t IsquareXmin = (i < 1 ? 0 : i-1);
          std::size_t IsquareXmax = std::min(dtm.width()-1,i+1); 
          std::size_t JsquareYmin = (j < 1 ? 0 : j-1);
          std::size_t JsquareYmax = std::min(dtm.height()-1,j+1); 
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
    for (std::size_t i = 0; i< (std::size_t)(end - begin); i++){
      int I = grid.x(i);
      int J = grid.y(i);
      bool test=false;
      if(get(point_pmap, begin[i]).z()-im_Zfront(I,J)<=0) {test=true; test_ground.push_back(test);}
      else{test_ground.push_back(test);}
    }

    for (int j=0;j<(int)dtm.height();j++){	
      for (int i=0;i<(int)dtm.width();i++){
        dtm(i,j)=0;
        im_Zfront(i,j)=0;
      }
    }

    for (std::size_t j = step/2+1; j < dtm.height(); j = j+step){	
      for (std::size_t i = step/2+1; i < dtm.width(); i = i+step){
		
        //stockage des nouveaux points
        std::vector < float > list_pointsZ;
        std::size_t squareXmin = (i < square ? 0 : i-square);
        std::size_t squareXmax = std::min(dtm.width()-1,i+square);
        std::size_t squareYmin = (j < square ? 0 : j-square);
        std::size_t squareYmax = std::min(dtm.height()-1,j+square);

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
        dtm(i,j)=G1;
        im_Zfront(i,j)=Gfront;

        //extension by duplication
        std::size_t IsquareXmin = (i < step/2 ? 0 : i-step/2);
        std::size_t IsquareXmax = std::min(dtm.width()-1,i+step/2);
        std::size_t JsquareYmin = (j < step/2 ? 0 : j-step/2);
        std::size_t JsquareYmax = std::min(dtm.height()-1,j+step/2);
		
        if(dtm.width()-1-IsquareXmax<step) IsquareXmax=dtm.width()-1;
        if(dtm.height()-1-JsquareYmax<step) JsquareYmax=dtm.height()-1;
        if(IsquareXmin<step) IsquareXmin=0;
        if(JsquareYmin<step) JsquareYmin=0;

        for(std::size_t k = IsquareXmin; k <= IsquareXmax; k++){
          for(std::size_t l = JsquareYmin; l <= JsquareYmax; l++){
            if(grid.mask(k,l)){
              dtm(k,l)=G1; 
              im_Zfront(k,l)=Gfront;
            }
          }
        }
      }
    }

    //Gaussian smoothness
    for(std::size_t ik = 0; ik < 10; ik++) {
      for(std::size_t j = 0; j < dtm.height(); j++){ 
        for (std::size_t i = 0; i < dtm.width(); i++) {
		
          std::size_t IsquareXmin = (i < 1 ? 0 : i - 1);
          std::size_t IsquareXmax = std::min(dtm.width()-1,i+1); 
          std::size_t JsquareYmin = (j < 1 ? 0 : j - 1);
          std::size_t JsquareYmax = std::min(dtm.height()-1,j+1); 
          int count=0; 

          for(std::size_t k=IsquareXmin; k<=IsquareXmax; k++){ 
            for(std::size_t l=JsquareYmin; l<=JsquareYmax; l++){ 
              if(grid.mask(k,l)) count++;
            }
          }	
          if(count==9) dtm(i,j)=(dtm(i-1,j-1)+2*dtm(i,j-1)+dtm(i+1,j-1)+2*dtm(i-1,j)+4*dtm(i,j)+2*dtm(i+1,j)+dtm(i-1,j+1)+2*dtm(i,j+1)+dtm(i+1,j+1))/16;
        }
      }
    }

    //ranger les valeurs dans scans lidars
    for (std::size_t i = 0; i < (std::size_t)(end - begin); i++){
      int I = grid.x(i);
      int J = grid.y(i);
      elevation_attribute.push_back ((double)(get(point_pmap, begin[i]).z()-dtm(I,J)));
    }

    this->compute_mean_max (elevation_attribute, this->mean, this->max);
    //    max *= 5;
  }

  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return elevation_attribute[pt_index];
  }
  
  virtual std::string id() { return "elevation"; }
  /// \endcond
};

} // namespace Data_classification


} // namespace CGAL

#endif // CGAL_DATA_CLASSIFICATION_ATTRIBUTE_ELEVATION_H
