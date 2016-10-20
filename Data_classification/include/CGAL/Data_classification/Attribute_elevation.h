// Copyright (c) 2016  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Simon Giraudot, Florent Lafarge

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

    \brief Attribute based on local elevation.

    The local position of the ground can be computed for urban
    scenes. This attribute computes the distance to the local
    estimation of the ground.

    It is useful to discriminate the ground from horizontal roofs.

    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointMap Property map to access the input points.

  */
template <typename Kernel, typename RandomAccessIterator, typename PointMap>
class Attribute_elevation : public Attribute
{
  typedef typename Kernel::Iso_cuboid_3 Iso_cuboid_3;

  typedef Image<float> Image_float;
  typedef Planimetric_grid<Kernel, RandomAccessIterator, PointMap> Grid;
   
  std::vector<double> elevation_attribute;
  
public:
  /*!
    \brief Constructs the attribute.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_map Property map to access the input points
    \param bbox Bounding box of the input range
    \param grid Precomputed `Planimetric_grid`
    \param grid_resolution Resolution of the planimetric grid
    \param radius_neighbors Radius of local neighborhoods
    \param radius_dtm Radius for digital terrain modeling (must be bigger than the size of a building)

  */
  Attribute_elevation (RandomAccessIterator begin,
                       RandomAccessIterator end,
                       PointMap point_map,
                       const Iso_cuboid_3& bbox,
                       const Grid& grid,
                       const double grid_resolution,
                       double radius_neighbors = -1.,
                       double radius_dtm = -1.)
  {
    this->weight = 1.;
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
        list_Z.reserve (grid.indices(i,j).size());
        for(std::size_t k=0;k<grid.indices(i,j).size();k++) list_Z.push_back(get(point_map, begin[grid.indices(i,j)[k]]).z());
			
        if(list_Z.size()>1){
			
          std::sort(list_Z.begin(),list_Z.end());
          int ind_k= (int)floor((double)(list_Z.size()*0.9)-0.5); // 9/10
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
          std::size_t squareXmax = (std::min)(dem.width() - 1, i + square);
          std::size_t squareYmin = (j < square ? 0 : j - square);
          std::size_t squareYmax = (std::min)(dem.height() - 1, j + square);
			
          for(std::size_t k = squareXmin; k <= squareXmax; k++){
            for(std::size_t l = squareYmin; l <= squareYmax; l++){
					
              double distance = (CGAL::square((double)i-k)+CGAL::square((double)j-l))*CGAL::square(grid_resolution);
					
              if((distance<=CGAL::square(radius_neighbors))&&(dem(k,l)>0)&&(distance!=0)){
                val=val+dem(k,l)/distance;
                distance_tot=distance_tot+1/distance;
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

    std::size_t step=15;
    square = (std::size_t)(radius_dtm/grid_resolution)+1;
    
    Image_float dtm (grid.width(), grid.height());
    Image_float im_Zfront(grid.width(),grid.height());

    for (std::size_t j = step/2+1; j < dtm.height(); j = j+step){	
      for (std::size_t i = step+1/2; i < dtm.width(); i = i+step){
		
        //storage of the points in the disk
        std::vector<float> list_pointsZ;
        std::size_t squareXmin = (i < square ? 0 : i - square);
        std::size_t squareXmax = (std::min)(dtm.width() - 1, i + square);
        std::size_t squareYmin = (j < square ? 0 : j - square);
        std::size_t squareYmax = (std::min)(dtm.height() - 1, j + square);

        for(std::size_t k = squareXmin; k <= squareXmax; k++){
          for(std::size_t l = squareYmin; l <= squareYmax; l++){
			
            double distance= (CGAL::square((double)i-k)+CGAL::square((double)j-l))*CGAL::square(grid_resolution);
            if(distance <= CGAL::square(radius_dtm)){
              for(int nb=0;nb<(int)grid.indices(k,l).size();nb++)
                list_pointsZ.push_back(get(point_map, begin[grid.indices(k,l)[nb]]).z());
            }
          }
        }

        //ordering of the points
        float G1=0; 
        float G2=0;

        if(list_pointsZ.size()>1){
          std::sort(list_pointsZ.begin(),list_pointsZ.end());
          int ind_k2= (int)floor((double)(list_pointsZ.size()*0.6)-0.5); // 6/10
          int ind_k1= (int)floor((double)(list_pointsZ.size()*0.3)-0.5); // 3/10
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
        std::size_t IsquareXmax = (std::min)(dtm.width()-1,i+step/2);
        std::size_t JsquareYmin = (j < step/2 ? 0 : j-step/2);
        std::size_t JsquareYmax = (std::min)(dtm.height()-1,j+step/2);
		
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
          std::size_t IsquareXmax = (std::min)(dtm.width()-1,i+1); 
          std::size_t JsquareYmin = (j < 1 ? 0 : j-1);
          std::size_t JsquareYmax = (std::min)(dtm.height()-1,j+1); 
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
    test_ground.reserve(end - begin);
    for (std::size_t i = 0; i< (std::size_t)(end - begin); i++){
      int I = grid.x(i);
      int J = grid.y(i);
      bool test=false;
      if(get(point_map, begin[i]).z()-im_Zfront(I,J)<=0) {test=true; test_ground.push_back(test);}
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
        std::size_t squareXmax = (std::min)(dtm.width()-1,i+square);
        std::size_t squareYmin = (j < square ? 0 : j-square);
        std::size_t squareYmax = (std::min)(dtm.height()-1,j+square);

        for(std::size_t k = squareXmin; k <= squareXmax; k++){
          for(std::size_t l = squareYmin; l <= squareYmax; l++){
				
            double distance = (CGAL::square((double)i-k)+CGAL::square((double)j-l))*CGAL::square(grid_resolution);
				
            if(distance <= CGAL::square(radius_dtm)){
					
              for(int nb=0;nb<(int)grid.indices(k,l).size();nb++){
                if(test_ground[grid.indices(k,l)[nb]])
                  list_pointsZ.push_back(get(point_map, begin[grid.indices(k,l)[nb]]).z());
              }
            }
          }
        }

        //ordering the new points
        float G1=0; 
        float G2=0;

        if(list_pointsZ.size()>1){
          std::sort(list_pointsZ.begin(),list_pointsZ.end());
          int ind_k2= (int)floor((double)(list_pointsZ.size()*0.6)-0.5); // 6/10
          int ind_k1= (int)floor((double)(list_pointsZ.size()*0.3)-0.5); // 3/10
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
        std::size_t IsquareXmax = (std::min)(dtm.width()-1,i+step/2);
        std::size_t JsquareYmin = (j < step/2 ? 0 : j-step/2);
        std::size_t JsquareYmax = (std::min)(dtm.height()-1,j+step/2);
		
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
          std::size_t IsquareXmax = (std::min)(dtm.width()-1,i+1); 
          std::size_t JsquareYmin = (j < 1 ? 0 : j - 1);
          std::size_t JsquareYmax = (std::min)(dtm.height()-1,j+1); 
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
    elevation_attribute.reserve(end - begin);
    for (std::size_t i = 0; i < (std::size_t)(end - begin); i++){
      int I = grid.x(i);
      int J = grid.y(i);
      elevation_attribute.push_back ((double)(get(point_map, begin[i]).z()-dtm(I,J)));
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
