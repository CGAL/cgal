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
// Author(s)     : Florent Lafarge, Simon Giraudot

#ifndef CGAL_POINT_SET_CLASSIFICATION_H
#define CGAL_POINT_SET_CLASSIFICATION_H

#include <cstdio>
#include <cassert>
#include <vector>
#include <list>
#include <set>
#include <string>

#include <CGAL/bounding_box.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Default_diagonalize_traits.h>
#include <CGAL/centroid.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#define CGAL_CLASSIFICATION_VERBOSE
#if defined(CGAL_CLASSIFICATION_VERBOSE)
#define CGAL_CLASSIFICATION_CERR std::cerr
#else
#define CGAL_CLASSIFICATION_CERR std::ostream(0)
#endif


namespace CGAL {

namespace internal {

namespace Classification {


template <typename Type>
class Image
{
  std::size_t m_width;
  std::size_t m_height;
  Type* m_raw;


public:

  Image () : m_width(0), m_height(0), m_raw (NULL)
  {
  }
  
  Image (std::size_t width, std::size_t height)
    : m_width (width),
      m_height (height)
  {
    if (m_width * m_height > 0)
      m_raw = new Type[width * height];
    else
      m_raw = NULL;
  }
  
  ~Image ()
  {
    if (m_raw != NULL)
      delete[] m_raw;
  }

  Image (const Image& other)
    : m_width (other.width()),
      m_height (other.height())

  {
    if (m_width * m_height > 0)
      {
        m_raw = new Type[width * height];
        std::copy (other.m_raw, other.m_raw + (m_width * m_height), this->m_raw);
      }
    else
      m_raw = NULL;
  }
  Image& operator= (const Image& other)
  {
    if (m_raw != NULL)
      delete m_raw;
    m_raw = NULL;
    m_width = other.width();
    m_height = other.height();
    if (m_width * m_height > 0)
      {
        m_raw = new Type[m_width * m_height];
        std::copy (other.m_raw, other.m_raw + (m_width * m_height), this->m_raw);
      }
    
    return *this;
  }
  
  std::size_t width() const { return m_width; }
  std::size_t height() const { return m_height; }

  Type& operator() (const std::size_t& x, const std::size_t& y)
  {
    return m_raw[y * m_width + x];
  }
  const Type& operator() (const std::size_t& x, const std::size_t& y) const
  {
    return m_raw[y * m_width + x];
  }
  

};

  
inline void compute_mean_max (std::vector<double>& vect, double& mean, double& max)
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


template <typename HSV_Color, typename RGB_Color>
HSV_Color rgb_to_hsv (const RGB_Color& c)
{
  double r = (double)(c[0]) / 255.;
  double g = (double)(c[1]) / 255.;
  double b = (double)(c[2]) / 255.;
  double Cmax = std::max (r, std::max (g, b));
  double Cmin = std::min (r, std::min (g, b));
  double delta = Cmax - Cmin;
  double H = 0.;
  
  if (delta != 0.)
    {
      if (Cmax == r)
        H = 60. * (((int)((g - b) / delta)) % 6);
      else if (Cmax == g)
        H = 60. * (((b - r) / delta) + 2.);
      else
        H = 60. * (((r - g) / delta) + 4.);
    }
  if (H < 0.) H += 360.;
  double S = (Cmax == 0. ? 0. : 100. * (delta / Cmax));
  double V = 100. * Cmax;
  HSV_Color out = {{ H, S, V }};
  return out;
}

} // namespace internal

} // namespace Classification

template <typename Kernel>
class Point_set_classification;

class Abstract_classification_type
{
public:

  virtual double data_term_computation (int pt_index) = 0;

};


template <typename Kernel>
class Point_set_classification
{
public:
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Segment_3 Segment;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Line_3 Line;
  typedef typename Kernel::Plane_3 Plane;
  typedef typename Kernel::Vector_3 Vector;

  typedef CGAL::cpp11::array<double, 3> Eigenvalues;

  typedef CGAL::Search_traits_3<Kernel> SearchTraits_3;
  typedef CGAL::Kd_tree<SearchTraits_3> Tree;
  typedef CGAL::Fuzzy_sphere<SearchTraits_3> Fuzzy_sphere;
  typedef typename Kernel::Iso_cuboid_3 Iso_cuboid_3;

  typedef CGAL::cpp11::array<unsigned char, 3> Color;
  typedef CGAL::cpp11::array<double, 3> HSV_Color;

  typedef internal::Classification::Image<std::vector<int> > Image_indices;
  typedef internal::Classification::Image<bool> Image_bool;
  typedef internal::Classification::Image<float> Image_float;
  
  struct HPoint {
    Point position;
    unsigned char echo; 
    int ind_x;
    int ind_y;
    int ultimate_label; 
    unsigned char AE_label;
    Color color;
  };

  std::vector<Abstract_classification_type*> segmentation_classes;
  
  bool is_echo_given;
  bool is_normal_given;
  Iso_cuboid_3 BBox_scan;

  std::vector<HPoint> HPS;
  Image_indices grid_HPS;

  //Hpoints attributes
  std::vector<Vector> normal;
  std::vector<Eigenvalues> eigenvalues;
  std::vector<Plane> planes;
  std::vector<std::vector<int> > spherical_neighborhood;
  std::vector<std::vector<int> > cylindrical_neighborhood;

  Image_bool Mask;                     //imagg
  Image_float DTM; //a enregistrer ?   //imagg           

  double m_grid_resolution;
  double m_radius_neighbors; 
  double m_radius_dtm; 

  template <typename InputIterator>
  Point_set_classification (InputIterator begin, InputIterator end,
                            double grid_resolution = -1.,
                            double radius_neighbors = -1.,
                            double radius_dtm = -1.)
    : m_grid_resolution (grid_resolution),
      m_radius_neighbors (radius_neighbors),
      m_radius_dtm (radius_dtm)
  {
    if (m_grid_resolution < 0.)
      m_grid_resolution = CGAL::compute_average_spacing<CGAL::Sequential_tag> (begin, end, 6);
    if (m_radius_neighbors < 0.)
      m_radius_neighbors = 5. * m_grid_resolution;
    if (m_radius_dtm < 0.)
      m_radius_dtm = 5 * m_grid_resolution;

    for (InputIterator it = begin; it != end; ++ it)
      {
        HPS.push_back (HPoint());
        HPS.back().position = *it;
        HPS.back().echo = -1;
        HPS.back().ind_x = -1;
        HPS.back().ind_y = -1;
        HPS.back().ultimate_label = -1;
        HPS.back().AE_label = -1;
        Color c = {{ 0, 0, 0 }};
        HPS.back().color = c;
      }
  }


  void set_parameters (double grid_resolution, double radius_neighbors, double radius_dtm)
  {
    m_grid_resolution = grid_resolution;
    m_radius_neighbors = radius_neighbors;
    m_radius_dtm = radius_dtm;
  }
  
  bool initialization(int phase = 0) // for training
  {
    clock_t t;
    t = clock();

    CGAL_CLASSIFICATION_CERR << std::endl << "Initialization: ";

    //1-Neighborhood computation and reset the attributes of the structure points
    CGAL_CLASSIFICATION_CERR<<"spherical neighborhood..";
    spherical_neighborhood.clear();
    eigenvalues.clear();
    planes.clear();
    
    std::map<Point,int> map_indice_point;
    for(std::size_t i=0;i<HPS.size();i++){
      Point pt=HPS[i].position;
      map_indice_point[pt]=i;
	
    }

    std::vector<Point> list_points;
    for(int i=0;i<(int)HPS.size();i++){
      Point pt=HPS[i].position;
      list_points.push_back(pt);
    }

    Tree tree(list_points.begin(), list_points.end());

    std::size_t nb_neigh = 0;
    for(int i=0;i<(int)HPS.size();i++){
      const Point& query=HPS[i].position;
      std::vector<Point> neighbors;
      
      if (phase == 1 && HPS[i].echo == 0)
        {
          compute_principal_curvature (query, neighbors);
          continue;
        }
      
      Fuzzy_sphere fs(query, m_radius_neighbors, 0);
      tree.search(std::back_inserter(neighbors), fs);
      nb_neigh += neighbors.size();

      compute_principal_curvature (query, neighbors);
    }

    CGAL_CLASSIFICATION_CERR<<"ok";

    if (phase == 1)
      return true;

    //2-creation of the bounding box
    CGAL_CLASSIFICATION_CERR<<", bounding box..";
    BBox_scan = CGAL::bounding_box(list_points.begin(), list_points.end());
    CGAL_CLASSIFICATION_CERR<<"ok";

    //3-creation grille_points
    CGAL_CLASSIFICATION_CERR<<", planimetric grid of HPS..";
    Image_indices tess((std::size_t)((BBox_scan.xmax()-BBox_scan.xmin())/m_grid_resolution)+1,
                       (std::size_t)((BBox_scan.ymax()-BBox_scan.ymin())/m_grid_resolution)+1);
    grid_HPS=tess;

    for(int i=0;i<(int)HPS.size();i++){

      //for each 3D point, its coordinates in the grid are inserted in its Hpoint structure
      HPS[i].ind_x=(int)((HPS[i].position.x()-BBox_scan.xmin())/m_grid_resolution);
      HPS[i].ind_y=(int)((HPS[i].position.y()-BBox_scan.ymin())/m_grid_resolution);

      //index of points are collected in grid_HPS
      std::vector < int > temp;
      temp=grid_HPS(HPS[i].ind_x,HPS[i].ind_y);
      temp.push_back(i);
      grid_HPS(HPS[i].ind_x,HPS[i].ind_y)=temp;
    }
    CGAL_CLASSIFICATION_CERR<<"ok";


    //4-Mask creation
    CGAL_CLASSIFICATION_CERR<<", planimetric mask..";
    Image_bool masktp((std::size_t)((BBox_scan.xmax()-BBox_scan.xmin())/m_grid_resolution)+1,
                      (std::size_t)((BBox_scan.ymax()-BBox_scan.ymin())/m_grid_resolution)+1);
    Mask=masktp;

    CGAL_CLASSIFICATION_CERR << "(" << Mask.height() << "x" << Mask.width() << ")" << std::endl;
    int square=(int)16*(m_radius_neighbors/m_grid_resolution)+1;
    int nb_true = 0;
    for (int j=0;j<(int)Mask.height();j++){	
      for (int i=0;i<(int)Mask.width();i++){	
        //Mask(i,j)=false;
        if(grid_HPS(i,j).size()>0) { Mask(i,j)=true; nb_true ++; }
        else{Mask(i,j)=false;}
      }
    }

    Image_bool Mask_tmp ((std::size_t)((BBox_scan.xmax()-BBox_scan.xmin())/m_grid_resolution)+1,
                         (std::size_t)((BBox_scan.ymax()-BBox_scan.ymin())/m_grid_resolution)+1);

    for (std::size_t i = 0; i < Mask.width(); ++ i)
      for (std::size_t j = 0; j < Mask.height(); ++ j)
        {
          if(!Mask(i,j))
            {
              int squareYmin=std::max(0,(int)j-square);
              int squareYmax=std::min((int)(Mask.height())-1,(int)j+square);

              for (int k = squareYmin; k <= squareYmax; ++ k)
                if (Mask(i,k))
                  {
                    Mask_tmp(i,j) = true;
                    break;
                  }
            }
          else
            Mask_tmp(i,j) = true;
        }

    for (std::size_t i = 0; i < Mask.width(); ++ i)
      for (std::size_t j = 0; j < Mask.height(); ++ j)
        {
          if(!Mask_tmp(i,j))
            {
              int squareXmin=std::max(0,(int)i-square);
              int squareXmax=std::min((int)(Mask.width())-1,(int)i+square);

              for (int k = squareXmin; k <= squareXmax; ++ k)
                if (Mask_tmp(k,j))
                  {
                    Mask(i,j) = true;
                    break;
                  }
            }
          else
            Mask(i,j) = true;
        }




    //5-normal computation
    //    if(!is_normal_given){CGAL_CLASSIFICATION_CERR<<", normals.."; compute_normal(); CGAL_CLASSIFICATION_CERR<<"ok";}

    CGAL_CLASSIFICATION_CERR<<std::endl<<"-> OK ( "<<((float)clock()-t)/CLOCKS_PER_SEC<<" sec )"<< std::endl;

    return true;
  }





  void compute_principal_curvature(const Point& point, std::vector<Point>& neighborhood)
  {
    if (neighborhood.size() == 0)
      {
        Eigenvalues v = {{ 0., 0., 0. }};
        eigenvalues.push_back (v);
        planes.push_back (Plane (point, Vector (0., 0., 1.)));
        return;
      }
    Point centroid = CGAL::centroid (neighborhood.begin(), neighborhood.end());

    CGAL::cpp11::array<double, 6> covariance = {{ 0., 0., 0., 0., 0., 0. }};
      
    for (std::size_t i = 0; i < neighborhood.size(); ++ i)
      {
        Vector d = neighborhood[i] - centroid;
        covariance[0] += d.x () * d.x ();
        covariance[1] += d.x () * d.y ();
        covariance[2] += d.x () * d.z ();
        covariance[3] += d.y () * d.y ();
        covariance[4] += d.y () * d.z ();
        covariance[5] += d.z () * d.z ();
      }

    Eigenvalues evalues = {{ 0., 0., 0. }};
    CGAL::cpp11::array<double, 9> eigenvectors = {{ 0., 0., 0.,
                                                    0., 0., 0.,
                                                    0., 0., 0. }};

    CGAL::Default_diagonalize_traits<double,3>::diagonalize_selfadjoint_covariance_matrix
      (covariance, evalues, eigenvectors);
    eigenvalues.push_back (evalues);

    Plane plane;
    CGAL::linear_least_squares_fitting_3 (neighborhood.begin(),neighborhood.end(),plane, CGAL::Dimension_tag<0>());
    planes.push_back (plane);
  }



  bool quick_labeling_PC()
  {

    // data term initialisation
    CGAL_CLASSIFICATION_CERR << "Labeling... ";


    int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
    for(int s=0;s<(int)HPS.size();s++){
			
      int nb_class_best=0; 

      double val_class_best = std::numeric_limits<double>::max();

      for(std::size_t k = 0; k < segmentation_classes.size(); ++ k)
        {
          double value = segmentation_classes[k]->data_term_computation(s);
          
          if(val_class_best > value)
            {
              val_class_best = value;
              nb_class_best=k;
            }
        }

      HPS[s].AE_label = nb_class_best;

      if(nb_class_best==0) count1++;
      else if(nb_class_best==1) count2++;
      else if(nb_class_best==2) count3++;
      else count4++;

    }
    
    CGAL_CLASSIFICATION_CERR<<" "<<(double)100*count1/HPS.size()<<"% vegetation, "<<(double)100*count2/HPS.size()<<"% ground, "<<(double)100*count3/HPS.size()<<"% roof, "<<(double)100*count4/HPS.size()<<"% clutter"<<std::endl;
	
    return true;
  }


  bool smoothed_labeling_PC()
  {

    // data term initialisation
    CGAL_CLASSIFICATION_CERR << "Labeling... ";

    std::vector<std::vector<double> > values;
    values.resize (segmentation_classes.size());
    
    std::map<Point, std::size_t> map_p2i;
    std::vector<Point> list_points;    
    for(int s=0;s<(int)HPS.size();s++)
      {
        map_p2i[HPS[s].position] = s;
        list_points.push_back(HPS[s].position);
        for(std::size_t k = 0; k < segmentation_classes.size(); ++ k)
          values[k].push_back(segmentation_classes[k]->data_term_computation(s));
      }

    Tree tree(list_points.begin(), list_points.end());

    std::vector<std::vector<double> > smoothed_values;
    smoothed_values.resize (segmentation_classes.size());

    for(int s=0;s<(int)HPS.size();s++)
      {
        const Point& query=HPS[s].position;
        std::vector<Point> neighbors;
      
        Fuzzy_sphere fs(query, m_radius_neighbors, 0);
        tree.search(std::back_inserter(neighbors), fs);

        std::vector<double> mean (values.size(), 0.);
        for (std::size_t n = 0; n < neighbors.size(); ++ n)
          {
            std::size_t index = map_p2i[neighbors[n]];
            for (std::size_t j = 0; j < values.size(); ++ j)
              mean[j] += values[j][index];
          }
        for (std::size_t j = 0; j < values.size(); ++ j)
          smoothed_values[j].push_back (mean[j] / neighbors.size());
      }
    
    int count1 = 0, count2 = 0, count3 = 0, count4 = 0;    
    for(int s=0;s<(int)HPS.size();s++){
			
      int nb_class_best=0; 

      double val_class_best = std::numeric_limits<double>::max();

      for(std::size_t k = 0; k < smoothed_values.size(); ++ k)
        {
          if(val_class_best > smoothed_values[k][s])
            {
              val_class_best = smoothed_values[k][s];
              nb_class_best=k;
            }
        }

      HPS[s].AE_label = nb_class_best;

      if(nb_class_best==0) count1++;
      else if(nb_class_best==1) count2++;
      else if(nb_class_best==2) count3++;
      else count4++;
    }

    CGAL_CLASSIFICATION_CERR<<" "<<(double)100*count1/HPS.size()<<"% vegetation, "<<(double)100*count2/HPS.size()<<"% ground, "<<(double)100*count3/HPS.size()<<"% roof, "<<(double)100*count4/HPS.size()<<"% clutter"<<std::endl;
    
	
    return true;
  }

  bool point_cloud_classification(int method = 1)
  {


    clock_t t;
    t = clock();
	
    CGAL_CLASSIFICATION_CERR<<std::endl<<"Classification of the point cloud: ";
    
    if (method == 1)
      quick_labeling_PC();
    else if (method == 2)
      smoothed_labeling_PC();
    else
      {
        std::cerr << "Unknown method number." << std::endl;
        abort();
      }

    CGAL_CLASSIFICATION_CERR<<"-> OK ( "<<((float)clock()-t)/CLOCKS_PER_SEC<<" sec )"<< std::endl;

    return true;
  }



protected: 



};


template <typename Kernel>
class Scatter_segmentation_attribute
{
  typedef Point_set_classification<Kernel> PSC;
  typedef typename PSC::Image_float Image_float;
  typedef typename PSC::Color Color;
  typedef typename PSC::HSV_Color HSV_Color;
  
  std::vector<double> vegetation_attribute;
  
public:
  double weight;
  double mean;
  double max;

  Scatter_segmentation_attribute (PSC& M,
                                  double weight,
                                  std::vector<Color> colors = std::vector<Color>()) : weight (weight)
  {
    Image_float Vegetation(M.grid_HPS.width(), M.grid_HPS.height());

    for (int j=0;j<(int)M.DTM.height();j++)	
      for (int i=0;i<(int)M.DTM.width();i++)
        Vegetation(i,j)=0;
		
    std::size_t square = (std::size_t)(0.5 * M.m_radius_neighbors / M.m_grid_resolution) + 1;

    if(M.is_echo_given){
				
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
    else if (colors.empty()){
      for(int i=0;i<(int)M.HPS.size();i++){
        if (M.eigenvalues[i][2] == 0)
          vegetation_attribute.push_back(0);
        else
          vegetation_attribute.push_back(std::fabs(M.eigenvalues[i][0] / M.eigenvalues[i][2]));
      }
    }
    else
      {
        std::cerr << "Using colors" << std::endl;
      for(int i=0;i<(int)M.HPS.size();i++)
        {
          HSV_Color c = internal::Classification::rgb_to_hsv<HSV_Color> (colors[i]);
          vegetation_attribute.push_back (std::exp (-(c[0] - 120) * (c[0] - 120) / (2 * 60 * 60))
                                          * std::exp (-(c[1] - 25) * (c[1] - 25) / (2 * 20 * 20))
                                          * std::exp (-(c[2] - 37) * (c[2] - 37) / (2 * 10 * 10)));
        }
      }

    internal::Classification::compute_mean_max (vegetation_attribute, mean, max);
    max *= 2;
  }

  double value (int pt_index)
  {
    return std::max (0., std::min (1., vegetation_attribute[pt_index] / weight));
  }

  const char* id() { return "scatter"; }
};


template <typename Kernel>
class Distance_to_plane_segmentation_attribute
{
  typedef Point_set_classification<Kernel> PSC;

  std::vector<double> distance_to_plane_attribute;
  
public:
  double weight;
  double mean;
  double max;
  
  Distance_to_plane_segmentation_attribute (PSC& M, double weight) : weight (weight)
  {
    for(int i=0; i<(int)M.HPS.size(); i++)
      distance_to_plane_attribute.push_back (std::sqrt (CGAL::squared_distance (M.HPS[i].position, M.planes[i])));

    internal::Classification::compute_mean_max (distance_to_plane_attribute, mean, max);
    max *= 2;
  }

  double value (int pt_index)
  {
    return std::max (0., std::min (1., distance_to_plane_attribute[pt_index] / weight));
  }

  const char* id() { return "planimetry"; }
};

template <typename Kernel>
class Elevation_segmentation_attribute
{
  typedef Point_set_classification<Kernel> PSC;
  typedef typename PSC::Image_float Image_float;
  
  std::vector<double> elevation_attribute;
  
public:
  double weight;
  double mean;
  double max;
  
  Elevation_segmentation_attribute (PSC& M, double weight) : weight (weight)
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

    internal::Classification::compute_mean_max (elevation_attribute, mean, max);
    max *= 5;
  }

  double value (int pt_index)
  {
    return std::max (0., std::min (1., elevation_attribute[pt_index] / weight));
  }

  const char* id() { return "elevetation"; }

};


template <typename Kernel>
class Vegetation_classification_type : public Abstract_classification_type
{
  Scatter_segmentation_attribute<Kernel>& scat_att;
  Distance_to_plane_segmentation_attribute<Kernel>& d2p_att;
  Elevation_segmentation_attribute<Kernel>& elev_att;
public:

  Vegetation_classification_type (Scatter_segmentation_attribute<Kernel>& scat_att,
                                  Distance_to_plane_segmentation_attribute<Kernel>& d2p_att,
                                  Elevation_segmentation_attribute<Kernel>& elev_att)
    : scat_att (scat_att)
    , d2p_att (d2p_att)
    , elev_att (elev_att)
  { }
  
  virtual double data_term_computation (int pt_index)
  {
    return (1. - scat_att.value (pt_index))
      + (1. - d2p_att.value (pt_index))
      + (1. - elev_att.value (pt_index));
  }

};

template <typename Kernel>
class Ground_classification_type : public Abstract_classification_type
{
  Scatter_segmentation_attribute<Kernel>& scat_att;
  Distance_to_plane_segmentation_attribute<Kernel>& d2p_att;
  Elevation_segmentation_attribute<Kernel>& elev_att;
public:

  Ground_classification_type (Scatter_segmentation_attribute<Kernel>& scat_att,
                              Distance_to_plane_segmentation_attribute<Kernel>& d2p_att,
                              Elevation_segmentation_attribute<Kernel>& elev_att)
    : scat_att (scat_att)
    , d2p_att (d2p_att)
    , elev_att (elev_att)
  { }

  virtual double data_term_computation (int pt_index)
  {
    return scat_att.value (pt_index)
      + d2p_att.value (pt_index)
      + elev_att.value (pt_index);
  }
};

template <typename Kernel>
class Roof_classification_type : public Abstract_classification_type
{
  Scatter_segmentation_attribute<Kernel>& scat_att;
  Distance_to_plane_segmentation_attribute<Kernel>& d2p_att;
  Elevation_segmentation_attribute<Kernel>& elev_att;
public:
  Roof_classification_type (Scatter_segmentation_attribute<Kernel>& scat_att,
                           Distance_to_plane_segmentation_attribute<Kernel>& d2p_att,
                           Elevation_segmentation_attribute<Kernel>& elev_att)
    : scat_att (scat_att)
    , d2p_att (d2p_att)
    , elev_att (elev_att)
  { }

  virtual double data_term_computation (int pt_index)
  {
    return scat_att.value (pt_index)
      + (1. - d2p_att.value (pt_index))
      + (1. - elev_att.value (pt_index));
  }
};

template <typename Kernel>
class Facade_classification_type : public Abstract_classification_type
{
  Scatter_segmentation_attribute<Kernel>& scat_att;
  Distance_to_plane_segmentation_attribute<Kernel>& d2p_att;
  Elevation_segmentation_attribute<Kernel>& elev_att;
public:
  Facade_classification_type (Scatter_segmentation_attribute<Kernel>& scat_att,
                              Distance_to_plane_segmentation_attribute<Kernel>& d2p_att,
                              Elevation_segmentation_attribute<Kernel>& elev_att)
    : scat_att (scat_att)
    , d2p_att (d2p_att)
    , elev_att (elev_att)
  { }

  virtual double data_term_computation (int pt_index)
  {
    return scat_att.value (pt_index)
      + d2p_att.value (pt_index)
      + (1. - elev_att.value (pt_index));
  }
};


template <typename Kernel>
class Building_classification_type : public Abstract_classification_type
{
  Scatter_segmentation_attribute<Kernel>& scat_att;
  Distance_to_plane_segmentation_attribute<Kernel>& d2p_att;
  Elevation_segmentation_attribute<Kernel>& elev_att;
public:
  Building_classification_type (Scatter_segmentation_attribute<Kernel>& scat_att,
                                Distance_to_plane_segmentation_attribute<Kernel>& d2p_att,
                                Elevation_segmentation_attribute<Kernel>& elev_att)
    : scat_att (scat_att)
    , d2p_att (d2p_att)
    , elev_att (elev_att)
  { }

  virtual double data_term_computation (int pt_index)
  {
    return scat_att.value (pt_index)
      + (1. - elev_att.value (pt_index));
  }
};




} // namespace CGAL

#endif // CGAL_POINT_SET_CLASSIFICATION_H

