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
#include <cstring>
#include <cassert>
#include <vector>
#include <list>
#include <set>

#include <CGAL/Default_diagonalize_traits.h>
#include <CGAL/centroid.h>
#include <CGAL/compute_average_spacing.h>


namespace CGAL {

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using I::misc::Ply;
typedef unsigned char byte;

//#define VERBOSE
#if defined(VERBOSE)
#define DEBUG_CERR std::cerr
#else
#define DEBUG_CERR std::ostream(0)
#endif

#if defined(VERBOSE)
inline void display_vect_stat (std::string name, std::vector<double> vect)
{
  std::sort (vect.begin(), vect.end());
  
  std::cerr << name << std::endl
            << " * Min = " << vect.front() << std::endl
            << " * 10% = " << vect[vect.size()/10] << std::endl
            << " * 50% = " << vect[vect.size()/2] << std::endl
            << " * 90% = " << vect[9 * vect.size()/10] << std::endl
            << " * Max = " << vect.back() << std::endl;

  double mean = 0.;
  for (std::size_t i = 0; i < vect.size(); ++ i)
    mean += vect[i];
  mean /= vect.size();
  std::cerr << " -> Mean = " << mean << std::endl;

  double standard_dev = 0.;
  for (std::size_t i = 0; i < vect.size(); ++ i)
    standard_dev += (vect[i] - mean) * (vect[i] - mean);
  standard_dev = std::sqrt (standard_dev / vect.size());
  std::cerr << " -> Standard deviation " << standard_dev << std::endl;

}
#else
inline void display_vect_stat (std::string, std::vector<double>) { }
#endif

inline void normalize (std::vector<double>& vect)
{
  // double mean = 0.;
  // for (std::size_t i = 0; i < vect.size(); ++ i)
  //   mean += vect[i];
  // mean /= vect.size();

  // double standard_dev = 0.;
  // for (std::size_t i = 0; i < vect.size(); ++ i)
  //   standard_dev += (vect[i] - mean) * (vect[i] - mean);
  // standard_dev = std::sqrt (standard_dev / vect.size());
  
  // for (std::size_t i = 0; i < vect.size(); ++ i)
  //   vect[i] = (vect[i] - mean) / (1.5 * standard_dev) + 1.5;
}

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


HSV_Color rgb_to_hsv (const Color& c)
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
  return {{ H, S, V }};
}



class ScanLidar;

class Abstract_segmentation_class
{
public:

  virtual double data_term_computation (int pt_index) = 0;

};

class ScanLidar
{
  struct HPoint {
    Point_d position;
    unsigned char echo; 
    int ind_x;
    int ind_y;
    int ultimate_label; 
    unsigned char AE_label;
    Color color;
  };

public:

  std::vector<Abstract_segmentation_class*> segmentation_classes;
  
  bool is_echo_given;
  bool is_normal_given;
  bounding_3 BBox_scan;

  vector<HPoint> HPS;
  Image< std::vector < int > ,2> grid_HPS;

  //Hpoints attributes
  vector < Vector_3 > normal;
  vector < CGAL::cpp11::array<double,3> > eigenvalues;
  std::vector<Plane_3> planes;
  std::vector < std::vector < int > > spherical_neighborhood;
  std::vector < std::vector < int > > cylindrical_neighborhood;


  Image<bool,2> Mask;                     //imagg
  Image<float,2> DTM; //a enregistrer ?   //imagg           

  /////    General parameters   //////
  double OUTLIER_PERCENTAGE_RM; 
  double RADIUS_NEIGHBORS_MAX; 
  double GRID_RESOLUTION;

  ////   Point cloud classification   /////
  double RADIUS_DTM; 
  double RADIUS_FILLING; 
  double SIGMA_NON_PLANARITY;  
  double SIGMA_GROUPING_1; 
  double SIGMA_GROUPING_2;
  double SEGMENTATION_SMOOTH_COEFFICIENT;
  double SCATTER_REGULARIZATION_VALUE; 
  int NB_LOCAL_NEIGHBORS; 

  ///////     file saving      ////////
  bool SAVE_POINT_ATTRIBUTES;
  bool SAVE_LABELING;



  double potts_computation(int l1, int l2) const{	
		
    double res=0;
    double smooth_seg=SEGMENTATION_SMOOTH_COEFFICIENT;

    if(l1!=l2) res=1; 

    return smooth_seg*res;
		
  }



  class Edge_score
  {
    const ScanLidar& M;

  public:
    Edge_score(const ScanLidar& _M) : M(_M) {}

    float compute(int, int, int l1, int l2)
    {
      return M.potts_computation(l1,l2); 
    }
  };


  class Facet_score
  {
    const ScanLidar& M;

  public:
    Facet_score(const ScanLidar& _M) : M(_M) {}

    double compute(int s, int l)
    {
      return M.segmentation_classes[l]->data_term_computation(s);
    }
  };


  void load_parameters(std::map<std::string, double>& values)
  {
    GRID_RESOLUTION = values[std::string("GRID_RESOLUTION")];

    OUTLIER_PERCENTAGE_RM = values[std::string("OUTLIER_PERCENTAGE_RM")];
    //energy weights
    SIGMA_NON_PLANARITY = values[std::string("SIGMA_NON_PLANARITY")];
        
    SIGMA_GROUPING_1 = values[std::string("SIGMA_GROUPING_1")];
    SIGMA_GROUPING_2 = values[std::string("SIGMA_GROUPING_2")];
    SEGMENTATION_SMOOTH_COEFFICIENT = values[std::string("SEGMENTATION_SMOOTH_COEFFICIENT")];
    //others
    RADIUS_DTM = values[std::string("RADIUS_DTM")];
    SCATTER_REGULARIZATION_VALUE = values[std::string("SCATTER_REGULARIZATION_VALUE")];
    //saving
    SAVE_POINT_ATTRIBUTES = values[std::string("SAVE_POINT_ATTRIBUTES")];
    SAVE_LABELING = values[std::string("SAVE_LABELING")];

    RADIUS_NEIGHBORS_MAX = values[std::string("RADIUS_NEIGHBORS_MAX")];
    NB_LOCAL_NEIGHBORS = values[std::string("NB_LOCAL_NEIGHBORS")];
  }

  
  bool load_points(string filename)
  {

    clock_t t;
    t = clock();

    is_normal_given=false;
    DEBUG_CERR << endl<<"Loading point cloud from '" << filename;
    /*if (binary) DEBUG_CERR << "' (bin,"; 
      else DEBUG_CERR << "' (ascii,"; 
      DEBUG_CERR.flush();*/

    HPS.clear();
    std::vector< Vector_3 > normals_temp;
    std::vector< HPoint > points_temp;
    Ply ply;

    if (!ply.open(filename)){ DEBUG_CERR << "Cannot read '" << filename << "'" << endl; return false;}

    for (Ply::ElementsIterator it = ply.elements_begin(); it != ply.elements_end(); ++it){
      const Ply::Element& element = *it;

      if (element.name() != "vertex"){
        if (!ply.skip(element)){ ply.close(); return false;}
        continue;
      }

      /* Check if the echo is available */
      if (element.num_properties()>3){
        if ((element.property(3).name() == "echo")){is_echo_given=true; DEBUG_CERR<<" with echo number) ";}else{is_echo_given=false; DEBUG_CERR<<" (WITHOUT ECHO NUMBER) ";}
      }else{is_echo_given=false; DEBUG_CERR<<" without echo number) ";}


      /* Check the properties are what we expect */
      if ((element.property(0).name() != "x")||(element.property(1).name() != "y")||(element.property(2).name() != "z")){DEBUG_CERR << "Unexpected vertex properties in the PLY file" << endl; ply.close(); return false;}


      size_t num_vertices = element.count();
      points_temp.resize(num_vertices);
	
      if (element.num_properties() == 3) {
        for (size_t i=0; i<num_vertices; i++){
          double x,y,z;

          if ((!ply.read(element.property(0), x))||(!ply.read(element.property(1), y))||(!ply.read(element.property(2), z))){DEBUG_CERR << "error while reading (pos) vertex " << i+1 << endl; ply.close(); return false;}	
				
          points_temp[i].position = Point_d(x,y,z);
          points_temp[i].echo = 1;
        }
      }

      else if (element.num_properties() == 4) {
        for (size_t i=0; i<num_vertices; i++){
          double x,y,z;
          double echo;

          if ((!ply.read(element.property(0), x))||(!ply.read(element.property(1), y))||(!ply.read(element.property(2), z))){DEBUG_CERR << "error while reading (pos) vertex " << i+1 << endl; ply.close(); return false;}	
          if (!ply.read(element.property(3), echo)){DEBUG_CERR << "error while reading attribut " << i+1 << endl; ply.close(); return false;}
				
          points_temp[i].position = Point_d(x,y,z);
          points_temp[i].echo = echo;
        }
      }

      else if (element.num_properties() == 5) {
        for (size_t i=0; i<num_vertices; i++){
          double x,y,z;
          double ampli;
          int echo;

          if ((!ply.read(element.property(0), x))||(!ply.read(element.property(1), y))||(!ply.read(element.property(2), z))){DEBUG_CERR << "error while reading (pos) vertex " << i+1 << endl; ply.close(); return false;}
          if ((!ply.read(element.property(3), ampli))||(!ply.read(element.property(4), echo))){DEBUG_CERR << "error while reading attribut " << i+1 << endl; ply.close(); return false;}
				
          points_temp[i].position = Point_d(x,y,z);
          points_temp[i].echo = echo;
        }
      }
		
      else if (element.num_properties() == 6) {
        for (size_t i=0; i<num_vertices; i++){
          double x,y,z;
          double nx,ny,nz;

          if ((!ply.read(element.property(0), x))||(!ply.read(element.property(1), y))||(!ply.read(element.property(2), z))){DEBUG_CERR << "error while reading (pos) vertex " << i+1 << endl; ply.close(); return false;}	
          if ((!ply.read(element.property(3), nx))||(!ply.read(element.property(4), ny))||(!ply.read(element.property(5), nz))){DEBUG_CERR << "error while reading attribut " << i+1 << endl; ply.close(); return false; }
				
          points_temp[i].position = Point_d(x,y,z);
          normals_temp.push_back(Vector_3(nx,ny,nz));
        }
        is_normal_given=true;
      }

      else if (element.num_properties()>6) {
        for (size_t i=0; i<num_vertices; i++){
          double x,y,z;
          double echo;
          double r, g, b;

          if ((!ply.read(element.property(0), x))||(!ply.read(element.property(1), y))||(!ply.read(element.property(2), z))){DEBUG_CERR << "error while reading (pos) vertex " << i+1 << endl; ply.close(); return false;}	
          if (!ply.read(element.property(3), echo)){DEBUG_CERR << "error while reading attribut " << i+1 << endl; ply.close(); return false;}
          if ((!ply.read(element.property(4), r))||(!ply.read(element.property(5), g))||(!ply.read(element.property(6), b))){DEBUG_CERR << "error while reading (color) vertex " << i+1 << endl; ply.close(); return false;}	
				
          points_temp[i].position = Point_d(x,y,z);
          points_temp[i].echo = echo;
          points_temp[i].color = {{ (unsigned char)r, (unsigned char)g, (unsigned char)b }};
        }
      }

      else {DEBUG_CERR << "Unexpected vertex properties in the PLY file" << endl; ply.close(); return false;}
    }


    /* Remove outliers */
    double Zmean=0; 
    double echo_mean=0;
    for(int i=0;i<(int)points_temp.size();i++){
      Zmean+=(double)points_temp[i].position.z(); 
      echo_mean+=(double)points_temp[i].echo;
    }
    Zmean/=(double)points_temp.size();
    DEBUG_CERR<<"echo: "<<echo_mean/points_temp.size()<<endl;
	
    double StDev=0;
    for(int i=0;i<(int)points_temp.size();i++) StDev+=pow((double)points_temp[i].position.z()-Zmean,2);
    StDev=sqrt(StDev/(double)points_temp.size());

    int count=0;
    double threshold_StDev=3;
    do {
      count=0;
      threshold_StDev+=0.5;
      for(std::size_t i=0;i<points_temp.size();i++){
        if(std::fabs((double)points_temp[i].position.z()-Zmean)>threshold_StDev*StDev) count++;
      }
    }while((double)count/points_temp.size()>OUTLIER_PERCENTAGE_RM);

    for(std::size_t i=0;i<points_temp.size();i++){
      if(std::fabs((double)points_temp[i].position.z()-Zmean)<=threshold_StDev*StDev){
        HPS.push_back(points_temp[i]);
        if(is_normal_given) normal.push_back(normals_temp[i]);
      }
    }
    DEBUG_CERR <<points_temp.size() << " initial points: "<<HPS.size()<<" inliers and " <<points_temp.size()-HPS.size()<<" outliers"<< endl<<"-> OK ( "<<((float)clock()-t)/CLOCKS_PER_SEC<<" sec )"<< endl;

    points_temp.clear();
    ply.close();

    return true;
  }


  double average_spacing () const
  {
    std::vector<Point_d> pts;
    for (std::size_t i = 0; i < HPS.size(); ++ i)
      pts.push_back (HPS[i].position);
    return CGAL::compute_average_spacing<CGAL::Sequential_tag> (pts.begin(), pts.end(), 6);
  }
  

  bool initialization(int phase = 0) // for training
  {
    clock_t t;
    t = clock();


    DEBUG_CERR<<endl<<"Initialization: ";

    //1-Neighborhood computation and reset the attributes of the structure points
    DEBUG_CERR<<"spherical neighborhood..";
    spherical_neighborhood.clear();
    eigenvalues.clear();
    planes.clear();
    
    std::map<Point_d,int> map_indice_point;
    for(std::size_t i=0;i<HPS.size();i++){
      Point_d pt=HPS[i].position;
      map_indice_point[pt]=i;
	
    }

    std::vector<Point_d> list_points;
    for(int i=0;i<(int)HPS.size();i++){
      Point_d pt=HPS[i].position;
      list_points.push_back(pt);
    }

    Tree tree(list_points.begin(), list_points.end());

    std::size_t nb_neigh = 0;
    for(int i=0;i<(int)HPS.size();i++){
      const Point_d& query=HPS[i].position;
      std::vector<Point_d> neighbors;
      
      if (phase == 1 && HPS[i].echo == 0)
        {
          compute_principal_curvature (query, neighbors);
          continue;
        }
      
      Fuzzy_sphere fs(query, RADIUS_NEIGHBORS_MAX, 0);
      tree.search(std::back_inserter(neighbors), fs);
      nb_neigh += neighbors.size();

      compute_principal_curvature (query, neighbors);
    }

    DEBUG_CERR<<"ok";

    if (phase == 1)
      return true;

    //2-creation of the bounding box
    DEBUG_CERR<<", bounding box..";
    BBox_scan = CGAL::bounding_box(list_points.begin(), list_points.end());
    DEBUG_CERR<<"ok";

    //3-creation grille_points
    DEBUG_CERR<<", planimetric grid of HPS..";
    Image< std::vector < int > ,2> tess((int)((BBox_scan.xmax()-BBox_scan.xmin())/GRID_RESOLUTION)+1,(int)((BBox_scan.ymax()-BBox_scan.ymin())/GRID_RESOLUTION)+1);
    grid_HPS=tess;

    for(int i=0;i<(int)HPS.size();i++){

      //for each 3D point, its coordinates in the grid are inserted in its Hpoint structure
      HPS[i].ind_x=(int)((HPS[i].position.x()-BBox_scan.xmin())/GRID_RESOLUTION);
      HPS[i].ind_y=(int)((HPS[i].position.y()-BBox_scan.ymin())/GRID_RESOLUTION);

      //index of points are collected in grid_HPS
      std::vector < int > temp;
      temp=grid_HPS(HPS[i].ind_x,HPS[i].ind_y);
      temp.push_back(i);
      grid_HPS(HPS[i].ind_x,HPS[i].ind_y)=temp;
    }
    DEBUG_CERR<<"ok";


    //4-Mask creation
    DEBUG_CERR<<", planimetric mask..";
    Image< bool ,2> masktp((int)((BBox_scan.xmax()-BBox_scan.xmin())/GRID_RESOLUTION)+1,(int)((BBox_scan.ymax()-BBox_scan.ymin())/GRID_RESOLUTION)+1);
    Mask=masktp;

    DEBUG_CERR << "(" << Mask.height() << "x" << Mask.width() << ")" << std::endl;
    int square=(int)16*(RADIUS_NEIGHBORS_MAX/GRID_RESOLUTION)+1;
    int nb_true = 0;
    for (int j=0;j<(int)Mask.height();j++){	
      for (int i=0;i<(int)Mask.width();i++){	
        //Mask(i,j)=false;
        if(grid_HPS(i,j).size()>0) { Mask(i,j)=true; nb_true ++; }
        else{Mask(i,j)=false;}
      }
    }

    Image< bool ,2> Mask_tmp ((int)((BBox_scan.xmax()-BBox_scan.xmin())/GRID_RESOLUTION)+1,(int)((BBox_scan.ymax()-BBox_scan.ymin())/GRID_RESOLUTION)+1);

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
    //    if(!is_normal_given){DEBUG_CERR<<", normals.."; compute_normal(); DEBUG_CERR<<"ok";}

    DEBUG_CERR<<endl<<"-> OK ( "<<((float)clock()-t)/CLOCKS_PER_SEC<<" sec )"<< endl;

    return true;
  }



  bool compute_normal()
  {

    for(int k=0;k<(int)HPS.size();k++){

      std::list<Point_d> list_pts;
      for(int n=0;n<(int)spherical_neighborhood[k].size();n++){
				
        if(n<NB_LOCAL_NEIGHBORS){
          Point_d pt_temp=HPS[spherical_neighborhood[k][n]].position;
          list_pts.push_back(pt_temp);
        }
        else break;
      }
			
      Plane_3 plane_temp;
      linear_least_squares_fitting_3(list_pts.begin(),list_pts.end(),plane_temp, CGAL::Dimension_tag<0>());

      Vector_3 normal_temp=plane_temp.orthogonal_vector();

      for(int n=0;n<(int)spherical_neighborhood[k].size();n++){
        if(k>spherical_neighborhood[k][n] && n<NB_LOCAL_NEIGHBORS){
          if(normal[spherical_neighborhood[k][n]] * normal_temp <0 ) normal_temp=-normal_temp;
          break;
        }
        else break;
      }

      FT normal_temp_norm=1/sqrt(normal_temp.squared_length());
      if(normal_temp.z()>0) normal.push_back(normal_temp_norm*normal_temp);
      else normal.push_back(-normal_temp_norm*normal_temp);
    }

    return true;
  }





  void compute_principal_curvature(const Point_d& point, std::vector<Point_d>& neighborhood)
  {
    if (neighborhood.size() == 0)
      {
        eigenvalues.push_back ({{ 0., 0., 0. }});
        planes.push_back (Plane_3 (point, Vector_3 (0., 0., 1.)));
        return;
      }
    Point_d centroid = CGAL::centroid (neighborhood.begin(), neighborhood.end());

    CGAL::cpp11::array<double, 6> covariance = {{ 0., 0., 0., 0., 0., 0. }};
      
    for (std::size_t i = 0; i < neighborhood.size(); ++ i)
      {
        Vector_3 d = neighborhood[i] - centroid;
        covariance[0] += d.x () * d.x ();
        covariance[1] += d.x () * d.y ();
        covariance[2] += d.x () * d.z ();
        covariance[3] += d.y () * d.y ();
        covariance[4] += d.y () * d.z ();
        covariance[5] += d.z () * d.z ();
      }


    eigenvalues.push_back ({{ 0., 0., 0. }});
    CGAL::cpp11::array<double, 9> eigenvectors = {{ 0., 0., 0.,
                                                    0., 0., 0.,
                                                    0., 0., 0. }};

    CGAL::Default_diagonalize_traits<double,3>::diagonalize_selfadjoint_covariance_matrix
      (covariance, eigenvalues.back(), eigenvectors);

    Plane_3 plane;
    CGAL::linear_least_squares_fitting_3 (neighborhood.begin(),neighborhood.end(),plane, CGAL::Dimension_tag<0>());
    planes.push_back (plane);
  }



  bool labeling_PC(){

    int number_of_alpha_expansions = 2;
	
    // graph description
    DEBUG_CERR << "Graph description... ";
    GCoptimizationGeneralGraph<float, float> *gc= new GCoptimizationGeneralGraph<float, float>
      ((int)HPS.size(),(int)(segmentation_classes.size()));
	
    for(int k=0;k<(int)HPS.size();k++){
      vector < int > index_of_neighboring_points;
      for(int n=0;n<(int)spherical_neighborhood[k].size();n++)
        if(n<NB_LOCAL_NEIGHBORS) index_of_neighboring_points.push_back(spherical_neighborhood[k][n]);
		
      for(int k_vois=0; k_vois<(int)index_of_neighboring_points.size(); k_vois++)
        gc->setNeighbors(k, index_of_neighboring_points[k_vois]);
    }

    gc->specializeDataCostFunctor(Facet_score(*this));
    gc->specializeSmoothCostFunctor(Edge_score(*this));
	
    
    // data term initialisation
    DEBUG_CERR << "Data term initialization... ";
    for(int s=0;s<(int)HPS.size();s++){
			
      int nb_class_best=0; 

      double val_class_best= segmentation_classes[0]->data_term_computation(s);

      for(std::size_t k = 1; k < segmentation_classes.size(); ++ k)
        {
          double value = segmentation_classes[k]->data_term_computation(s);
          if(val_class_best > value)
            {
              val_class_best = value;
              nb_class_best=k;
            }
        }
      gc->setLabel(s,nb_class_best);

    }

    // Optimize labeling
    DEBUG_CERR << "Optimize labeling... ";
    for (int iter=0; iter<number_of_alpha_expansions; iter++){
			
      for (std::size_t i=0; i< segmentation_classes.size();i++){
				
        // Compute vector of active sites
        vector<int> sites;
        for(int s=0;s<(int)HPS.size();s++) sites.push_back(s);
			
        // Compute alpha expansion for this label on these sites
        gc->alpha_expansion((int)i, &(sites[0]), sites.size());
      }
    }

    // Store labeling
    DEBUG_CERR << "Store labeling... ";
    int count1=0; 
    int count2=0; 
    int count3=0; 
    int count4=0;
    for(int s=0;s<(int)HPS.size();s++){
      HPS[s].AE_label=gc->whatLabel(s); 
      if(gc->whatLabel(s)==0) count1++;
      else if(gc->whatLabel(s)==1) count2++;
      else if(gc->whatLabel(s)==2) count3++;
      else count4++;
    }

    delete gc;

    DEBUG_CERR<<" "<<(double)100*count1/HPS.size()<<"% vegetation, "<<(double)100*count2/HPS.size()<<"% ground, "<<(double)100*count3/HPS.size()<<"% building, "<<(double)100*count4/HPS.size()<<"% clutter"<<endl;
    
    if(SAVE_LABELING)
      {
        CGAL::Color bblue(130,130,255);
        CGAL::Color rred(220,140,140);
        CGAL::Color orange(226,247,130);
        CGAL::Color white(225,225,225);

        std::vector < CGAL::Color > vecteur_col; 
        vecteur_col.push_back(rred);
        vecteur_col.push_back(orange);
        vecteur_col.push_back(bblue);
        vecteur_col.push_back(white);


        std::string nom_classif("PC_class_all.ply");
        std::string nom_classif1("PC_class_vegetation.ply");
        std::string nom_classif2("PC_class_ground.ply");
        std::string nom_classif3("PC_class_building.ply");
        std::string nom_classif4("PC_class_clutter.ply");

        std::vector< Point_d > pointtts0; std::vector< CGAL::Color > colors000; 
        std::vector< Point_d > pointAE1; std::vector< CGAL::Color > colorsAE1;
        std::vector< Point_d > pointAE2; std::vector< CGAL::Color > colorsAE2;
        std::vector< Point_d > pointAE3; std::vector< CGAL::Color > colorsAE3;
        std::vector< Point_d > pointAE4; std::vector< CGAL::Color > colorsAE4;

        for(int i=0; i<(int)HPS.size();i++){
          Point_d pt=HPS[i].position; 
          pointtts0.push_back(pt); 
          colors000.push_back(vecteur_col[HPS[i].AE_label]);

          if(HPS[i].AE_label==0){pointAE1.push_back(pt); colorsAE1.push_back(vecteur_col[HPS[i].AE_label]);}
          else if(HPS[i].AE_label==1){pointAE2.push_back(pt); colorsAE2.push_back(vecteur_col[HPS[i].AE_label]);}
          else if(HPS[i].AE_label==2){pointAE3.push_back(pt); colorsAE3.push_back(vecteur_col[HPS[i].AE_label]);}
          else{pointAE4.push_back(pt); colorsAE4.push_back(vecteur_col[HPS[i].AE_label]);}
        }
        colorpointset2ply(nom_classif.c_str(),pointtts0,colors000);
        colorpointset2ply(nom_classif1.c_str(),pointAE1,colorsAE1);
        colorpointset2ply(nom_classif2.c_str(),pointAE2,colorsAE2);
        colorpointset2ply(nom_classif3.c_str(),pointAE3,colorsAE3);
        colorpointset2ply(nom_classif4.c_str(),pointAE4,colorsAE4);
      }
	
    return true;
  }


  bool quick_labeling_PC(bool log = false)
  {

    // data term initialisation
    DEBUG_CERR << "Labeling... ";

    std::vector<std::vector<double> > values;
    if (log)
      values.resize (segmentation_classes.size());

    
    int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
    for(int s=0;s<(int)HPS.size();s++){
			
      int nb_class_best=0; 

      double val_class_best = std::numeric_limits<double>::max();

      for(std::size_t k = 0; k < segmentation_classes.size(); ++ k)
        {
          double value = segmentation_classes[k]->data_term_computation(s);
          
          if (log)
            values[k].push_back(value);
              
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

    if (log)
      {
        for (std::size_t i = 0; i < values.size(); ++ i)
          {
            std::string name;
            if (i == 0)
              name = "vegetation";
            else if (i == 1)
              name = "ground";
            else if (i == 2)
              name = "roof";
            else if (i == 3)
              name = "facade";
            save_PC_attribute (name, values[i]);
          }
      }
    
    DEBUG_CERR<<" "<<(double)100*count1/HPS.size()<<"% vegetation, "<<(double)100*count2/HPS.size()<<"% ground, "<<(double)100*count3/HPS.size()<<"% roof, "<<(double)100*count4/HPS.size()<<"% clutter"<<endl;
    
    if(SAVE_LABELING)
      {
        CGAL::Color green(165, 226, 104);
        CGAL::Color grey(169, 152, 126);
        CGAL::Color yellow(255, 254, 106);
        CGAL::Color orange(255, 135, 77);

        std::vector < CGAL::Color > vecteur_col; 
        vecteur_col.push_back(green);
        vecteur_col.push_back(grey);
        vecteur_col.push_back(orange);
        vecteur_col.push_back(yellow);


        std::string nom_classif("PC_class_all.ply");
        std::string nom_classif1("PC_class_vegetation.ply");
        std::string nom_classif2("PC_class_ground.ply");
        std::string nom_classif3("PC_class_roof.ply");
        std::string nom_classif4("PC_class_facade.ply");

        std::vector< Point_d > pointtts0; std::vector< CGAL::Color > colors000; 
        std::vector< Point_d > pointAE1; std::vector< CGAL::Color > colorsAE1;
        std::vector< Point_d > pointAE2; std::vector< CGAL::Color > colorsAE2;
        std::vector< Point_d > pointAE3; std::vector< CGAL::Color > colorsAE3;
        std::vector< Point_d > pointAE4; std::vector< CGAL::Color > colorsAE4;

        for(int i=0; i<(int)HPS.size();i++){
          Point_d pt=HPS[i].position; 
          pointtts0.push_back(pt); 
          colors000.push_back(vecteur_col[HPS[i].AE_label]);

          if(HPS[i].AE_label==0)
            {
              pointAE1.push_back(pt);
              colorsAE1.push_back(vecteur_col[HPS[i].AE_label]);
            }
          else if(HPS[i].AE_label==1)
            {
              pointAE2.push_back(pt);
              colorsAE2.push_back(vecteur_col[HPS[i].AE_label]);
            }
          else if(HPS[i].AE_label==2)
            {
              pointAE3.push_back(pt);
              colorsAE3.push_back(vecteur_col[HPS[i].AE_label]);
            }
          else
            {
              pointAE4.push_back(pt);
              colorsAE4.push_back(vecteur_col[HPS[i].AE_label]);
            }
        }
        colorpointset2ply(nom_classif.c_str(),pointtts0,colors000);
        colorpointset2ply(nom_classif1.c_str(),pointAE1,colorsAE1);
        colorpointset2ply(nom_classif2.c_str(),pointAE2,colorsAE2);
        colorpointset2ply(nom_classif3.c_str(),pointAE3,colorsAE3);
        colorpointset2ply(nom_classif4.c_str(),pointAE4,colorsAE4);
      }
	
    return true;
  }


  template <typename SegmentationAttribute>
  bool save_PC_attribute (SegmentationAttribute& att)
  {
    std::vector< Point_d > pointtts;
    std::vector< CGAL::Color > colors;

    std::stringstream ss;
    ss << "F_" << att.id() << ".ply";

    for(int i=0; i<(int)HPS.size();i++){
      Point_d pt=HPS[i].position; pointtts.push_back(pt);
      double v = att.value(i);
      double cc = std::max (0., std::min(1., v))*255;
      CGAL::Color colo (255-(int)cc, 255-(int)cc, 255); 
      colors.push_back(colo);
    }

    colorpointset2ply(ss.str().c_str(),pointtts,colors);
    return true;
  }


  bool save_PC_attribute (std::string name, std::vector<double>& val)
  {
    std::vector< Point_d > pointtts;
    std::vector< CGAL::Color > colors;

    std::stringstream ss;
    ss << "PROBA_" << name << ".ply";

    for(int i=0; i<(int)HPS.size();i++){
      Point_d pt=HPS[i].position; pointtts.push_back(pt);
      double v = val[i];
      double cc = std::max (0., std::min(1., v))*255;
      CGAL::Color colo (255-(int)cc, 255-(int)cc, 255); 
      colors.push_back(colo);
    }

    colorpointset2ply(ss.str().c_str(),pointtts,colors);
    return true;
  }

  bool smoothed_labeling_PC()
  {

    // data term initialisation
    DEBUG_CERR << "Labeling... ";

    std::vector<std::vector<double> > values;
    values.resize (segmentation_classes.size());
    
    std::map<Point_d, std::size_t> map_p2i;
    std::vector<Point_d> list_points;    
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
        const Point_d& query=HPS[s].position;
        std::vector<Point_d> neighbors;
      
        Fuzzy_sphere fs(query, RADIUS_NEIGHBORS_MAX, 0);
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

    for (std::size_t i = 0; i < values.size(); ++ i)
      {
        std::string name;
        if (i == 0)
          name = "vegetation";
        else if (i == 1)
          name = "ground";
        else if (i == 2)
          name = "roof";
        else if (i == 3)
          name = "facade";
        save_PC_attribute (name, smoothed_values[i]);
      }

    DEBUG_CERR<<" "<<(double)100*count1/HPS.size()<<"% vegetation, "<<(double)100*count2/HPS.size()<<"% ground, "<<(double)100*count3/HPS.size()<<"% roof, "<<(double)100*count4/HPS.size()<<"% clutter"<<endl;
    
    if(SAVE_LABELING)
      {
        CGAL::Color green(165, 226, 104);
        CGAL::Color grey(169, 152, 126);
        CGAL::Color yellow(255, 254, 106);
        CGAL::Color orange(255, 135, 77);

        std::vector < CGAL::Color > vecteur_col; 
        vecteur_col.push_back(green);
        vecteur_col.push_back(grey);
        vecteur_col.push_back(orange);
        vecteur_col.push_back(yellow);


        std::string nom_classif("PC_class_all.ply");
        std::string nom_classif1("PC_class_vegetation.ply");
        std::string nom_classif2("PC_class_ground.ply");
        std::string nom_classif3("PC_class_roof.ply");
        std::string nom_classif4("PC_class_facade.ply");

        std::vector< Point_d > pointtts0; std::vector< CGAL::Color > colors000; 
        std::vector< Point_d > pointAE1; std::vector< CGAL::Color > colorsAE1;
        std::vector< Point_d > pointAE2; std::vector< CGAL::Color > colorsAE2;
        std::vector< Point_d > pointAE3; std::vector< CGAL::Color > colorsAE3;
        std::vector< Point_d > pointAE4; std::vector< CGAL::Color > colorsAE4;

        for(int i=0; i<(int)HPS.size();i++){
          Point_d pt=HPS[i].position; 
          pointtts0.push_back(pt); 
          colors000.push_back(vecteur_col[HPS[i].AE_label]);

          if(HPS[i].AE_label==0)
            {
              pointAE1.push_back(pt);
              colorsAE1.push_back(vecteur_col[HPS[i].AE_label]);
            }
          else if(HPS[i].AE_label==1)
            {
              pointAE2.push_back(pt);
              colorsAE2.push_back(vecteur_col[HPS[i].AE_label]);
            }
          else if(HPS[i].AE_label==2)
            {
              pointAE3.push_back(pt);
              colorsAE3.push_back(vecteur_col[HPS[i].AE_label]);
            }
          else
            {
              pointAE4.push_back(pt);
              colorsAE4.push_back(vecteur_col[HPS[i].AE_label]);
            }
        }
        colorpointset2ply(nom_classif.c_str(),pointtts0,colors000);
        colorpointset2ply(nom_classif1.c_str(),pointAE1,colorsAE1);
        colorpointset2ply(nom_classif2.c_str(),pointAE2,colorsAE2);
        colorpointset2ply(nom_classif3.c_str(),pointAE3,colorsAE3);
        colorpointset2ply(nom_classif4.c_str(),pointAE4,colorsAE4);
      }
	
    return true;
  }

  bool point_cloud_classification(int method = 1, bool log = false){


    clock_t t;
    t = clock();
	
    DEBUG_CERR<<endl<<"Classification of the point cloud: ";
    if (method == 0)
      labeling_PC();
    else if (method == 1)
      quick_labeling_PC(log);
    else if (method == 2)
      smoothed_labeling_PC();
    else
      {
        std::cerr << "Unknown method number." << std::endl;
        abort();
      }

    DEBUG_CERR<<"-> OK ( "<<((float)clock()-t)/CLOCKS_PER_SEC<<" sec )"<< endl;

    return true;
  }



protected: 



};


class Scatter_segmentation_attribute
{
  std::vector<double> vegetation_attribute;
  
public:
  double weight;
  double mean;
  double max;

  Scatter_segmentation_attribute (ScanLidar& M, double weight, bool colors = false) : weight (weight)
  {
    Image<float,2> Vegetation(M.grid_HPS.width(), M.grid_HPS.height());

    for (int j=0;j<(int)M.DTM.height();j++)	
      for (int i=0;i<(int)M.DTM.width();i++)
        Vegetation(i,j)=0;
		
    int square = (int)0.5 * M.RADIUS_NEIGHBORS_MAX / M.GRID_RESOLUTION + 1;

    if(M.is_echo_given){
				
      for (int j=0;j<(int)M.grid_HPS.height();j++){	
        for (int i=0;i<(int)M.grid_HPS.width();i++){
						
          if(M.Mask(i,j)){

            int squareXmin=std::max(0,i-square);
            int squareXmax=std::min(M.grid_HPS.width()-1,i+square);
            int squareYmin=std::max(0,j-square);
            int squareYmax=std::min(M.grid_HPS.height()-1,j+square);
			
            int NB_echo_sup=0;
            int NB_echo_total=0;

            for(int k=squareXmin; k<=squareXmax; k++){
              for(int l=squareYmin; l<=squareYmax; l++){
									
                if(sqrt(pow((double)k-i,2)+pow((double)l-j,2))<=(double)0.5*M.RADIUS_NEIGHBORS_MAX/M.GRID_RESOLUTION){
										
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
    else if (!colors){
      for(int i=0;i<(int)M.HPS.size();i++){
        int I= M.HPS[i].ind_x;
        int J= M.HPS[i].ind_y;
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
          HSV_Color c = rgb_to_hsv (M.HPS[i].color);
          vegetation_attribute.push_back (std::exp (-(c[0] - 120) * (c[0] - 120) / (2 * 60 * 60))
                                          * std::exp (-(c[1] - 25) * (c[1] - 25) / (2 * 20 * 20))
                                          * std::exp (-(c[2] - 37) * (c[2] - 37) / (2 * 10 * 10)));
        }
      }
    normalize (vegetation_attribute);
    display_vect_stat ("Distance to plane:", vegetation_attribute);
    compute_mean_max (vegetation_attribute, mean, max);
    max *= 2;
  }

  double value (int pt_index)
  {
    return std::max (0., std::min (1., vegetation_attribute[pt_index] / weight));
  }

  const char* id() { return "scatter"; }
};

class Distance_to_plane_segmentation_attribute
{

  std::vector<double> distance_to_plane_attribute;
  
public:
  double weight;
  double mean;
  double max;
  
  Distance_to_plane_segmentation_attribute (ScanLidar& M, double weight) : weight (weight)
  {
    for(int i=0; i<(int)M.HPS.size(); i++)
      distance_to_plane_attribute.push_back (std::sqrt (CGAL::squared_distance (M.HPS[i].position, M.planes[i])));

    normalize (distance_to_plane_attribute);
    display_vect_stat ("Distance to plane:", distance_to_plane_attribute);
    compute_mean_max (distance_to_plane_attribute, mean, max);
    max *= 2;
  }

  double value (int pt_index)
  {
    return std::max (0., std::min (1., distance_to_plane_attribute[pt_index] / weight));
  }

  const char* id() { return "planimetry"; }
};

class Elevation_segmentation_attribute
{

  std::vector<double> elevation_attribute;
  
public:
  double weight;
  double mean;
  double max;
  
  Elevation_segmentation_attribute (ScanLidar& M, double weight) : weight (weight)
  {
    //DEM
    Image<float,2> DEM(M.grid_HPS.width(),M.grid_HPS.height());
    Image<float,2> DEMtemp(M.grid_HPS.width(),M.grid_HPS.height());

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

    int square=(int)M.RADIUS_NEIGHBORS_MAX/M.GRID_RESOLUTION+1;

    for (int j=0;j<DEM.height();j++){	
      for (int i=0;i<DEM.width();i++){
		
        if((M.Mask(i,j))&&(DEM(i,j)==0)){
				
          double distance_tot=0;	
          double val=0;
          int squareXmin=std::max(0,i-square);
          int squareXmax=std::min(DEM.width()-1,i+square);
          int squareYmin=std::max(0,j-square);
          int squareYmax=std::min(DEM.height()-1,j+square);
			
          for(int k=squareXmin; k<=squareXmax; k++){
            for(int l=squareYmin; l<=squareYmax; l++){
					
              double distance=sqrt(pow((double)i-k,2)+pow((double)j-l,2))*M.GRID_RESOLUTION;
					
              if((distance<=M.RADIUS_NEIGHBORS_MAX)&&(DEM(k,l)>0)&&(distance!=0)){
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

    for (int j=0;j<DEM.height();j++){	
      for (int i=0;i<DEM.width();i++){
		
        DEM(i,j)=DEMtemp(i,j);

        if(M.Mask(i,j)){
          if(DEM(i,j)>scale_byte_max_dem) scale_byte_max_dem=DEM(i,j);
          if(DEM(i,j)<scale_byte_min_dem) scale_byte_min_dem=DEM(i,j);
        }
      }
    }

    //TODO: smooth the DEM
    for (int j=1;j<DEM.height()-1;j++){	
      for (int i=1;i<DEM.width()-1;i++){
        if(M.Mask(i,j)){
          //DEM(i,j)=
        }
      }
    }

    //DTM computation
    int step=15;
    square=(int)(M.RADIUS_DTM/M.GRID_RESOLUTION)+1;
    Image<float,2> toto(M.grid_HPS.width(),M.grid_HPS.height());
    Image<float,2> im_Zfront(M.grid_HPS.width(),M.grid_HPS.height());
    M.DTM=toto;

    //round 1
    for (int j=0;j<(int)M.DTM.height();j++){	
      for (int i=0;i<(int)M.DTM.width();i++){
        M.DTM(i,j)=0;
        im_Zfront(i,j)=0;
      }
    }


    for (int j=step/2+1;j<(int)M.DTM.height();j=j+step){	
      for (int i=step+1/2;i<(int)M.DTM.width();i=i+step){
		
        //storage of the points in the disk
        std::vector < float > list_pointsZ;
        int squareXmin=std::max(0,i-square);
        int squareXmax=std::min(M.DTM.width()-1,i+square);
        int squareYmin=std::max(0,j-square);
        int squareYmax=std::min(M.DTM.height()-1,j+square);

        for(int k=squareXmin; k<=squareXmax; k++){
          for(int l=squareYmin; l<=squareYmax; l++){
			
            double distance=sqrt(pow((double)i-k,2)+pow((double)j-l,2))*M.GRID_RESOLUTION;
            if(distance<=M.RADIUS_DTM){
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
        int IsquareXmin=std::max(0,i-step/2);
        int IsquareXmax=std::min(M.DTM.width()-1,i+step/2);
        int JsquareYmin=std::max(0,j-step/2);
        int JsquareYmax=std::min(M.DTM.height()-1,j+step/2);
		
        if(M.DTM.width()-1-IsquareXmax<step) IsquareXmax=M.DTM.width()-1;
        if(M.DTM.height()-1-JsquareYmax<step) JsquareYmax=M.DTM.height()-1;
        if(IsquareXmin<step) IsquareXmin=0;
        if(JsquareYmin<step) JsquareYmin=0;

        for(int k=IsquareXmin; k<=IsquareXmax; k++){
          for(int l=JsquareYmin; l<=JsquareYmax; l++){
            if(M.Mask(k,l)){
              M.DTM(k,l)=G1;
              im_Zfront(k,l)=Gfront;
            }
          }
        }
      }
    }

    //Gaussian smoothness
    for(int ik=0;ik<10;ik++) {
      for(int j=0;j<(int)M.DTM.height();j++){ 
        for (int i=0;i<(int)M.DTM.width();i++) {
		
          int IsquareXmin=std::max(0,i-1); 
          int IsquareXmax=std::min(M.DTM.width()-1,i+1); 
          int JsquareYmin=std::max(0,j-1); 
          int JsquareYmax=std::min(M.DTM.height()-1,j+1); 
          int count=0; 

          for(int k=IsquareXmin; k<=IsquareXmax; k++){ 
            for(int l=JsquareYmin; l<=JsquareYmax; l++){ 
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

    for (int j=step/2+1;j<(int)M.DTM.height();j=j+step){	
      for (int i=step/2+1;i<(int)M.DTM.width();i=i+step){
		
        //stockage des nouveaux points
        std::vector < float > list_pointsZ;
        int squareXmin=std::max(0,i-square);
        int squareXmax=std::min(M.DTM.width()-1,i+square);
        int squareYmin=std::max(0,j-square);
        int squareYmax=std::min(M.DTM.height()-1,j+square);

        for(int k=squareXmin; k<=squareXmax; k++){
          for(int l=squareYmin; l<=squareYmax; l++){
				
            double distance=sqrt(pow((double)i-k,2)+pow((double)j-l,2))*M.GRID_RESOLUTION;
				
            if(distance<=M.RADIUS_DTM){
					
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
        int IsquareXmin=std::max(0,i-step/2);
        int IsquareXmax=std::min(M.DTM.width()-1,i+step/2);
        int JsquareYmin=std::max(0,j-step/2);
        int JsquareYmax=std::min(M.DTM.height()-1,j+step/2);
		
        if(M.DTM.width()-1-IsquareXmax<step) IsquareXmax=M.DTM.width()-1;
        if(M.DTM.height()-1-JsquareYmax<step) JsquareYmax=M.DTM.height()-1;
        if(IsquareXmin<step) IsquareXmin=0;
        if(JsquareYmin<step) JsquareYmin=0;

        for(int k=IsquareXmin; k<=IsquareXmax; k++){
          for(int l=JsquareYmin; l<=JsquareYmax; l++){
            if(M.Mask(k,l)){
              M.DTM(k,l)=G1; 
              im_Zfront(k,l)=Gfront;
            }
          }
        }
      }
    }

    //Gaussian smoothness
    for(int ik=0;ik<10;ik++) {
      for(int j=0;j<(int)M.DTM.height();j++){ 
        for (int i=0;i<(int)M.DTM.width();i++) {
		
          int IsquareXmin=std::max(0,i-1); 
          int IsquareXmax=std::min(M.DTM.width()-1,i+1); 
          int JsquareYmin=std::max(0,j-1); 
          int JsquareYmax=std::min(M.DTM.height()-1,j+1); 
          int count=0; 

          for(int k=IsquareXmin; k<=IsquareXmax; k++){ 
            for(int l=JsquareYmin; l<=JsquareYmax; l++){ 
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
    normalize (elevation_attribute);
    display_vect_stat ("Elevation:", elevation_attribute);
    compute_mean_max (elevation_attribute, mean, max);
    max *= 5;
  }

  double value (int pt_index)
  {
    return std::max (0., std::min (1., elevation_attribute[pt_index] / weight));
  }

  const char* id() { return "elevetation"; }

};


class Vegetation_segmentation_class : public Abstract_segmentation_class
{
  Scatter_segmentation_attribute& scat_att;
  Distance_to_plane_segmentation_attribute& d2p_att;
  Elevation_segmentation_attribute& elev_att;
public:

  Vegetation_segmentation_class (Scatter_segmentation_attribute& scat_att,
                                 Distance_to_plane_segmentation_attribute& d2p_att,
                                 Elevation_segmentation_attribute& elev_att)
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

class Ground_segmentation_class : public Abstract_segmentation_class
{
  Scatter_segmentation_attribute& scat_att;
  Distance_to_plane_segmentation_attribute& d2p_att;
  Elevation_segmentation_attribute& elev_att;
public:

  Ground_segmentation_class (Scatter_segmentation_attribute& scat_att,
                             Distance_to_plane_segmentation_attribute& d2p_att,
                             Elevation_segmentation_attribute& elev_att)
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

class Roof_segmentation_class : public Abstract_segmentation_class
{
  Scatter_segmentation_attribute& scat_att;
  Distance_to_plane_segmentation_attribute& d2p_att;
  Elevation_segmentation_attribute& elev_att;
public:
  Roof_segmentation_class (Scatter_segmentation_attribute& scat_att,
                           Distance_to_plane_segmentation_attribute& d2p_att,
                           Elevation_segmentation_attribute& elev_att)
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

class Facade_segmentation_class : public Abstract_segmentation_class
{
  Scatter_segmentation_attribute& scat_att;
  Distance_to_plane_segmentation_attribute& d2p_att;
  Elevation_segmentation_attribute& elev_att;
public:
  Facade_segmentation_class (Scatter_segmentation_attribute& scat_att,
                             Distance_to_plane_segmentation_attribute& d2p_att,
                             Elevation_segmentation_attribute& elev_att)
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


class Building_segmentation_class : public Abstract_segmentation_class
{
  Scatter_segmentation_attribute& scat_att;
  Distance_to_plane_segmentation_attribute& d2p_att;
  Elevation_segmentation_attribute& elev_att;
public:
  Building_segmentation_class (Scatter_segmentation_attribute& scat_att,
                               Distance_to_plane_segmentation_attribute& d2p_att,
                               Elevation_segmentation_attribute& elev_att)
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


class Clutter_segmentation_class : public Abstract_segmentation_class
{
public:

  virtual double data_term_computation (int pt_index)
  {
    // double PD=std::min(1., scan.distance_to_plane_attribute[pt_index] / scan.SIGMA_NON_PLANARITY); 
    // double PNB=0.; if(scan.spherical_neighborhood[pt_index].size() < scan.SIGMA_GROUPING_1) PNB=1;
    // double PZ2=0; if(scan.elevation_attribute[pt_index]<-10){PZ2=1;} 
    // double PDL=scan.distance_to_line_attribute[pt_index]; if(PDL<scan.SIGMA_GROUPING_2){PDL=1;} else{PDL=0;}
    // double result=std::max(PNB,PD*PZ2); result=std::max(result,PDL);
    return 0;
  }
};


} // namespace CGAL

#endif // CGAL_POINT_SET_CLASSIFICATION_H
