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
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Default_diagonalize_traits.h>
#include <CGAL/centroid.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <CGAL/Data_classification/Image.h>
#include <CGAL/Data_classification/Color.h>

#include <boost/iterator/counting_iterator.hpp>

#define CGAL_CLASSIFICATION_VERBOSE
#if defined(CGAL_CLASSIFICATION_VERBOSE)
#define CGAL_CLASSIFICATION_CERR std::cerr
#else
#define CGAL_CLASSIFICATION_CERR std::ostream(0)
#endif

namespace CGAL {


class Abstract_segmentation_attribute
{
public:

  virtual ~Abstract_segmentation_attribute() { }
  
  virtual double value (int pt_index) = 0;

  virtual double favored (int pt_index) { return (1. - value (pt_index)); }
  virtual double penalized (int pt_index) { return value (pt_index); }
  //  virtual double ignored (int pt_index) { return std::min (favored(pt_index), penalized(pt_index)); }
  virtual double ignored (int) { return 0.5; }

  virtual std::string id() { return "abstract_attribute"; }

  void compute_mean_max (std::vector<double>& vect, double& mean, double& max)
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

};

class Classification_type
{
public:
  enum Attribute_side { FAVORED_ATT = 0,
                        PENALIZED_ATT = 2,
                        NEUTRAL_ATT = 1};

private:
  std::string m_id;
  std::map<std::string, Attribute_side> m_attribute_effects;

public:
  Classification_type (std::string id) : m_id (id) { }

  void change_attribute_effect (Abstract_segmentation_attribute* att, Attribute_side effect)
  {
    m_attribute_effects[att->id()] = effect;
  }

  Attribute_side attribute_effect (Abstract_segmentation_attribute* att) 
  {
    typename std::map<std::string, Attribute_side>::iterator
      search = m_attribute_effects.find (att->id());
    return (search == m_attribute_effects.end () ? NEUTRAL_ATT : search->second);
  }
  
  const std::string& id() const { return m_id; }

  void info()
  {
    std::cerr << "Attribute " << m_id << ": ";
    for (typename std::map<std::string, Attribute_side>::iterator it = m_attribute_effects.begin();
         it != m_attribute_effects.end(); ++ it)
      {
        if (it->second == NEUTRAL_ATT)
          continue;
        
        std::cerr << it->first;
        if (it->second == FAVORED_ATT) std::cerr << " (favored), ";
        else if (it->second == PENALIZED_ATT) std::cerr << " (penalized), ";
      }
    std::cerr << std::endl;
  }

  // Convenience functions
  void make_vegetation ()
  {
    m_attribute_effects.clear();
    m_attribute_effects["scatter"] = FAVORED_ATT;
    m_attribute_effects["distance_to_plane"] = FAVORED_ATT;
    m_attribute_effects["elevation"] = FAVORED_ATT;
    m_id = "vegetation";
  }
  void make_ground ()
  {
    m_attribute_effects.clear();
    m_attribute_effects["scatter"] = PENALIZED_ATT;
    m_attribute_effects["distance_to_plane"] = PENALIZED_ATT;
    m_attribute_effects["horizontality"] = PENALIZED_ATT;
    m_attribute_effects["elevation"] = PENALIZED_ATT;
    m_id = "ground";
  }
  void make_road ()
  {
    m_attribute_effects.clear();
    m_attribute_effects["scatter"] = PENALIZED_ATT;
    m_attribute_effects["distance_to_plane"] = PENALIZED_ATT;
    m_attribute_effects["horizontality"] = PENALIZED_ATT;
    m_attribute_effects["elevation"] = PENALIZED_ATT;
    m_attribute_effects["color"] = FAVORED_ATT;
    m_id = "road";
  }
  void make_roof ()
  {
    m_attribute_effects.clear();
    m_attribute_effects["scatter"] = PENALIZED_ATT;
    m_attribute_effects["horizontality"] = PENALIZED_ATT;
    m_attribute_effects["elevation"] = FAVORED_ATT;
    m_id = "roof";
  }
  void make_facade ()
  {
    m_attribute_effects.clear();
    m_attribute_effects["distance_to_plane"] = PENALIZED_ATT;
    m_attribute_effects["horizontality"] = FAVORED_ATT;
    m_attribute_effects["elevation"] = FAVORED_ATT;
    m_id = "facade";
  }
  void make_building ()
  {
    m_attribute_effects.clear();
    m_attribute_effects["scatter"] = PENALIZED_ATT;
    m_attribute_effects["elevation"] = FAVORED_ATT;
    m_id = "building";
  }
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
  
  typedef typename Kernel::Iso_cuboid_3 Iso_cuboid_3;

  typedef Data_classification::Image<std::vector<int> > Image_indices;
  typedef Data_classification::Image<bool> Image_bool;
  typedef Data_classification::Image<float> Image_float;
  typedef Data_classification::RGB_Color RGB_Color;
  typedef Data_classification::HSV_Color HSV_Color;
  
  struct HPoint {
    Point position;
    unsigned char echo; 
    int ind_x;
    int ind_y;
    std::size_t group; 
    unsigned char AE_label;
    double confidence;
    RGB_Color color;
  };

  class My_point_property_map{
    const std::vector<HPoint>& points;
  public:
    typedef Point value_type;
    typedef const value_type& reference;
    typedef std::size_t key_type;
    typedef boost::lvalue_property_map_tag category;  
    My_point_property_map (const std::vector<HPoint>& pts) : points (pts) {}
    reference operator[] (key_type k) const { return points[k].position; }
    friend inline reference get (const My_point_property_map& ppmap, key_type i) 
    { return ppmap[i]; }
  };
  
  typedef CGAL::Search_traits_3<Kernel> SearchTraits_3;
  typedef Search_traits_adapter <std::size_t, My_point_property_map, SearchTraits_3> Search_traits;
  typedef CGAL::Kd_tree<Search_traits> Tree;
  typedef CGAL::Fuzzy_sphere<Search_traits> Fuzzy_sphere;
  typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> Neighbor_search;


  std::vector<Classification_type*> segmentation_classes;
  std::vector<Abstract_segmentation_attribute*> segmentation_attributes;

  typedef Classification_type::Attribute_side Attribute_side;
  std::vector<std::vector<Attribute_side> > effect_table;

  bool has_colors;
  bool is_echo_given;
  Iso_cuboid_3 BBox_scan;

  std::vector<HPoint> HPS;
  Image_indices grid_HPS;

  //Hpoints attributes
  std::vector<Vector> normal;
  std::vector<Eigenvalues> eigenvalues;
  std::vector<Plane> planes;
  std::vector<std::vector<int> > spherical_neighborhood;
  std::vector<std::vector<int> > cylindrical_neighborhood;

  std::vector<Plane> groups;
  
  Image_bool Mask;                     //imagg
  Image_float DTM; //a enregistrer ?   //imagg           

  double m_grid_resolution;
  double m_radius_neighbors; 
  double m_radius_dtm;
  bool m_multiplicative;

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
        HPS.back().group = (std::size_t)(-1);
        HPS.back().AE_label = (unsigned char)(-1);
        HPS.back().confidence = 0;
        RGB_Color c = {{ 0, 0, 0 }};
        HPS.back().color = c;
      }
    has_colors = false;
    is_echo_given = false;
    m_multiplicative = false;
  }



  void change_hue (RGB_Color& color, const RGB_Color& hue)
  {
    HSV_Color hcolor = Data_classification::rgb_to_hsv (color);
    HSV_Color hhue = Data_classification::rgb_to_hsv (hue);
    hcolor[0] = hhue[0];
    //    hcolor[1] = hhue[1];
    color = Data_classification::hsv_to_rgb (hcolor);
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
    
    My_point_property_map pmap (HPS);
    Tree tree (boost::counting_iterator<std::size_t> (0),
               boost::counting_iterator<std::size_t> (HPS.size()),
               typename Tree::Splitter(),
               Search_traits (pmap));

    std::size_t nb_neigh = 0;
    for(int i=0;i<(int)HPS.size();i++){
      const Point& query=HPS[i].position;
      std::vector<std::size_t> neighbors;
      
      Fuzzy_sphere fs(i, m_radius_neighbors, 0, tree.traits());
      tree.search(std::back_inserter(neighbors), fs);
      nb_neigh += neighbors.size();
      std::vector<Point> neighborhood;
      for (std::size_t j = 0; j < neighbors.size(); ++ j)
        neighborhood.push_back (HPS[neighbors[j]].position);
        
      compute_principal_curvature (query, neighborhood);
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

  double classification_value (std::size_t segmentation_class, int pt_index)
  {
    double out = 0.;
    if (m_multiplicative)
      {
        out = 1.;
        for (std::size_t i = 0; i < effect_table[segmentation_class].size(); ++ i)
          {
            if (effect_table[segmentation_class][i] == Classification_type::FAVORED_ATT)
              out *= segmentation_attributes[i]->favored (pt_index);
            else if (effect_table[segmentation_class][i] == Classification_type::PENALIZED_ATT)
              out *= segmentation_attributes[i]->penalized (pt_index);
            else if (effect_table[segmentation_class][i] == Classification_type::NEUTRAL_ATT)
              out *= segmentation_attributes[i]->ignored (pt_index);
          }
      }
    else
      {
        for (std::size_t i = 0; i < effect_table[segmentation_class].size(); ++ i)
          {
            if (effect_table[segmentation_class][i] == Classification_type::FAVORED_ATT)
              out += segmentation_attributes[i]->favored (pt_index);
            else if (effect_table[segmentation_class][i] == Classification_type::PENALIZED_ATT)
              out += segmentation_attributes[i]->penalized (pt_index);
            else if (effect_table[segmentation_class][i] == Classification_type::NEUTRAL_ATT)
              out += segmentation_attributes[i]->ignored (pt_index);
          }
      }
    return out;
  }


  bool quick_labeling_PC()
  {

    // data term initialisation
    CGAL_CLASSIFICATION_CERR << "Labeling... ";


    int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
    for(int s=0;s<(int)HPS.size();s++){
			
      int nb_class_best=0; 

      double val_class_best = std::numeric_limits<double>::max();
      std::vector<double> values;
      
      for(std::size_t k = 0; k < effect_table.size(); ++ k)
        {
          double value = classification_value (k, s);
          values.push_back (value);
          
          if(val_class_best > value)
            {
              val_class_best = value;
              nb_class_best=k;
            }
        }

      HPS[s].AE_label = nb_class_best;

      std::sort (values.begin(), values.end());
      if (m_multiplicative)
        HPS[s].confidence = (values[1] - values[0]) / values[1];
      else
        HPS[s].confidence = values[1] - values[0];

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

    std::vector<std::vector<double> > values
      (segmentation_classes.size(),
       std::vector<double> (HPS.size(), -1.));

    My_point_property_map pmap (HPS);

    Tree tree (boost::counting_iterator<std::size_t> (0),
               boost::counting_iterator<std::size_t> (HPS.size()),
               typename Tree::Splitter(),
               Search_traits (pmap));

    for (std::size_t s=0; s < HPS.size(); ++ s)
      {
        std::vector<std::size_t> neighbors;
      
        Fuzzy_sphere fs(s, m_radius_neighbors, 0, tree.traits());
        tree.search(std::back_inserter(neighbors), fs);

        std::vector<double> mean (values.size(), 0.);
        for (std::size_t n = 0; n < neighbors.size(); ++ n)
          {
            if (values[0][neighbors[n]] < 0.)
              for(std::size_t k = 0; k < effect_table.size(); ++ k)
                {
                  values[k][neighbors[n]] = classification_value (k, neighbors[n]);
                  mean[k] += values[k][neighbors[n]];
                }
            else
              for (std::size_t j = 0; j < values.size(); ++ j)
                mean[j] += values[j][neighbors[n]];
          }

        int nb_class_best=0; 
        double val_class_best = std::numeric_limits<double>::max();
        for(std::size_t k = 0; k < mean.size(); ++ k)
          {
            mean[k] /= neighbors.size();
            if(val_class_best > mean[k])
              {
                val_class_best = mean[k];
                nb_class_best = k;
              }
          }

        HPS[s].AE_label = nb_class_best;

        std::sort (mean.begin(), mean.end());
        HPS[s].confidence = mean[1] - mean[0];      

      }

	
    return true;
  }

  void reset_groups()
  {
    groups.clear();
    for (std::size_t i = 0; i < HPS.size(); ++ i)
      HPS[i].group = (std::size_t)(-1);
  }

  bool regularized_labeling_PC()
  {
    std::vector<std::vector<std::size_t> > groups;
    for (std::size_t i = 0; i < HPS.size(); ++ i)
      {
        std::size_t index = HPS[i].group;
        if (index == (std::size_t)(-1))
          continue;

        if (groups.size() <= index)
          groups.resize (index + 1);
        
        groups[index].push_back (i);
      }

    if (groups.empty())
      return false;
        
    // data term initialisation
    CGAL_CLASSIFICATION_CERR << "Labeling... ";

    std::vector<std::vector<double> > values;
    values.resize (segmentation_classes.size());
    
    std::map<Point, std::size_t> map_p2i;
    for(int s=0;s<(int)HPS.size();s++)
      {
        map_p2i[HPS[s].position] = s;

        int nb_class_best=0; 
        double val_class_best = std::numeric_limits<double>::max();
        for(std::size_t k = 0; k < effect_table.size(); ++ k)
          {
            double v = classification_value (k, s);

            values[k].push_back(v);
            if (v < val_class_best)
              {
                nb_class_best = k;
                val_class_best = v;
              }
          }
        HPS[s].AE_label = nb_class_best;
      }

    for(std::size_t i = 0; i < groups.size(); ++ i)
      {
        std::vector<double> mean (values.size(), 0.);

        for (std::size_t n = 0; n < groups[i].size(); ++ n)
          {
            for (std::size_t j = 0; j < values.size(); ++ j)
              mean[j] += values[j][groups[i][n]];
          }
        
        int nb_class_best=0; 

        double val_class_best = std::numeric_limits<double>::max();

        for (std::size_t j = 0; j < mean.size(); ++ j)
          if (mean[j] < val_class_best)
            {
              nb_class_best = j;
              val_class_best = mean[j];
            }

        for (std::size_t n = 0; n < groups[i].size(); ++ n)
          {
            HPS[groups[i][n]].AE_label = nb_class_best;
            for (std::size_t j = 0; j < mean.size(); ++ j)
              values[j][groups[i][n]] = mean[j] / groups[i].size();
          }
      }

    My_point_property_map pmap (HPS);

    Tree tree (boost::counting_iterator<std::size_t> (0),
               boost::counting_iterator<std::size_t> (HPS.size()),
               typename Tree::Splitter(),
               Search_traits (pmap));

    for (std::size_t s=0; s < HPS.size(); ++ s)
      {
        if (HPS[s].group != (std::size_t)(-1))
          continue;
        
        std::vector<std::size_t> neighbors;
      
        Fuzzy_sphere fs(s, m_radius_neighbors, 0, tree.traits());
        tree.search(std::back_inserter(neighbors), fs);

        std::vector<double> mean (values.size(), 0.);
        for (std::size_t n = 0; n < neighbors.size(); ++ n)
          for (std::size_t j = 0; j < values.size(); ++ j)
            mean[j] += values[j][neighbors[n]];

        int nb_class_best=0; 
        double val_class_best = std::numeric_limits<double>::max();
        for(std::size_t k = 0; k < mean.size(); ++ k)
          {
            mean[k] /= neighbors.size();
            if(val_class_best > mean[k])
              {
                val_class_best = mean[k];
                nb_class_best = k;
              }
          }

        HPS[s].AE_label = nb_class_best;

        std::sort (mean.begin(), mean.end());
        HPS[s].confidence = mean[1] - mean[0];      

      }
    
	
    return true;
  }

  bool point_cloud_classification(int method)
  {

    clock_t t;
    t = clock();

    build_effect_table ();
	
    CGAL_CLASSIFICATION_CERR<<std::endl<<"Classification of the point cloud: ";
    
    if (method == 0)
      quick_labeling_PC();
    else if (method == 1)
      smoothed_labeling_PC();
    else if (method == 2)
      regularized_labeling_PC();
    else
      {
        std::cerr << "Unknown method number." << std::endl;
        abort();
      }

    CGAL_CLASSIFICATION_CERR<<"-> OK ( "<<((float)clock()-t)/CLOCKS_PER_SEC<<" sec )"<< std::endl;

    return true;
  }


  void build_effect_table ()
  {
    effect_table = std::vector<std::vector<Attribute_side> >
      (segmentation_classes.size(), std::vector<Attribute_side> (segmentation_attributes.size(),
                                                                 Classification_type::NEUTRAL_ATT));
    
    for (std::size_t i = 0; i < effect_table.size (); ++ i)
      for (std::size_t j = 0; j < effect_table[i].size (); ++ j)
        effect_table[i][j] = segmentation_classes[i]->attribute_effect (segmentation_attributes[j]);

  }

protected: 



};


template <typename Kernel>
class Scatter_segmentation_attribute : public Abstract_segmentation_attribute
{
  typedef Point_set_classification<Kernel> PSC;
  typedef typename PSC::Image_float Image_float;
  typedef typename Data_classification::RGB_Color RGB_Color;
  typedef typename Data_classification::HSV_Color HSV_Color;
  
  std::vector<double> vegetation_attribute;
  
public:
  double weight;
  double mean;
  double max;

  Scatter_segmentation_attribute (PSC& M,
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
    else //if (!(M.has_colors))
      {
        std::cerr << "No colors" << std::endl;
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
					
            //              Vegetation(i,j)= (float)(1. - (1. / (float)nb_layers));
            Vegetation(i,j)= 1.f - (nb_occ / (float)(occupy.size()));
			
          }
		
        }
        for(int i=0;i<(int)M.HPS.size();i++){
          int I= M.HPS[i].ind_x;
          int J= M.HPS[i].ind_y;
          vegetation_attribute.push_back((double)Vegetation(I,J));
        }
	
        // for(int i=0;i<(int)M.HPS.size();i++){
        //   if (M.eigenvalues[i][2] == 0)
        //     vegetation_attribute.push_back(0);
        //   else
        //     vegetation_attribute.push_back(std::fabs(M.eigenvalues[i][0] / M.eigenvalues[i][2]));
        // }
    }
    // else
    //   {
    //     std::cerr << "Using colors" << std::endl;
    //     for(int i=0;i<(int)M.HPS.size();i++)
    //       {
    //         HSV_Color c = internal::Classification::rgb_to_hsv<HSV_Color> (M.HPS[i].color);
    //         vegetation_attribute.push_back (std::exp (-(c[0] - 120) * (c[0] - 120) / (2 * 60 * 60))
    //                                         * std::exp (-(c[1] - 25) * (c[1] - 25) / (2 * 20 * 20))
    //                                         * std::exp (-(c[2] - 37) * (c[2] - 37) / (2 * 10 * 10)));
    //       }
    //   }

    this->compute_mean_max (vegetation_attribute, mean, max);
    //    max *= 2;
  }

  virtual double value (int pt_index)
  {
    return std::max (0., std::min (1., vegetation_attribute[pt_index] / weight));
  }

  virtual std::string id() { return "scatter"; }
};




template <typename Kernel>
class Color_segmentation_attribute : public Abstract_segmentation_attribute
{
  typedef Point_set_classification<Kernel> PSC;
  typedef typename PSC::Image_float Image_float;
  typedef typename Data_classification::RGB_Color RGB_Color;
  typedef typename Data_classification::HSV_Color HSV_Color;
  
  std::vector<double> color_attribute;
  
public:
  double weight;
  double mean;
  double max;

  Color_segmentation_attribute (PSC& M,
                                double weight) : weight (weight)
  {

    std::cerr << "Using colors" << std::endl;
    for(int i=0;i<(int)M.HPS.size();i++)
      {
        HSV_Color c = Data_classification::rgb_to_hsv (M.HPS[i].color);
        color_attribute.push_back (std::exp (-(c[0] - 156.) * (c[0] - 156.) / (2. * 81. * 81.))
                                   * std::exp (-(c[1] - 5.) * (c[1] - 5.) / (2. * 4. * 4.))
                                   * std::exp (-(c[2] - 76.) * (c[2] - 76.) / (2. * 6.8 * 6.8)));
      }
    this->compute_mean_max (color_attribute, mean, max);
  }

  virtual double value (int pt_index)
  {
    return std::max (0., std::min (1., color_attribute[pt_index] / weight));
  }

  virtual std::string id() { return "color"; }
};


template <typename Kernel>
class Distance_to_plane_segmentation_attribute : public Abstract_segmentation_attribute
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

    this->compute_mean_max (distance_to_plane_attribute, mean, max);
    //    max *= 2;
  }

  virtual double value (int pt_index)
  {
    return std::max (0., std::min (1., distance_to_plane_attribute[pt_index] / weight));
  }

  virtual std::string id() { return "distance_to_plane"; }
};

template <typename Kernel>
class Horizontality_segmentation_attribute : public Abstract_segmentation_attribute
{
  typedef Point_set_classification<Kernel> PSC;

  std::vector<double> horizontality_attribute;
  
public:
  double weight;
  double mean;
  double max;
  
  Horizontality_segmentation_attribute (PSC& M, double weight) : weight (weight)
  {
    typename Kernel::Vector_3 vertical (0., 0., 1.);
    
    for(int i=0; i<(int)M.HPS.size(); i++)
      {
        typename Kernel::Vector_3 normal = M.planes[i].orthogonal_vector();
        normal = normal / std::sqrt (normal * normal);
        horizontality_attribute.push_back (1. - std::fabs(normal * vertical));
      }
    
    this->compute_mean_max (horizontality_attribute, mean, max);
    //    max *= 2;
  }

  virtual double value (int pt_index)
  {
    return std::max (0., std::min (1., horizontality_attribute[pt_index] / weight));
  }

  virtual std::string id() { return "horizontality"; }
};

template <typename Kernel>
class Elevation_segmentation_attribute : public Abstract_segmentation_attribute
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

    this->compute_mean_max (elevation_attribute, mean, max);
    //    max *= 5;
  }

  virtual double value (int pt_index)
  {
    return std::max (0., std::min (1., elevation_attribute[pt_index] / weight));
  }
  
  virtual std::string id() { return "elevation"; }

};
  



} // namespace CGAL

#endif // CGAL_POINT_SET_CLASSIFICATION_H

