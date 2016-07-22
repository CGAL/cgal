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
#include <queue>

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
#include <CGAL/Data_classification/Segmentation_attribute.h>
#include <CGAL/gco/GCoptimization.h>
#include <boost/iterator/counting_iterator.hpp>

#define CGAL_CLASSIFICATION_VERBOSE
#if defined(CGAL_CLASSIFICATION_VERBOSE)
#define CGAL_CLASSIFICATION_CERR std::cerr
#else
#define CGAL_CLASSIFICATION_CERR std::ostream(0)
#endif

namespace CGAL {



/*!
\ingroup PkgDataClassification

\brief Definition of a classification type based on an ID and a set of
relationship with attributes.

A classification type is used to segment the input data set. Usual
classification types are ground, vegetation and buildings (but other
can be defined).

*/
class Classification_type
{
public:
  
  enum Attribute_side /// Defines the effect of the values of an attribute on the classification type.
    {
      FAVORED_ATT = 0, ///< High values of the attribute favor this type
      NEUTRAL_ATT = 1, ///< The attribute has no effect on this type
      PENALIZED_ATT = 2 ///< Low values of the attribute favor this type
    };

private:
  std::string m_id;
  std::map<Segmentation_attribute*, Attribute_side> m_attribute_effects;

public:

  /// \name Main methods 
  /// @{
  /*! 
    \param id The name of the classification type
    (e.g. vegetation). Two different classification types must have
    different IDs.
  */ 
  Classification_type (std::string id) : m_id (id) { }

  /*! 
    \brief Sets how an attribute affects the classification type.

    \param att Attribute whose effect on the classification type will be set

    \param effect The effect the attribute will have on the classification type

  */ 
  void set_attribute_effect (Segmentation_attribute* att, Attribute_side effect)
  {
    m_attribute_effects[att] = effect;
  }

  /*!
    \brief Get the effects of an attribute on the classification type.

    \param att Attribute

    \return The effect of the attribute on the classification type.
   */
  Attribute_side attribute_effect (Segmentation_attribute* att) 
  {
    std::map<Segmentation_attribute*, Attribute_side>::iterator
      search = m_attribute_effects.find (att);
    return (search == m_attribute_effects.end () ? NEUTRAL_ATT : search->second);
  }

  /*!
    \brief Get the ID of the classification type.

    \return The ID of the classification type.
  */
  const std::string& id() const { return m_id; }

  /// @}
  
  /// \cond SKIP_IN_MANUAL
  void info()
  {
    std::cerr << "Attribute " << m_id << ": ";
    for (std::map<Segmentation_attribute*, Attribute_side>::iterator it = m_attribute_effects.begin();
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
  /// \endcond

};



/*!
\ingroup PkgDataClassification

\brief Classifies a point set based on a set of attribute and a set of classification types.

This class implement the core of the algorithm. It uses a point set as
input. Based on a set of segmentation attributes and a set of
classification types, it segments the point set into the different
types given. The output can be regularized with different smoothing
methods.

\tparam Kernel The geometric kernel used

*/
template <typename Kernel>
class Point_set_classification
{

  
public:
  /// \cond SKIP_IN_MANUAL
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
    unsigned char neighbor;
    double confidence;
    RGB_Color color;
  };

  struct Cluster
  {
    Point centroid;
    std::vector<std::size_t> indices;
    std::set<std::size_t> neighbors;
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
  typedef typename Neighbor_search::Tree KTree;
  typedef typename Neighbor_search::Distance Distance;

  std::vector<Classification_type*> segmentation_classes; 
  std::vector<Segmentation_attribute*> segmentation_attributes; 

  typedef Classification_type::Attribute_side Attribute_side;
  std::vector<std::vector<Attribute_side> > effect_table;

  bool has_colors;
  bool is_echo_given;
  Iso_cuboid_3 BBox_scan;

  std::vector<HPoint> HPS;
  Image_indices grid_HPS;

  //Hpoints attributes
  std::vector<Eigenvalues> eigenvalues;
  std::vector<Plane> planes;

  std::vector<Plane> groups;
  std::vector<Cluster> clusters;
  
  Image_bool Mask;                     //imagg
  Image_float DTM; //a enregistrer ?   //imagg           

  double m_grid_resolution;
  double m_radius_neighbors; 
  double m_radius_dtm;
  bool m_multiplicative;
  /// \endcond

  /// \name Main methods
  /// @{
  
  /*! 
    \brief Constructs a classification object based on the input range.

    \tparam InputIterator Iterator on the input point. Value type must be `Point_3<Kernel>`.

    \param begin Iterator to the first input point

    \param end Past-the-end iterator

    \param grid_resolution Resolution of the 2D map of the ground. If
    the default value is used, it is computed as the average spacing
    of the point set.

    \param radius_neighbors Size used for neighborhood computation. If
    the default value is used, it is computed as 5 times
    `grid_resolution`.

    \param radius_dtm Size used to estimate the ground (should be
    higher than the thickest non-ground object of the scene). If the
    default value is used, it is computed as 5 times `radius_neighbors`.
  */ 
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
        HPS.back().neighbor = (unsigned char)(-1);
        HPS.back().confidence = 0;
        RGB_Color c = {{ 0, 0, 0 }};
        HPS.back().color = c;
      }
    has_colors = false;
    is_echo_given = false;
    m_multiplicative = false;
  }

  /// \cond SKIP_IN_MANUAL
  void change_hue (RGB_Color& color, const RGB_Color& hue)
  {
    HSV_Color hcolor = Data_classification::rgb_to_hsv (color);
    HSV_Color hhue = Data_classification::rgb_to_hsv (hue);
    hcolor[0] = hhue[0];
    //    hcolor[1] = hhue[1];
    color = Data_classification::hsv_to_rgb (hcolor);
  }
  /// \endcond

  /*!
    \brief Change the parameters of the classification algorithm.

    \param grid_resolution Resolution of the 2D map of the ground. 

    \param radius_neighbors Size used for neighborhood computation. 

    \param radius_dtm Size used to estimate the ground (should be
    higher than the thickest non-ground object of the scene). 
  */
  
  void set_parameters (double grid_resolution, double radius_neighbors, double radius_dtm)
  {
    m_grid_resolution = grid_resolution;
    m_radius_neighbors = radius_neighbors;
    m_radius_dtm = radius_dtm;
  }

  /*!
    Computes all the interlate structures needed by the algorithm.
   */
  void initialization()
  {
    clock_t t;
    t = clock();

    CGAL_CLASSIFICATION_CERR << std::endl << "Initialization: ";

    //1-Neighborhood computation and reset the attributes of the structure points
    CGAL_CLASSIFICATION_CERR<<"spherical neighborhood..";
    eigenvalues.clear();
    planes.clear();
    
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

    CGAL_CLASSIFICATION_CERR<<std::endl<<"-> OK ( "<<((float)clock()-t)/CLOCKS_PER_SEC<<" sec )"<< std::endl;
  }


  /// \cond SKIP_IN_MANUAL
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

  double classification_value (std::size_t segmentation_class, int pt_index) const
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


  class Edge_score
  {
    const Point_set_classification& M;

  public:
    Edge_score(const Point_set_classification& _M) : M(_M) {}

    float compute(int, int, int l1, int l2)
    {
      double res=0;
      double smooth_seg= M.m_grid_resolution;

      if(l1!=l2) res=1; 

      return smooth_seg*res;
    }
  };


  class Facet_score
  {
    const Point_set_classification& M;

  public:
    Facet_score(const Point_set_classification& _M) : M(_M) {}

    double compute(int s, int l)
    {
      return M.classification_value (l, s);
    }
  };
  
  bool graphcut_labeling_PC()
  {

    std::size_t nb_alpha_exp = 2;
    
    // data term initialisation
    CGAL_CLASSIFICATION_CERR << "Labeling... ";

    std::vector<std::vector<double> > values
      (segmentation_classes.size(),
       std::vector<double> (HPS.size(), -1.));

    My_point_property_map pmap (HPS);

    KTree tree (boost::counting_iterator<std::size_t> (0),
                boost::counting_iterator<std::size_t> (HPS.size()),
                typename KTree::Splitter(),
                Search_traits (pmap));

    GCoptimizationGeneralGraph<float, float> *gc= new GCoptimizationGeneralGraph<float, float>
      ((int)HPS.size(),(int)(segmentation_classes.size()));

    gc->specializeDataCostFunctor(Facet_score(*this));
    gc->specializeSmoothCostFunctor(Edge_score(*this));

    Distance tr_dist(pmap);
    for (std::size_t s=0; s < HPS.size(); ++ s)
      {
        std::vector<std::size_t> neighbors;

        Neighbor_search search (tree, HPS[s].position, 12, 0, true, tr_dist);
        for (typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++ it)
          gc->setNeighbors (s, it->first);
        
        int nb_class_best=0; 
        double val_class_best = std::numeric_limits<double>::max();
        for(std::size_t k = 0; k < effect_table.size(); ++ k)
          {
            double value = classification_value (k, s);

            if(val_class_best > value)
              {
                val_class_best = value;
                nb_class_best=k;
              }
          }
        gc->setLabel (s, nb_class_best);

      }

    CGAL_CLASSIFICATION_CERR << "Graph cut... ";
    for (std::size_t iter = 0; iter < nb_alpha_exp; ++ iter)
      {
        for (std::size_t i = 0; i< segmentation_classes.size(); ++ i)
          {
				
            // Compute vector of active sites
            std::vector<int> sites;
            for (std::size_t s = 0; s < HPS.size(); ++ s)
              sites.push_back ((int)s);
            
            // Compute alpha expansion for this label on these sites
            gc->alpha_expansion((int)i, &(sites[0]), sites.size());
          }
      }
    
    for (std::size_t s=0; s < HPS.size(); ++ s)
      HPS[s].AE_label = gc->whatLabel (s);

    delete gc;
      
    return true;
  }
  
  void reset_groups()
  {
    groups.clear();
    for (std::size_t i = 0; i < HPS.size(); ++ i)
      HPS[i].group = (std::size_t)(-1);
  }

  void cluster_points (const double& tolerance)
  {
    if (clusters.empty())
      {
        My_point_property_map pmap (HPS);

        Tree tree (boost::counting_iterator<std::size_t> (0),
                   boost::counting_iterator<std::size_t> (HPS.size()),
                   typename Tree::Splitter(),
                   Search_traits (pmap));

        std::vector<std::size_t> done (HPS.size(), (std::size_t)(-1));

        for (std::size_t s=0; s < HPS.size(); ++ s)
          HPS[s].neighbor = (unsigned char)(-1);
    
        for (std::size_t s=0; s < HPS.size(); ++ s)
          {
            if (done[s] != (std::size_t)(-1))
              continue;
            unsigned char label = HPS[s].AE_label;
        
            clusters.push_back (Cluster());

            std::queue<std::size_t> todo;
            todo.push (s);
            done[s] = clusters.size()-1;

            while (!(todo.empty()))
              {
                std::size_t current = todo.front();
                todo.pop();
                clusters.back().indices.push_back (current);
            
                std::vector<std::size_t> neighbors;
                Fuzzy_sphere fs(current, tolerance, 0, tree.traits());
                tree.search(std::back_inserter(neighbors), fs);
            
                for (std::size_t n = 0; n < neighbors.size(); ++ n)
                  {
                    if (done[neighbors[n]] == (std::size_t)(-1))
                      {
                        if (HPS[neighbors[n]].AE_label == label)
                          {
                            todo.push (neighbors[n]);
                            done[neighbors[n]] = clusters.size()-1;
                          }
                        else
                          {
                            HPS[current].neighbor = HPS[neighbors[n]].AE_label;
                            HPS[neighbors[n]].neighbor = HPS[current].AE_label;
                            continue;
                          }
                      }
                    else if (done[neighbors[n]] != clusters.size()-1)
                      {
                        clusters.back().neighbors.insert (done[neighbors[n]]);
                        clusters[done[neighbors[n]]].neighbors.insert (clusters.size()-1);
                      }
                  }
              }
          }
        std::cerr << "Found " << clusters.size() << " cluster(s)" << std::endl;

        for (std::size_t i = 0; i < clusters.size(); ++ i)
          {
            std::vector<Point> pts;
            for (std::size_t j = 0; j < clusters[i].indices.size(); ++ j)
              pts.push_back (HPS[clusters[i].indices[j]].position);
            clusters[i].centroid = CGAL::centroid (pts.begin(), pts.end());
        
          }
      }
    else
      {
        std::size_t index_ground = 0;
        std::size_t index_roof = 0;
        std::size_t index_facade = 0;
        std::size_t index_vege = 0;
  
        for (std::size_t i = 0; i < segmentation_classes.size(); ++ i)
          if (segmentation_classes[i]->id() == "ground")
            index_ground = i;
          else if (segmentation_classes[i]->id() == "roof")
            index_roof = i;
          else if (segmentation_classes[i]->id() == "facade")
            index_facade = i;
          else if (segmentation_classes[i]->id() == "vegetation")
            index_vege = i;

        std::size_t min_ind = 0;
        double min = 1.;
        
        for (std::size_t i = 0; i < clusters.size(); ++ i)
          {
            std::size_t nb_border = 0;
            std::size_t nb_good = 0;

            for (std::size_t n = 0; n < clusters[i].indices.size(); ++ n)
              {
                if (HPS[clusters[i].indices[n]].neighbor == (unsigned char)(-1))
                  continue;

                std::size_t c0 = HPS[clusters[i].indices[n]].AE_label;
                std::size_t c1 = HPS[clusters[i].indices[n]].neighbor;
                nb_border ++;
                if ((c0 == index_ground && c1 == index_facade) || (c0 == index_facade && c1 == index_ground) ||
                    (c0 == index_ground && c1 == index_vege) || (c0 == index_vege && c1 == index_ground) ||
                    (c0 == index_facade && c1 == index_roof) || (c0 == index_roof && c1 == index_facade))
                  nb_good ++;
              }

            if (nb_border == 0)
              continue;

            double ratio = nb_good / (double)nb_border;
            if (ratio < min)
              {
                min = ratio;
                min_ind = i;
              }
          }

        std::cerr << "Found cluster " << min_ind << " with only " << (int)(100. * min) << "% of correct border." << std::endl;

        std::vector<double> values (segmentation_classes.size(), 0.);
        
        for (std::size_t n = 0; n < clusters[min_ind].indices.size(); ++ n)
          {
            for (std::size_t i = 0; i < values.size(); ++ i)
              values[i] += classification_value (i, clusters[min_ind].indices[n]);
          }

        std::size_t nb_class_best = 0;
        double val_class_best = std::numeric_limits<double>::max();
        
        for (std::size_t i = 0; i < values.size(); ++ i)
          {
            if (i == HPS[clusters[min_ind].indices[0]].AE_label)
              continue;

            if (values[i] < val_class_best)
              {
                nb_class_best = i;
                val_class_best = values[i];
              }
          }

        for (std::size_t n = 0; n < clusters[min_ind].indices.size(); ++ n)
          HPS[clusters[min_ind].indices[n]].AE_label = nb_class_best;
        
        clusters.clear();        
      }

  }

  /// \cond SKIP_IN_MANUAL
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
  /// \endcond

  /*!
    Performs classification using the user-chosen regularization method.

    \param method Regularization method. 0 = no regularization; 1 =
    global regularization through graphcut; 2 = local regularization
    based on pre-computed groups.

    \note Classification without regularization (method = 0) is very
    quick: almost instantaneous or up to a few seconds for very large
    point clouds. Regularization improves the quality of the output at
    the cost of longer computation time.
  */
  void classify (int method)
  {

    clock_t t;
    t = clock();

    build_effect_table ();
	
    CGAL_CLASSIFICATION_CERR<<std::endl<<"Classification of the point cloud: ";
    
    if (method == 0)
      quick_labeling_PC();
    else if (method == 1)
      //      smoothed_labeling_PC();
      graphcut_labeling_PC();
    else if (method == 2)
      regularized_labeling_PC();
    else
      {
        std::cerr << "Unknown method number." << std::endl;
        abort();
      }

    CGAL_CLASSIFICATION_CERR<<"-> OK ( "<<((float)clock()-t)/CLOCKS_PER_SEC<<" sec )"<< std::endl;

  }

  /// @}

  
  /// \cond SKIP_IN_MANUAL
  void build_effect_table ()
  {
    effect_table = std::vector<std::vector<Attribute_side> >
      (segmentation_classes.size(), std::vector<Attribute_side> (segmentation_attributes.size(),
                                                                 Classification_type::NEUTRAL_ATT));
    
    for (std::size_t i = 0; i < effect_table.size (); ++ i)
      for (std::size_t j = 0; j < effect_table[i].size (); ++ j)
        effect_table[i][j] = segmentation_classes[i]->attribute_effect (segmentation_attributes[j]);

  }
  /// \endcond



#ifdef DOXYGEN_RUNNING

  /// \name Types and attributes
  /// @{
  
  /*!
    \brief Add a classification type
    \param type Pointer to the classification type object
   */
  void add_classification_type (Classification_type* type);

  void clear_classification_types ();

  /*!
    \brief Add a segmentation attribute
    \param attribute Pointer to the attribute object
   */
  void add_segmentation_attribute (Segmentation_attribute* attribute);

  void clear_segmentation_attributes ();
  
  /// @}

  /// \name Groups
  /// @{

  /*!
    \brief Add a point to a group

    Grouping points can be used for regularization (for example, to
    apply the same classification type to all inliers of a detected
    RANSAC primitive).

    \param point_index Index of the point in the input range
    \param group_index Index of the group
   */
  void add_to_group (std::size_t point_index, std::size_t group_index);

  /*!
    \brief Reset all groups attributes of points
   */
  void clear_groups();

  /// @}
  

  /// \name Output
  /// @{

  /*!
    \brief Get the classification type of indexed point.
    \param index Index of the input point
    \return Pointer to the classification type 
  */
  Classification_type* classification_type_of (std::size_t index) const;

  /// @}

#endif


};





} // namespace CGAL

#endif // CGAL_POINT_SET_CLASSIFICATION_H

