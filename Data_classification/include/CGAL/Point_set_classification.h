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
#include <CGAL/centroid.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <CGAL/Data_classification/Planimetric_grid.h>
#include <CGAL/Data_classification/Neighborhood.h>
#include <CGAL/Data_classification/Image.h>
#include <CGAL/Data_classification/Color.h>
#include <CGAL/Data_classification/Segmentation_attribute.h>
#include <CGAL/gco/GCoptimization.h>

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
template <typename Kernel,
          typename RandomAccessIterator,
          typename PointPMap>
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

  typedef Data_classification::Planimetric_grid<Kernel,
                                                RandomAccessIterator,
                                                PointPMap> Grid;
  typedef Data_classification::Neighborhood<Kernel,
                                            RandomAccessIterator,
                                            PointPMap> Neighborhood;
  typedef Data_classification::Image<float> Image_float;
  typedef Data_classification::RGB_Color RGB_Color;
  typedef Data_classification::HSV_Color HSV_Color;

  class Point_range
  {
    std::size_t m_size;
    RandomAccessIterator m_begin;
    PointPMap m_point_pmap;
    
  public:

    Point_range (RandomAccessIterator begin, RandomAccessIterator end,
                 PointPMap point_pmap)
      : m_size (end - begin), m_begin (begin), m_point_pmap (point_pmap)
    { }

    std::size_t size() const { return m_size; }
    
    const Point& operator[] (std::size_t index) const { return get (m_point_pmap, *(m_begin + index)); }

    RandomAccessIterator begin() { return m_begin; }
    RandomAccessIterator end() { return m_begin + m_size; }

    friend RandomAccessIterator operator+ (Point_range& range, std::size_t index)
    {
      return range.begin() + index;
    }
  };
  

  Point_range m_input;
  
  std::vector<unsigned char> echo;
  std::vector<std::size_t> m_group;
  std::vector<unsigned char> m_assigned_type;
  std::vector<unsigned char> m_neighbor;
  std::vector<double> m_confidence;
  std::vector<RGB_Color> color;

  struct Cluster
  {
    Point centroid;
    std::vector<std::size_t> indices;
    std::set<std::size_t> neighbors;
  };

  
  std::vector<Classification_type*> segmentation_classes; 
  std::vector<Segmentation_attribute*> segmentation_attributes; 

  typedef Classification_type::Attribute_side Attribute_side;
  std::vector<std::vector<Attribute_side> > effect_table;

  Grid grid;

  //Hpoints attributes
  std::vector<Plane> groups;
  std::vector<Cluster> clusters;
  
  double m_grid_resolution;
  double m_radius_neighbors; 
  bool m_multiplicative;
  /// \endcond

  /// \name Main methods
  /// @{
  
  /*! 
    \brief Constructs a classification object based on the input range.

    \param begin Iterator to the first input point

    \param end Past-the-end iterator

    \param grid_resolution Resolution of the 2D map of the ground. If
    the default value is used, it is computed as the average spacing
    of the point set.

    \param radius_neighbors Size used for neighborhood computation. If
    the default value is used, it is computed as 5 times
    `grid_resolution`.

  */ 
  Point_set_classification (RandomAccessIterator begin,
                            RandomAccessIterator end,
                            PointPMap point_pmap,
                            double grid_resolution = -1.,
                            double radius_neighbors = -1.)
    : m_input (begin, end, point_pmap),
      m_grid_resolution (grid_resolution),
      m_radius_neighbors (radius_neighbors)
  {
    if (m_grid_resolution < 0.)
      m_grid_resolution = CGAL::compute_average_spacing<CGAL::Sequential_tag> (begin, end, 6);
    if (m_radius_neighbors < 0.)
      m_radius_neighbors = 5. * m_grid_resolution;

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

  */
  
  void set_parameters (double grid_resolution, double radius_neighbors)
  {
    m_grid_resolution = grid_resolution;
    m_radius_neighbors = radius_neighbors;
  }

  /// \cond SKIP_IN_MANUAL

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
    for (std::size_t s = 0; s < m_input.size(); s++)
      {
			
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

        m_assigned_type[s] = nb_class_best;

        std::sort (values.begin(), values.end());
        if (m_multiplicative)
          m_confidence[s] = (values[1] - values[0]) / values[1];
        else
          m_confidence[s] = values[1] - values[0];

        if(nb_class_best==0) count1++;
        else if(nb_class_best==1) count2++;
        else if(nb_class_best==2) count3++;
        else count4++;

      }
    
    CGAL_CLASSIFICATION_CERR<<"ok"<<std::endl;
	
    return true;
  }


  bool smoothed_labeling_PC (const Neighborhood& neighborhood)
  {

    // data term initialisation
    CGAL_CLASSIFICATION_CERR << "Labeling... ";

    std::vector<std::vector<double> > values
      (segmentation_classes.size(),
       std::vector<double> (m_input.size(), -1.));

    for (std::size_t s=0; s < m_input.size(); ++ s)
      {
        std::vector<std::size_t> neighbors;
        neighborhood.get (s, neighbors);

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

        m_assigned_type[s] = nb_class_best;

        std::sort (mean.begin(), mean.end());
        m_confidence[s] = mean[1] - mean[0];      

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
  
  bool graphcut_labeling_PC(const Neighborhood& neighborhood)
  {

    std::size_t nb_alpha_exp = 2;
    
    // data term initialisation
    CGAL_CLASSIFICATION_CERR << "Labeling... ";

    std::vector<std::vector<double> > values
      (segmentation_classes.size(),
       std::vector<double> (m_input.size(), -1.));

    GCoptimizationGeneralGraph<float, float> *gc= new GCoptimizationGeneralGraph<float, float>
      ((int)(m_input.size()),(int)(segmentation_classes.size()));

    gc->specializeDataCostFunctor(Facet_score(*this));
    gc->specializeSmoothCostFunctor(Edge_score(*this));

    for (std::size_t s=0; s < m_input.size(); ++ s)
      {
        std::vector<std::size_t> neighbors;

        neighborhood.k_neighbors (s, 12, std::back_inserter (neighbors));

        for (std::size_t i = 0; i < neighbors.size(); ++ i)
          gc->setNeighbors (s, neighbors[i]);
        
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
            for (std::size_t s = 0; s < m_input.size(); ++ s)
              sites.push_back ((int)s);
            
            // Compute alpha expansion for this label on these sites
            gc->alpha_expansion((int)i, &(sites[0]), sites.size());
          }
      }
    
    for (std::size_t s=0; s < m_input.size(); ++ s)
      m_assigned_type[s] = gc->whatLabel (s);

    delete gc;
      
    return true;
  }
  
  void reset_groups()
  {
    groups.clear();
    std::vector<std::size_t>(m_input.size(), (std::size_t)(-1)).swap (m_group);
  }

  void cluster_points (const Neighborhood& neighborhood, const double& tolerance)
  {
    std::vector<unsigned char>(m_input.size(), (unsigned char)(-1)).swap (m_neighbor);
    
    std::vector<std::size_t> done (m_input.size(), (std::size_t)(-1));
    
    for (std::size_t s=0; s < m_input.size(); ++ s)
      {
        if (done[s] != (std::size_t)(-1))
          continue;
        unsigned char label = m_assigned_type[s];
        
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
            neighborhood.range_neighbors (current, tolerance,
                                          std::back_inserter (neighbors));

            for (std::size_t n = 0; n < neighbors.size(); ++ n)
              {
                if (done[neighbors[n]] == (std::size_t)(-1))
                  {
                    if (m_assigned_type[neighbors[n]] == label)
                      {
                        todo.push (neighbors[n]);
                        done[neighbors[n]] = clusters.size()-1;
                      }
                    else
                      {
                        m_neighbor[current] = m_assigned_type[neighbors[n]];
                        m_neighbor[neighbors[n]] = m_assigned_type[current];
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
          pts.push_back (m_input[clusters[i].indices[j]]);
        clusters[i].centroid = CGAL::centroid (pts.begin(), pts.end());
        
      }

  }

  /// \cond SKIP_IN_MANUAL
  bool regularized_labeling_PC(const Neighborhood& neighborhood,
                               const FT radius_neighbors)
  {
    std::vector<std::vector<std::size_t> > groups;
    for (std::size_t i = 0; i < m_input.size(); ++ i)
      {
        std::size_t index = m_group[i];
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
    for (std::size_t s = 0; s < m_input.size(); s++)
      {
        map_p2i[m_input[s]] = s;

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
        m_assigned_type[s] = nb_class_best;
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
            m_assigned_type[groups[i][n]] = nb_class_best;
            for (std::size_t j = 0; j < mean.size(); ++ j)
              values[j][groups[i][n]] = mean[j] / groups[i].size();
          }
      }

    for (std::size_t s=0; s < m_input.size(); ++ s)
      {
        if (m_group[s] != (std::size_t)(-1))
          continue;
        
        std::vector<std::size_t> neighbors;
        neighborhood.range_neighbors (s, radius_neighbors,
                                      std::back_inserter (neighbors));

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

        m_assigned_type[s] = nb_class_best;

        std::sort (mean.begin(), mean.end());
        m_confidence[s] = mean[1] - mean[0];      

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
  void classify (int method, const Neighborhood& neighborhood = Neighborhood(),
                 const FT radius_neighbors = 0.1)
  {

    clock_t t;
    t = clock();

    build_effect_table ();
	
    CGAL_CLASSIFICATION_CERR<<std::endl<<"Classification of the point cloud: ";

    // Reset data structure
    std::vector<unsigned char>(m_input.size()).swap (m_assigned_type);
    std::vector<double>(m_input.size()).swap (m_confidence);
    
    if (method == 0)
      quick_labeling_PC();
    else if (method == 1)
      graphcut_labeling_PC(neighborhood);
    else if (method == 2)
      regularized_labeling_PC(neighborhood, radius_neighbors);
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


  /// \name Types and attributes
  /// @{
  
  /*!
    \brief Add a classification type
    \param type Pointer to the classification type object
   */
  void add_classification_type (Classification_type* type)
  {
    segmentation_classes.push_back (type);
  }

  void clear_classification_types ()
  {
    segmentation_classes.clear();
  }

  /*!
    \brief Add a segmentation attribute
    \param attribute Pointer to the attribute object
   */
  void add_segmentation_attribute (Segmentation_attribute* attribute)
  {
    segmentation_attributes.push_back (attribute);
  }

  void clear_segmentation_attributes ()
  {
    segmentation_attributes.clear();
  }
  
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
  void add_to_group (std::size_t point_index, std::size_t group_index)
  {
    m_group[point_index] = group_index;
  }

  /*!
    \brief Reset all groups attributes of points
   */
  void clear_groups()
  {
    m_group.clear();
    groups.clear();
  }

  /// @}
  

  /// \name Output
  /// @{

  /*!
    \brief Get the classification type of indexed point.
    \param index Index of the input point
    \return Pointer to the classification type 
  */
  Classification_type* classification_type_of (std::size_t index) const
  {
    return segmentation_classes[m_assigned_type[index]];
  }

  /// @}



};





} // namespace CGAL

#endif // CGAL_POINT_SET_CLASSIFICATION_H

