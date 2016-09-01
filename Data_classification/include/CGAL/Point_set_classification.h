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
#include <CGAL/Data_classification/Segmentation_attribute.h>

#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>

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

  typedef Data_classification::Neighborhood<Kernel,
                                            RandomAccessIterator,
                                            PointPMap> Neighborhood;

private:
  
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

  std::vector<std::size_t> m_group;
  std::vector<std::size_t> m_assigned_type;
  std::vector<std::size_t> m_neighbor;
  std::vector<double> m_confidence;

  struct Cluster
  {
    Point centroid;
    std::vector<std::size_t> indices;
    std::set<std::size_t> neighbors;
  };

  
  std::vector<Classification_type*> m_segmentation_classes; 
  std::vector<Segmentation_attribute*> m_segmentation_attributes; 

  typedef Classification_type::Attribute_side Attribute_side;
  std::vector<std::vector<Attribute_side> > m_effect_table;

  //Hpoints attributes
  std::vector<Plane> m_planes;
  std::vector<Cluster> m_clusters;
  
  bool m_multiplicative;
  /// \endcond

public:


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
                            PointPMap point_pmap)
    : m_input (begin, end, point_pmap)
  {
    m_multiplicative = false;
  }

  /// \cond SKIP_IN_MANUAL

  double classification_value (std::size_t segmentation_class, int pt_index) const
  {
    double out = 0.;
    if (m_multiplicative)
      {
        out = 1.;
        for (std::size_t i = 0; i < m_effect_table[segmentation_class].size(); ++ i)
          {
            if (m_effect_table[segmentation_class][i] == Classification_type::FAVORED_ATT)
              out *= m_segmentation_attributes[i]->favored (pt_index);
            else if (m_effect_table[segmentation_class][i] == Classification_type::PENALIZED_ATT)
              out *= m_segmentation_attributes[i]->penalized (pt_index);
            else if (m_effect_table[segmentation_class][i] == Classification_type::NEUTRAL_ATT)
              out *= m_segmentation_attributes[i]->ignored (pt_index);
          }
      }
    else
      {
        for (std::size_t i = 0; i < m_effect_table[segmentation_class].size(); ++ i)
          {
            if (m_effect_table[segmentation_class][i] == Classification_type::FAVORED_ATT)
              out += m_segmentation_attributes[i]->favored (pt_index);
            else if (m_effect_table[segmentation_class][i] == Classification_type::PENALIZED_ATT)
              out += m_segmentation_attributes[i]->penalized (pt_index);
            else if (m_effect_table[segmentation_class][i] == Classification_type::NEUTRAL_ATT)
              out += m_segmentation_attributes[i]->ignored (pt_index);
          }
      }
    return out;
  }

  void set_multiplicative (bool mult)
  {
    m_multiplicative = mult;
  }

  void run_quick()
  {
    prepare_classification ();
    
    // data term initialisation
    CGAL_CLASSIFICATION_CERR << "Labeling... ";

    int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
    for (std::size_t s = 0; s < m_input.size(); s++)
      {
			
        int nb_class_best=0; 

        double val_class_best = std::numeric_limits<double>::max();
        std::vector<double> values;
      
        for(std::size_t k = 0; k < m_effect_table.size(); ++ k)
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
  }


  void run_with_local_smoothing (const Neighborhood& neighborhood)
  {
    prepare_classification ();
    
    // data term initialisation
    CGAL_CLASSIFICATION_CERR << "Labeling... ";

    std::vector<std::vector<double> > values
      (m_segmentation_classes.size(),
       std::vector<double> (m_input.size(), -1.));

    for (std::size_t s=0; s < m_input.size(); ++ s)
      {
        std::vector<std::size_t> neighbors;
        neighborhood.get (s, neighbors);

        std::vector<double> mean (values.size(), 0.);
        for (std::size_t n = 0; n < neighbors.size(); ++ n)
          {
            if (values[0][neighbors[n]] < 0.)
              for(std::size_t k = 0; k < m_effect_table.size(); ++ k)
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

  }


  void run_with_graphcut (const Neighborhood& neighborhood,
                          const double smoothing)
  {
    prepare_classification ();
    
    // data term initialisation
    CGAL_CLASSIFICATION_CERR << "Labeling... ";

    std::vector<std::pair<std::size_t, std::size_t> > edges;
    std::vector<double> edge_weights;
    std::vector<std::vector<double> > probability_matrix
      (m_effect_table.size(), std::vector<double>(m_input.size(), 0.));
    std::vector<std::size_t>(m_input.size()).swap(m_assigned_type);

    for (std::size_t s = 0; s < m_input.size(); ++ s)
      {
        std::vector<std::size_t> neighbors;

        neighborhood.k_neighbors (s, 12, std::back_inserter (neighbors));

        for (std::size_t i = 0; i < neighbors.size(); ++ i)
          {
            edges.push_back (std::make_pair (s, neighbors[i]));
            if (s != neighbors[i])
              edge_weights.push_back (smoothing);
            else
              edge_weights.push_back (0.);
          }
        
        std::size_t nb_class_best = 0;
        double val_class_best = std::numeric_limits<double>::max();
        for(std::size_t k = 0; k < m_effect_table.size(); ++ k)
          {
            double value = classification_value (k, s);
            probability_matrix[k][s] = value;
            
            if(val_class_best > value)
              {
                val_class_best = value;
                nb_class_best = k;
              }
          }
        m_assigned_type.push_back (nb_class_best);
      }

    
    internal::Alpha_expansion_graph_cut_boost graphcut;
    graphcut(edges, edge_weights, probability_matrix, m_assigned_type);
  }
  
  void reset_groups()
  {
    m_planes.clear();
    std::vector<std::size_t>(m_input.size(), (std::size_t)(-1)).swap (m_group);
  }

  void add_plane (const Plane& p) { m_planes.push_back (p); }
  
  std::size_t group_of (std::size_t idx) const
  {
    if (m_group.size() <= idx)
      return (std::size_t)(-1);
    return m_group[idx];
  }
  void set_group_of (std::size_t idx, std::size_t idx_group) { m_group[idx] = idx_group; }

  const std::vector<Cluster>& clusters() const { return m_clusters; }

  std::size_t neighbor_of (std::size_t idx) const { return m_neighbor[idx]; }
  
  void cluster_points (const Neighborhood& neighborhood, const double& tolerance)
  {
    std::vector<std::size_t>(m_input.size(), (std::size_t)(-1)).swap (m_neighbor);
    
    std::vector<std::size_t> done (m_input.size(), (std::size_t)(-1));
    
    for (std::size_t s=0; s < m_input.size(); ++ s)
      {
        if (done[s] != (std::size_t)(-1))
          continue;
        std::size_t label = m_assigned_type[s];
        
        m_clusters.push_back (Cluster());

        std::queue<std::size_t> todo;
        todo.push (s);
        done[s] = m_clusters.size()-1;

        while (!(todo.empty()))
          {
            std::size_t current = todo.front();
            todo.pop();
            m_clusters.back().indices.push_back (current);
            
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
                        done[neighbors[n]] = m_clusters.size()-1;
                      }
                    else
                      {
                        m_neighbor[current] = m_assigned_type[neighbors[n]];
                        m_neighbor[neighbors[n]] = m_assigned_type[current];
                        continue;
                      }
                  }
                else if (done[neighbors[n]] != m_clusters.size()-1)
                  {
                    m_clusters.back().neighbors.insert (done[neighbors[n]]);
                    m_clusters[done[neighbors[n]]].neighbors.insert (m_clusters.size()-1);
                  }
              }
          }
      }
    std::cerr << "Found " << m_clusters.size() << " cluster(s)" << std::endl;

    for (std::size_t i = 0; i < m_clusters.size(); ++ i)
      {
        std::vector<Point> pts;
        for (std::size_t j = 0; j < m_clusters[i].indices.size(); ++ j)
          pts.push_back (m_input[m_clusters[i].indices[j]]);
        m_clusters[i].centroid = CGAL::centroid (pts.begin(), pts.end());
        
      }

  }

  /// \cond SKIP_IN_MANUAL
  bool run_with_groups (const Neighborhood& neighborhood,
                        const FT radius_neighbors)
  {
    prepare_classification ();
    
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
    values.resize (m_segmentation_classes.size());
    
    std::map<Point, std::size_t> map_p2i;
    for (std::size_t s = 0; s < m_input.size(); s++)
      {
        map_p2i[m_input[s]] = s;

        int nb_class_best=0; 
        double val_class_best = std::numeric_limits<double>::max();
        for(std::size_t k = 0; k < m_effect_table.size(); ++ k)
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


  /// @}

  
  /// \cond SKIP_IN_MANUAL
  void prepare_classification ()
  {
    // Reset data structure
    std::vector<std::size_t>(m_input.size()).swap (m_assigned_type);
    std::vector<double>(m_input.size()).swap (m_confidence);
    

    m_effect_table = std::vector<std::vector<Attribute_side> >
      (m_segmentation_classes.size(), std::vector<Attribute_side> (m_segmentation_attributes.size(),
                                                                 Classification_type::NEUTRAL_ATT));
    
    for (std::size_t i = 0; i < m_effect_table.size (); ++ i)
      for (std::size_t j = 0; j < m_effect_table[i].size (); ++ j)
        m_effect_table[i][j] = m_segmentation_classes[i]->attribute_effect (m_segmentation_attributes[j]);

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
    m_segmentation_classes.push_back (type);
  }

  void number_of_classification_types () const
  {
    return m_segmentation_classes.size();
  }

  Classification_type* get_classification_type (std::size_t idx)
  {
    return m_segmentation_classes[idx];
  }

  void clear_classification_types ()
  {
    m_segmentation_classes.clear();
  }

  void clear_and_delete_classification_types ()
  {
    for (std::size_t i = 0; i < m_segmentation_classes.size(); ++ i)
      delete m_segmentation_classes[i];
    m_segmentation_classes.clear();
  }

  /*!
    \brief Add a segmentation attribute
    \param attribute Pointer to the attribute object
   */
  void add_segmentation_attribute (Segmentation_attribute* attribute)
  {
    m_segmentation_attributes.push_back (attribute);
  }

  void clear_segmentation_attributes ()
  {
    m_segmentation_attributes.clear();
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
    m_planes.clear();
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
    if (m_assigned_type.size() <= index)
      return NULL;
    return m_segmentation_classes[m_assigned_type[index]];
  }

  double confidence_of (std::size_t index) const
  {
    if (m_confidence.size() <= index)
      return 0.;
    return m_confidence[index];
  }

  /// @}



};





} // namespace CGAL

#endif // CGAL_POINT_SET_CLASSIFICATION_H

