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
#include <CGAL/Data_classification/Attribute.h>
#include <CGAL/Data_classification/Type.h>

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

\brief Classifies a point set based on a set of attribute and a set of classification types.

This class implement the core of the algorithm. It uses a point set as
input. Based on a set of segmentation attributes and a set of
classification types, it segments the point set into the different
types given. The output can be regularized with different smoothing
methods.

\tparam Kernel The geometric kernel used
\tparam RandomAccessIterator Iterator over the input
\tparam PointPMap Property map to access the input points

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

  
  std::vector<Data_classification::Type_handle> m_types; 
  std::vector<Data_classification::Attribute_handle> m_attributes; 

  typedef Data_classification::Type::Attribute_effect Attribute_effect;
  std::vector<std::vector<Attribute_effect> > m_effect_table;

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

    \param begin Iterator to the first input object

    \param end Past-the-end iterator

    \param point_pmap Property map to access the input points

  */

  Point_set_classification (RandomAccessIterator begin,
                            RandomAccessIterator end,
                            PointPMap point_pmap)
    : m_input (begin, end, point_pmap)
  {
    m_multiplicative = false;
  }

  /// \cond SKIP_IN_MANUAL

  double classification_value (std::size_t class_type, int pt_index) const
  {
    double out = 0.;
    if (m_multiplicative)
      {
        out = 1.;
        for (std::size_t i = 0; i < m_effect_table[class_type].size(); ++ i)
          {
            if (m_effect_table[class_type][i] == Data_classification::Type::FAVORED_ATT)
              out *= 1. + m_attributes[i]->favored (pt_index);
            else if (m_effect_table[class_type][i] == Data_classification::Type::PENALIZED_ATT)
              out *= 1. + m_attributes[i]->penalized (pt_index);
            else if (m_effect_table[class_type][i] == Data_classification::Type::NEUTRAL_ATT)
              out *= 1. + m_attributes[i]->ignored (pt_index);
          }
      }
    else
      {
        for (std::size_t i = 0; i < m_effect_table[class_type].size(); ++ i)
          {
            if (m_effect_table[class_type][i] == Data_classification::Type::FAVORED_ATT)
              out += m_attributes[i]->favored (pt_index);
            else if (m_effect_table[class_type][i] == Data_classification::Type::PENALIZED_ATT)
              out += m_attributes[i]->penalized (pt_index);
            else if (m_effect_table[class_type][i] == Data_classification::Type::NEUTRAL_ATT)
              out += m_attributes[i]->ignored (pt_index);
          }
      }
    return out;
  }

  void set_multiplicative (bool mult)
  {
    m_multiplicative = mult;
  }
  /// \endcond

  
  /*! 
    \brief Runs the classification algorithm without any regularization.

    There is no relationship between points, the classification energy
    is only minimized pointwise. This method is quick but produce
    suboptimal results.
  */
  void run()
  {
    prepare_classification ();
    
    // data term initialisation

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
      }
  }


  /*! 
    \brief Runs the classification algorithm without local smoothing.

    The computed classification energy is smoothed on a user defined
    local neighborhood of points. This method is a compromise between
    efficiency and reliability.

    \param neighborhood Object used to access neighborhoods of points
  */
  void run_with_local_smoothing (const Neighborhood& neighborhood)
  {
    prepare_classification ();
    
    // data term initialisation
    CGAL_CLASSIFICATION_CERR << "Labeling... ";

    std::vector<std::vector<double> > values
      (m_types.size(),
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


  /*! 

    \brief Runs the classification algorithm without a global
    regularizarion based on a graphcut.

    The computed classification energy is globally regularized through
    and alpha-expansion algorithm. This method is slow but provides
    the user with good quality results.

    \param neighborhood Object used to access neighborhoods of points

    \param weight Weight of the regularization with respect to the
    classification energy. Higher values produce more regularized
    output but may result in a loss of details.
  */
  void run_with_graphcut (const Neighborhood& neighborhood,
                          const double weight = 0.5)
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
              edge_weights.push_back (weight);
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
  

  /// \cond SKIP_IN_MANUAL
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
    values.resize (m_types.size());
    
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
    std::vector<std::size_t>(m_input.size(), (std::size_t)(-1)).swap (m_assigned_type);
    std::vector<double>(m_input.size()).swap (m_confidence);

    m_effect_table = std::vector<std::vector<Attribute_effect> >
      (m_types.size(), std::vector<Attribute_effect> (m_attributes.size(),
                                                                 Data_classification::Type::NEUTRAL_ATT));
    
    for (std::size_t i = 0; i < m_effect_table.size (); ++ i)
      for (std::size_t j = 0; j < m_effect_table[i].size (); ++ j)
        m_effect_table[i][j] = m_types[i]->attribute_effect (m_attributes[j]);

  }
  /// \endcond


  /// \name Types and attributes
  /// @{
  
  /*!
    \brief Adds a classification type

   */
  Data_classification::Type_handle add_classification_type (const char* name)
  {
    Data_classification::Type_handle out (new Data_classification::Type (name));
    m_types.push_back (out);
    return out;
  }

  
  /*!
    \brief Adds a classification type

   */
  void add_classification_type (Data_classification::Type_handle type)
  {
    m_types.push_back (type);
  }

  /// \cond SKIP_IN_MANUAL
  std::size_t number_of_classification_types () const
  {
    return m_types.size();
  }

  Data_classification::Type_handle get_classification_type (std::size_t idx)
  {
    return m_types[idx];
  }
  /// \endcond

  /*!
    \brief Removes all classification types
   */
  void clear_classification_types ()
  {
    m_types.clear();
  }

  /*!
    \brief Adds an attribute
    \param attribute Pointer to the attribute object
   */
  void add_attribute (Data_classification::Attribute_handle attribute)
  {
    m_attributes.push_back (attribute);
  }

  /*!
    \brief Removes all attributes
   */
  void clear_attributes ()
  {
    m_attributes.clear();
  }
  
  /// @}

  /// \cond SKIP_IN_MANUAL

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

  /// \endcond
  

  /// \name Output
  /// @{

  /*!
    \brief Gets the classification type of an indexed point.
    \param index Index of the input point
    \return Pointer to the classification type 
  */
  Data_classification::Type_handle classification_type_of (std::size_t index) const
  {
    if (m_assigned_type.size() <= index
        || m_assigned_type[index] == (std::size_t)(-1))
      return Data_classification::Type_handle();
    return m_types[m_assigned_type[index]];
  }

  /// \cond SKIP_IN_MANUAL
  bool classification_prepared() const
  {
    return !(m_assigned_type.empty());
  }
  void set_classification_type_of (std::size_t index, Data_classification::Type_handle class_type)
  {
    for (std::size_t i = 0; i < m_types.size(); ++ i)
      if (m_types[i] == class_type)
        {
          m_assigned_type[index] = i;
          return;
        }
    m_assigned_type[index] = (std::size_t)(-1);
  }
  /// \endcond

  /*!
    \brief Gets the confidence of the classification type of an indexed point.
    \param index Index of the input point
    \return Confidence ranging from 0 (not confident at all) to 1 (very confident).
  */
  double confidence_of (std::size_t index) const
  {
    if (m_confidence.size() <= index)
      return 0.;
    return m_confidence[index];
  }

  /// @}


  /// \name Training
  /// @{

  /*!
    \brief Runs the training algorithm.

    The object must have been filled with the `Data_classification::Type`
    and `Data_classification::Attribute` objects before running this
    function. Each classification type must be given a small set of
    user-defined inliers to provide the training algorithm with a
    ground truth.
  */
  void training (std::size_t nb_tests = 3000)
  {
    train_parameters(nb_tests);
  }
  /// @}

  
  /// \cond SKIP_IN_MANUAL
  void train_parameters(std::size_t nb_tests)
  {
    std::vector<double> best_weights (m_attributes.size(), 1.);
    
    double best_score = -1.;
    double best_confidence = -1.;
    
    std::cerr << "TRAINING GLOBALLY: Best score evolution: " << std::endl;

    std::size_t current_att_changed = 0;
    std::ofstream f("score.plot");
    for (std::size_t i = 0; i < nb_tests; ++ i)
      {
        if (!(i % 10))
          {
            for (std::size_t j = 0; j < m_attributes.size(); ++ j)
              {
                Data_classification::Attribute_handle att = m_attributes[j];
                att->weight = (rand() / (double)RAND_MAX) * 3. * att->max;
              }
          }
        else
          {
            for (std::size_t j = 0; j < m_attributes.size(); ++ j)
              {
                Data_classification::Attribute_handle att = m_attributes[j];
                att->weight = best_weights[j];
              }
            Data_classification::Attribute_handle att = m_attributes[current_att_changed];
            att->weight = (rand() / (double)RAND_MAX) * 3. * att->max;
          }
        
        estimate_attribute_effects();
        std::vector<Data_classification::Attribute_handle> used_attributes;
        for (std::size_t j = 0; j < m_attributes.size(); ++ j)
          {
            Data_classification::Attribute_handle att = m_attributes[j];
            Data_classification::Type::Attribute_effect side = m_types[0]->attribute_effect(att);

            for (std::size_t k = 1; k < m_types.size(); ++ k)
              if (m_types[k]->attribute_effect(att) != side)
                {
                  used_attributes.push_back (att);
                  break;
                }
          }

        used_attributes.swap (m_attributes);
        
        prepare_classification();
        double worst_confidence = training_compute_worst_confidence();

        double worst_score = training_compute_worst_score();
        used_attributes.swap (m_attributes);
        
        if (worst_score > best_score
            && worst_confidence > best_confidence)
          {
            best_score = worst_score;
            best_confidence = worst_confidence;
            std::cerr << 100. * best_score << "% (found at iteration "
                      << i+1 << "/" << nb_tests << ") "
                      << used_attributes.size() << std::endl;
            for (std::size_t j = 0; j < m_attributes.size(); ++ j)
              {
                Data_classification::Attribute_handle att = m_attributes[j];
                best_weights[j] = att->weight;
                f << best_weights[j] << " ";
              }
            f << std::endl;
          }

        if (!(i % 10))
          {
            ++ current_att_changed;
            if (current_att_changed == m_attributes.size())
              current_att_changed = 0;
          }
      }

    for (std::size_t i = 0; i < best_weights.size(); ++ i)
      {
        Data_classification::Attribute_handle att = m_attributes[i];
        att->weight = best_weights[i];
      }

    estimate_attribute_effects();
    
    std::cerr << std::endl << "Best score found is at least " << 100. * best_score
              << "% of correct classification" << std::endl;
    std::vector<Data_classification::Attribute_handle> to_keep;
    
    for (std::size_t i = 0; i < best_weights.size(); ++ i)
      {
        Data_classification::Attribute_handle att = m_attributes[i];
        std::cerr << "ATTRIBUTE " << att->id() << ": " << best_weights[i] << std::endl;
        att->weight = best_weights[i];

        Data_classification::Type::Attribute_effect side = m_types[0]->attribute_effect(att);
        bool to_remove = true;
        for (std::size_t j = 0; j < m_types.size(); ++ j)
          {
            Data_classification::Type_handle ctype = m_types[j];
            if (ctype->attribute_effect(att) == Data_classification::Type::FAVORED_ATT)
              CGAL_CLASSIFICATION_CERR << " * Favored for ";
            else if (ctype->attribute_effect(att) == Data_classification::Type::PENALIZED_ATT)
              CGAL_CLASSIFICATION_CERR << " * Penalized for ";
            else
              CGAL_CLASSIFICATION_CERR << " * Neutral for ";
            if (ctype->attribute_effect(att) != side)
              to_remove = false;
            CGAL_CLASSIFICATION_CERR << ctype->id() << std::endl;
          }
        if (to_remove)
          {
            std::cerr << "   -> Useless! Should be removed" << std::endl;
            //            delete att;
          }
        else
          to_keep.push_back (att);
      }
    std::cerr << "Removing " << (m_attributes.size() - to_keep.size())
              << " attribute(s) out of " << m_attributes.size() << std::endl;
    m_attributes.swap (to_keep);
  }

  void estimate_attribute_effects()
  {
    for (std::size_t i = 0; i < m_attributes.size(); ++ i)
      {
        Data_classification::Attribute_handle att = m_attributes[i];

        std::vector<double> mean (m_types.size(), 0.);
                                  
        for (std::size_t j = 0; j < m_types.size(); ++ j)
          {
            Data_classification::Type_handle ctype = m_types[j];
            
            for (std::size_t k = 0; k < ctype->training_set().size(); ++ k)
              {
                double val = att->normalized(ctype->training_set()[k]);
                mean[j] += val;
              }
            mean[j] /= ctype->training_set().size();
          }

        std::vector<double> sd (m_types.size(), 0.);
        
        for (std::size_t j = 0; j < m_types.size(); ++ j)
          {
            Data_classification::Type_handle ctype = m_types[j];
            
            for (std::size_t k = 0; k < ctype->training_set().size(); ++ k)
              {
                double val = att->normalized(ctype->training_set()[k]);
                sd[j] += (val - mean[j]) * (val - mean[j]);
              }
            sd[j] = std::sqrt (sd[j] / ctype->training_set().size());
            if (mean[j] - sd[j] > 0.5)
              ctype->set_attribute_effect (att, Data_classification::Type::FAVORED_ATT);
            else if (mean[j] + sd[j] < 0.5)
              ctype->set_attribute_effect (att, Data_classification::Type::PENALIZED_ATT);
            else
              ctype->set_attribute_effect (att, Data_classification::Type::NEUTRAL_ATT);
          }
      }
  }

  
  double training_compute_worst_score()
  {
    double worst_score = 1.;
    for (std::size_t j = 0; j < m_types.size(); ++ j)
      {
        Data_classification::Type_handle ctype = m_types[j];
        std::size_t nb_okay = 0;
        for (std::size_t k = 0; k < ctype->training_set().size(); ++ k)
          {
            std::size_t nb_class_best=0; 
            double val_class_best = std::numeric_limits<double>::max();
      
            for(std::size_t l = 0; l < m_effect_table.size(); ++ l)
              {
                double value = classification_value (l, ctype->training_set()[k]);
          
                if(val_class_best > value)
                  {
                    val_class_best = value;
                    nb_class_best = l;
                  }
              }
                
            if (nb_class_best == j)
              nb_okay ++;

          }

        double score = nb_okay / (double)(ctype->training_set().size());
        if (score < worst_score)
          worst_score = score;
      }
    return worst_score;
  }

  double training_compute_worst_confidence()
  {
    double worst_confidence = std::numeric_limits<double>::max();
    for (std::size_t j = 0; j < m_types.size(); ++ j)
      {
        Data_classification::Type_handle ctype = m_types[j];
        
        double confidence = 0.;
        
        for (std::size_t k = 0; k < ctype->training_set().size(); ++ k)
          {
            std::vector<std::pair<double, std::size_t> > values;
      
            for(std::size_t l = 0; l < m_effect_table.size(); ++ l)
              values.push_back (std::make_pair (classification_value (l, ctype->training_set()[k]),
                                                l));
            std::sort (values.begin(), values.end());

            if (values[0].second == j)
              confidence += values[1].first - values[0].first;
            else
              {
                // for(std::size_t l = 0; l < values.size(); ++ l)
                //   if (values[l].second == j)
                //     {
                //       confidence += values[0].first - values[l].first;
                //       break;
                //     }
              }
            
          }

        confidence /= (double)(ctype->training_set().size() * m_attributes.size());

        if (confidence < worst_confidence)
          worst_confidence = confidence;
      }
    return worst_confidence;
  }
  /// \endcond


};





} // namespace CGAL

#endif // CGAL_POINT_SET_CLASSIFICATION_H

