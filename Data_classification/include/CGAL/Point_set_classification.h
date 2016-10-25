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

#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>

#define CGAL_CLASSIFICATION_VERBOSE
#if defined(CGAL_CLASSIFICATION_VERBOSE)
#define CGAL_CLASSIFICATION_CERR std::cerr
#else
#define CGAL_CLASSIFICATION_CERR std::ostream(0)
#endif

//#define CGAL_CLASSTRAINING_VERBOSE
#if defined(CGAL_CLASSTRAINING_VERBOSE)
#define CGAL_CLASSTRAINING_CERR std::cerr
#else
#define CGAL_CLASSTRAINING_CERR std::ostream(0)
#endif

namespace CGAL {

/*!
\ingroup PkgDataClassification

\brief Classifies a point set based on a set of attribute and a set of classification types.

This class implement the core of the algorithm. It uses a point set as
input and assign each input point to a classification type among a set
of user defined classification types.

To achieve this classification algorithm, a set of local geometric
attributes are used, such as:

- planarity
- elevation
- vertical dispersion

The user must define a set of classification types such as:

- building
- ground
- vegetation

Each pair of attribute/type must be assigned an effect (for example,
vegetation has a low planarity and a high vertical dispersion) and
each attribute must be assigned a weight. These parameters can be set
up by hand or by providing a training set for each classification
type.

\tparam Kernel The geometric kernel used
\tparam RandomAccessIterator Iterator over the input
\tparam PointMap Property map to access the input points

*/
template <typename Kernel,
          typename RandomAccessIterator,
          typename PointMap>
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
                                            PointMap> Neighborhood;

  
#ifdef CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
  typedef internal::Alpha_expansion_graph_cut_boost             Alpha_expansion;
#else
  typedef internal::Alpha_expansion_graph_cut_boykov_kolmogorov Alpha_expansion;
#endif
  
private:
  
  class Point_range
  {
    std::size_t m_size;
    RandomAccessIterator m_begin;
    PointMap m_point_map;
    
  public:

    Point_range (RandomAccessIterator begin, RandomAccessIterator end,
                 PointMap point_map)
      : m_size (end - begin), m_begin (begin), m_point_map (point_map)
    { }

    std::size_t size() const { return m_size; }
    
    const Point& operator[] (std::size_t index) const { return get (m_point_map, *(m_begin + index)); }

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
  std::vector<std::size_t> m_training_type;
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


  /// \name Constructor
  /// @{
  
  /*! 
    \brief Constructs a classification object based on the input iterators.

    This method just initializes the structure and does not compute
    anything.

    \param begin Iterator to the first input object

    \param end Past-the-end iterator

    \param point_map Property map to access the input points

  */
  Point_set_classification (RandomAccessIterator begin,
                            RandomAccessIterator end,
                            PointMap point_map)
    : m_input (begin, end, point_map)
  {
    m_multiplicative = false;
  }

  /// @}


  /// \cond SKIP_IN_MANUAL

  double classification_value (std::size_t class_type, int pt_index) const
  {
    double out = 0.;
    if (m_multiplicative)
      {
        out = 1.;
        for (std::size_t i = 0; i < m_effect_table[class_type].size(); ++ i)
          {
            if (m_attributes[i]->weight == 0.)
              continue;
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
            if (m_attributes[i]->weight == 0.)
              continue;
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


  /// \name Classification
  /// @{

  
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

        double val_class_best = (std::numeric_limits<double>::max)();
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
    \brief Runs the classification algorithm with a local smoothing.

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
        double val_class_best = (std::numeric_limits<double>::max)();
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
    \brief Runs the classification algorithm with a global
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
    CGAL_CLASSIFICATION_CERR << "Labeling with regularization weight " << weight << "... ";

    std::vector<std::pair<std::size_t, std::size_t> > edges;
    std::vector<double> edge_weights;
    std::vector<std::vector<double> > probability_matrix
      (m_effect_table.size(), std::vector<double>(m_input.size(), 0.));
    std::vector<std::size_t>(m_input.size()).swap(m_assigned_type);

    for (std::size_t s = 0; s < m_input.size(); ++ s)
      {
        std::vector<std::size_t> neighbors;

        neighborhood.k_neighbors (m_input[s], 12, std::back_inserter (neighbors));

        for (std::size_t i = 0; i < neighbors.size(); ++ i)
          if (s != neighbors[i])
            {
              edges.push_back (std::make_pair (s, neighbors[i]));
              edge_weights.push_back (weight);
            }
        
        std::size_t nb_class_best = 0;
        double val_class_best = (std::numeric_limits<double>::max)();
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
        m_assigned_type[s] = nb_class_best;
      }
    
    Alpha_expansion graphcut;
    graphcut(edges, edge_weights, probability_matrix, m_assigned_type);
  }
  
  /// @}
  
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
            neighborhood.range_neighbors (m_input[current], tolerance,
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
        double val_class_best = (std::numeric_limits<double>::max)();
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

        double val_class_best = (std::numeric_limits<double>::max)();

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
        neighborhood.range_neighbors (m_input[s], radius_neighbors,
                                      std::back_inserter (neighbors));

        std::vector<double> mean (values.size(), 0.);
        for (std::size_t n = 0; n < neighbors.size(); ++ n)
          for (std::size_t j = 0; j < values.size(); ++ j)
            mean[j] += values[j][neighbors[n]];

        int nb_class_best=0; 
        double val_class_best = (std::numeric_limits<double>::max)();
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


  /// \name Classification Types
  /// @{
  
  /*!
    \brief Instanciates and adds a classification type.

    \param name ID of the classification type.

    \return A handle to the newly added classification type.
   */
  Data_classification::Type_handle add_classification_type (const char* name)
  {
    Data_classification::Type_handle out (new Data_classification::Type (name));
    m_types.push_back (out);
    return out;
  }

  
  /*!
    \brief Adds a classification type.

    \param type The handle to the classification type that must be added.
   */
  void add_classification_type (Data_classification::Type_handle type)
  {
    m_types.push_back (type);
  }

  /*!
    \brief Removes a classification type.

    \param type The handle to the classification type that must be removed.

    \return `true` if the classification type was correctly removed,
    `false` if its handle was not found inside the object.
   */ 
 bool remove_classification_type (Data_classification::Type_handle type)
  {
    for (std::size_t i = 0; i < m_types.size(); ++ i)
      if (m_types[i] == type)
        {
          m_types.erase (m_types.begin() + i);
          return true;
        }
    return false;
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
    \brief Removes all classification types.
   */
  void clear_classification_types ()
  {
    m_types.clear();
  }

  /// @}

  /// \name Attributes
  /// @{

  /*!
    \brief Adds an attribute.

    \param attribute Handle of the attribute to add.
   */
  void add_attribute (Data_classification::Attribute_handle attribute)
  {
    m_attributes.push_back (attribute);
  }

  /*!
    \brief Removes all attributes.
   */
  void clear_attributes ()
  {
    m_attributes.clear();
  }

  /// \cond SKIP_IN_MANUAL
  std::size_t number_of_attributes() const
  {
    return m_attributes.size();
  }

  Data_classification::Attribute_handle get_attribute(std::size_t idx)
  {
    return m_attributes[idx];
  }
  /// \endcond
  
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

    \note If classification was not performed (using `run()`,
    `run_with_local_smoothing()` or `run_with_graphcut()`), this
    function always returns the default empty `Type_handle`.

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

    \note If classification was not performed (using `run()`,
    `run_with_local_smoothing()` or `run_with_graphcut()`), this
    function always returns 0.

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
    function.

    Each classification type must be given a small set of user-defined
    inliers to provide the training algorithm with a ground truth.

    \param nb_tests Number of tests to perform. Higher values may
    provide the user with better results at the cost of higher
    computation time. Using a value of at least 10 times the number of
    attributes is advised.

    \return Minimum ratio (over all classification types) of provided
    ground truth points correctly classified using the best
    configuration found.
  */
  double training (std::size_t nb_tests = 300)
  {
    if (m_training_type.empty())
      return 0.;

    std::vector<std::vector<std::size_t> > training_sets (m_types.size());
    for (std::size_t i = 0; i < m_training_type.size(); ++ i)
      if (m_training_type[i] != (std::size_t)(-1))
        training_sets[m_training_type[i]].push_back (i);

    for (std::size_t i = 0; i < training_sets.size(); ++ i)
      if (training_sets[i].empty())
        std::cerr << "WARNING: \"" << m_types[i]->id() << "\" doesn't have a training set." << std::endl;

    std::vector<double> best_weights (m_attributes.size(), 1.);

    struct Attribute_training
    {
      bool skipped;
      double wmin;
      double wmax;
    };
    std::vector<Attribute_training> att_train;
    std::size_t nb_trials = 100;
    double wmin = 1e-5, wmax = 1e5;
    double factor = std::pow (wmax/wmin, 1. / (double)nb_trials);
    for (std::size_t j = 0; j < m_attributes.size(); ++ j)
      {
        Data_classification::Attribute_handle att = m_attributes[j];
        best_weights[j] = att->weight;

        std::size_t nb_useful = 0;
        double min = (std::numeric_limits<double>::max)();
        double max = -(std::numeric_limits<double>::max)();

        att->weight = wmin;
        for (std::size_t i = 0; i < 100; ++ i)
          {
            estimate_attribute_effect(training_sets, att);
            if (attribute_useful(att))
              {
                CGAL_CLASSTRAINING_CERR << "#";
                nb_useful ++;
                min = (std::min) (min, att->weight);
                max = (std::max) (max, att->weight);
              }
            else
              CGAL_CLASSTRAINING_CERR << "-";
            att->weight *= factor;
          }
        CGAL_CLASSTRAINING_CERR << std::endl;
        CGAL_CLASSTRAINING_CERR << att->id() << " useful in "
                  << nb_useful << "% of the cases, in interval [ "
                  << min << " ; " << max << " ]" << std::endl;
        att_train.push_back (Attribute_training());
        att_train.back().skipped = false;
        att_train.back().wmin = min / factor;
        att_train.back().wmax = max * factor;
        if (nb_useful < 2)
          {
            att_train.back().skipped = true;
            att->weight = 0.;
            best_weights[j] = att->weight;
          }
        else if (best_weights[j] == 1.)
          {
            att->weight = att_train.back().wmin
              + (rand() / (double)RAND_MAX) * (att_train.back().wmax - att_train.back().wmin);
            best_weights[j] = att->weight;
          }
        else
          {
            att->weight = best_weights[j];
          }
        estimate_attribute_effect(training_sets, att);
      }

    prepare_classification();
    
    double best_score = training_compute_worst_score(training_sets, 0.);
    double best_confidence = training_compute_worst_confidence(training_sets, 0.);
    
    std::cerr << "TRAINING GLOBALLY: Best score evolution: " << std::endl;

    std::cerr << 100. * best_score << "% (found at initialization)" << std::endl;

    bool first_round = true;    
    std::size_t current_att_changed = 0;
    std::ofstream f("score.plot");
    for (std::size_t i = 0; i < nb_tests; ++ i)
      {
        for (std::size_t j = 0; j < m_attributes.size(); ++ j)
          {
            Data_classification::Attribute_handle att = m_attributes[j];
            att->weight = best_weights[j];
          }
        Data_classification::Attribute_handle att = m_attributes[current_att_changed];            
        const Attribute_training& tr = att_train[current_att_changed];
        if (first_round)
          att->weight = 0.;
        else
          att->weight = tr.wmin + (rand() / (double)RAND_MAX) * (tr.wmax - tr.wmin);
        
        estimate_attributes_effects(training_sets);
        std::size_t nb_used = 0;
        for (std::size_t j = 0; j < m_attributes.size(); ++ j)
          {
            if (attribute_useful(m_attributes[j]))
              nb_used ++;
            else
              m_attributes[j]->weight = 0;
          }
        
        prepare_classification();
        double worst_confidence = training_compute_worst_confidence(training_sets,
                                                                    best_confidence);

        double worst_score = training_compute_worst_score(training_sets,
                                                          best_score);

        if (worst_score > best_score
            && worst_confidence > best_confidence)
          {
            best_score = worst_score;
            best_confidence = worst_confidence;
            std::cerr << 100. * best_score << "% (found at iteration "
                      << i+1 << "/" << nb_tests << ", "
                      << nb_used
                      << "/" << m_attributes.size() << " attribute(s) used)" << std::endl;
            for (std::size_t j = 0; j < m_attributes.size(); ++ j)
              {
                Data_classification::Attribute_handle att = m_attributes[j];
                best_weights[j] = att->weight;
                if (!(att_train[j].skipped))
                  f << att->weight << " ";
              }
            f << std::endl;
          }
        // if (worst_score > best_score)
        //   worst_score = best_score - (best_confidence - worst_confidence);
        // f << worst_score << " " << best_score << std::endl;

        do
          {
            ++ current_att_changed;
            if (current_att_changed == m_attributes.size())
              {
                first_round = false;
                current_att_changed = 0;
              }
          }
        while (att_train[current_att_changed].skipped);
      }

    for (std::size_t i = 0; i < best_weights.size(); ++ i)
      {
        Data_classification::Attribute_handle att = m_attributes[i];
        att->weight = best_weights[i];
      }

    estimate_attributes_effects(training_sets);
    
    std::cerr << std::endl << "Best score found is at least " << 100. * best_score
              << "% of correct classification" << std::endl;

    std::size_t nb_removed = 0;
    for (std::size_t i = 0; i < best_weights.size(); ++ i)
      {
        Data_classification::Attribute_handle att = m_attributes[i];
        CGAL_CLASSTRAINING_CERR << "ATTRIBUTE " << att->id() << ": " << best_weights[i] << std::endl;
        att->weight = best_weights[i];

        Data_classification::Type::Attribute_effect side = m_types[0]->attribute_effect(att);
        bool to_remove = true;
        for (std::size_t j = 0; j < m_types.size(); ++ j)
          {
            Data_classification::Type_handle ctype = m_types[j];
            if (ctype->attribute_effect(att) == Data_classification::Type::FAVORED_ATT)
              CGAL_CLASSTRAINING_CERR << " * Favored for ";
            else if (ctype->attribute_effect(att) == Data_classification::Type::PENALIZED_ATT)
              CGAL_CLASSTRAINING_CERR << " * Penalized for ";
            else
              CGAL_CLASSTRAINING_CERR << " * Neutral for ";
            if (ctype->attribute_effect(att) != side)
              to_remove = false;
            CGAL_CLASSTRAINING_CERR << ctype->id() << std::endl;
          }
        if (to_remove)
          {
            CGAL_CLASSTRAINING_CERR << "   -> Useless! Should be removed" << std::endl;
            ++ nb_removed;
          }
      }
    std::cerr << nb_removed
              << " attribute(s) out of " << m_attributes.size() << " are useless" << std::endl;
    
    return best_score;
  }

  /*!
    \brief Resets training sets.
  */
  void reset_training_sets()
  {
    std::vector<std::size_t>(m_input.size(), (std::size_t)(-1)).swap (m_training_type);
  }

  /*!
    \brief Adds the input point specified by index `idx` as an inlier
    of `class_type` for the training algorithm.

    \param class_type Handle to the classification type.

    \param idx Index of the input point.
  */
  bool add_training_index (Data_classification::Type_handle class_type,
                           std::size_t idx)
  {
    std::size_t type_idx = (std::size_t)(-1);
    for (std::size_t i = 0; i < m_types.size(); ++ i)
      if (m_types[i] == class_type)
        {
          type_idx = i;
          break;
        }
    if (type_idx == (std::size_t)(-1))
      return false;

    if (m_training_type.empty())
      reset_training_sets();

    m_training_type[idx] = type_idx;
    return true;
  }

  /*!
    \brief Adds input points specified by a range of indices as
    inliers of `class_type` for the training algorithm.

    \param class_type Handle to the classification type.

    \tparam IndexIterator Iterator with `std::size_t` as a
    `value_type`. \cgalModels InputIterator
  */
  template <class IndexIterator>
  bool add_training_set (Data_classification::Type_handle class_type,
                         IndexIterator first,
                         IndexIterator beyond)
  {
    std::size_t type_idx = (std::size_t)(-1);
    for (std::size_t i = 0; i < m_types.size(); ++ i)
      if (m_types[i] == class_type)
        {
          type_idx = i;
          break;
        }
    if (type_idx == (std::size_t)(-1))
      return false;

    if (m_training_type.empty())
      reset_training_sets();
    
    for (IndexIterator it = first; it != beyond; ++ it)
      m_training_type[*it] = type_idx;

    return true;
  }
  
  /// @}

  
  /// \cond SKIP_IN_MANUAL
  Data_classification::Type_handle training_type_of (std::size_t index) const
  {
    if (m_training_type.size() <= index
        || m_training_type[index] == (std::size_t)(-1))
      return Data_classification::Type_handle();
    return m_types[m_training_type[index]];
  }

  void estimate_attributes_effects
  (const std::vector<std::vector<std::size_t> >& training_sets)
  {
    for (std::size_t i = 0; i < m_attributes.size(); ++ i)
      estimate_attribute_effect (training_sets, m_attributes[i]);
  }

  void estimate_attribute_effect
  (const std::vector<std::vector<std::size_t> >& training_sets,
   Data_classification::Attribute_handle att)
  {
    std::vector<double> mean (m_types.size(), 0.);
                                  
    for (std::size_t j = 0; j < m_types.size(); ++ j)
      {
        for (std::size_t k = 0; k < training_sets[j].size(); ++ k)
          {
            double val = att->normalized(training_sets[j][k]);
            mean[j] += val;
          }
        mean[j] /= training_sets[j].size();
      }

    std::vector<double> sd (m_types.size(), 0.);
        
    for (std::size_t j = 0; j < m_types.size(); ++ j)
      {
        Data_classification::Type_handle ctype = m_types[j];
            
        for (std::size_t k = 0; k < training_sets[j].size(); ++ k)
          {
            double val = att->normalized(training_sets[j][k]);
            sd[j] += (val - mean[j]) * (val - mean[j]);
          }
        sd[j] = std::sqrt (sd[j] / training_sets[j].size());
        if (mean[j] - sd[j] > 0.5)
          ctype->set_attribute_effect (att, Data_classification::Type::FAVORED_ATT);
        else if (mean[j] + sd[j] < 0.5)
          ctype->set_attribute_effect (att, Data_classification::Type::PENALIZED_ATT);
        else
          ctype->set_attribute_effect (att, Data_classification::Type::NEUTRAL_ATT);
      }
  }

  
  double training_compute_worst_score
  (const std::vector<std::vector<std::size_t> >& training_sets,
   double lower_bound)
  {
    double worst_score = 1.;
    for (std::size_t j = 0; j < m_types.size(); ++ j)
      {
        std::size_t nb_okay = 0;
        for (std::size_t k = 0; k < training_sets[j].size(); ++ k)
          {
            std::size_t nb_class_best=0; 
            double val_class_best = (std::numeric_limits<double>::max)();
      
            for(std::size_t l = 0; l < m_effect_table.size(); ++ l)
              {
                double value = classification_value (l, training_sets[j][k]);
          
                if(val_class_best > value)
                  {
                    val_class_best = value;
                    nb_class_best = l;
                  }
              }
                
            if (nb_class_best == j)
              nb_okay ++;

          }

        double score = nb_okay / (double)(training_sets[j].size());
        if (score < worst_score)
          worst_score = score;
        if (worst_score < lower_bound)
          return worst_score;
      }
    return worst_score;
  }

  double training_compute_worst_confidence
  (const std::vector<std::vector<std::size_t> >& training_sets,
   double lower_bound)
  {
    double worst_confidence = (std::numeric_limits<double>::max)();
    for (std::size_t j = 0; j < m_types.size(); ++ j)
      {
        double confidence = 0.;
        
        for (std::size_t k = 0; k < training_sets[j].size(); ++ k)
          {
            std::vector<std::pair<double, std::size_t> > values;
      
            for(std::size_t l = 0; l < m_effect_table.size(); ++ l)
              values.push_back (std::make_pair (classification_value (l, training_sets[j][k]),
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

        confidence /= (double)(training_sets[j].size() * m_attributes.size());

        if (confidence < worst_confidence)
          worst_confidence = confidence;
        if (worst_confidence < lower_bound)
          return worst_confidence;
      }
    return worst_confidence;
  }

  bool attribute_useful (Data_classification::Attribute_handle att)
  {
    Data_classification::Type::Attribute_effect side = m_types[0]->attribute_effect(att);
    for (std::size_t k = 1; k < m_types.size(); ++ k)
      if (m_types[k]->attribute_effect(att) != side)
        return true;
    return false;
  }

  /// \endcond


};





} // namespace CGAL

#endif // CGAL_POINT_SET_CLASSIFICATION_H

