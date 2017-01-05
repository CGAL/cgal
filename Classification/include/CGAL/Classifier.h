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

#ifndef CGAL_CLASSIFIER_H
#define CGAL_CLASSIFIER_H

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

#include <CGAL/Classification/Planimetric_grid.h>
#include <CGAL/Classification/Attribute_base.h>
#include <CGAL/Classification/Type.h>

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
\ingroup PkgClassification

\brief Classifies a data set based on a set of attribute and a set of classification types.

This class implements the core of the classification algorithm
\cgalCite{cgal:lm-clscm-12}. It uses a data set as input and assignes
each input item to a classification type among a set of user
defined classification types.

To achieve this classification algorithm, a set of local geometric
attributes are used, such as planarity, elevation or vertical
dispersion.

The user must define a set of classification types such as building,
ground or vegetation.

Each pair of attribute and type must be assigned an
[Attribute::Effect](@ref CGAL::Classification::Attribute::Effect) (for
example, vegetation has a low planarity and a high vertical
dispersion) and each attribute must be assigned a weight. These
parameters can be set up by hand or by providing a training set for
each classification type.

\tparam Range range of items, model of `ConstRange`. Its iterator type
is `RandomAccessIterator`.

\tparam ItemMap model of `ReadablePropertyMap` whose key
type is the value type of the iterator of `Range` and value type is
`Point_3<Kernel>`.


*/
template <typename Range,
          typename ItemMap>
class Classifier
{

  
public:
  /// \cond SKIP_IN_MANUAL
  typedef typename ItemMap::value_type Item;

  typedef typename Classification::Type_handle      Type_handle;
  typedef typename Classification::Attribute_handle Attribute_handle;
  
#ifdef CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
  typedef internal::Alpha_expansion_graph_cut_boost             Alpha_expansion;
#else
  typedef internal::Alpha_expansion_graph_cut_boykov_kolmogorov Alpha_expansion;
#endif
  
private:
  
  const Range& m_input;
  ItemMap m_item_map;

  std::vector<std::size_t> m_assigned_type;
  std::vector<std::size_t> m_training_type;
  std::vector<double> m_confidence;

  std::vector<Type_handle> m_types; 
  std::vector<Attribute_handle> m_attributes; 

  typedef Classification::Attribute::Effect Attribute_effect;
  std::vector<std::vector<Attribute_effect> > m_effect_table;

  /// \endcond

public:


  /// \name Constructor
  /// @{
  
  /*! 
    \brief Initializes a classification object.

    \param input input range

    \param item_map property map to access the input items

  */
  Classifier (const Range& input,
              ItemMap item_map)
    : m_input (input), m_item_map (item_map)
  {
  }
  /// @}

  /// \name Classification
  /// @{

  
  /*! 
    \brief Runs the classification algorithm without any regularization.

    There is no relationship between items, the classification energy
    is only minimized itemwise. This method is quick but produce
    suboptimal results.
  */
  void run()
  {
    prepare_classification ();
    
    // data term initialisation

    for (std::size_t s = 0; s < m_input.size(); s++)
      {

        std::size_t nb_class_best=0; 

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
        m_confidence[s] = values[1] - values[0];
      }
  }


  /*! 
    \brief Runs the classification algorithm with a local smoothing.

    The computed classification energy is smoothed on a user defined
    local neighborhood of items. This method is a compromise between
    efficiency and reliability.

    \tparam NeighborQuery model of `NeighborQuery`
    \param neighbor_query used to access neighborhoods of items
  */
  template <typename NeighborQuery>
  void run_with_local_smoothing (const NeighborQuery& neighbor_query)
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
        neighbor_query (get (m_item_map, m_input[s]), std::back_inserter (neighbors));

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

        std::size_t nb_class_best=0; 
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
    regularization based on a graphcut.

    The computed classification energy is globally regularized through
    an alpha-expansion algorithm. This method is slow but provides
    the user with good quality results.

    \tparam NeighborQuery model of `NeighborQuery`
    \param neighbor_query used to access neighborhoods of items
    \param weight weight of the regularization with respect to the
    classification energy. Higher values produce more regularized
    output but may result in a loss of details.

  */
  template <typename NeighborQuery>
  void run_with_graphcut (const NeighborQuery& neighbor_query,
                          const double& weight)
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

        neighbor_query (get(m_item_map, m_input[s]), std::back_inserter (neighbors));

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
  


  /// \name Classification Types
  /// @{
  
  /*!
    \brief Instantiates and adds a classification type.

    \param name name of the classification type.

    \return a handle to the newly added classification type.
   */
  Type_handle add_classification_type (const char* name)
  {
    Type_handle out (new Classification::Type (name));
    m_types.push_back (out);
    return out;
  }

  
  /*!
    \brief Adds a classification type.

    \param type the handle to the classification type that must be added.
   */
  void add_classification_type (Type_handle type)
  {
    m_types.push_back (type);
  }

  /*!
    \brief Removes a classification type.

    \param type the handle to the classification type that must be removed.

    \return `true` if the classification type was correctly removed,
    `false` if its handle was not found inside the object.
   */ 
 bool remove_classification_type (Type_handle type)
  {
    std::size_t idx = (std::size_t)(-1);
    for (std::size_t i = 0; i < m_types.size(); ++ i)
      if (m_types[i] == type)
        {
          m_types.erase (m_types.begin() + i);
          idx = i;
          break;
        }
    if (idx == (std::size_t)(-1))
      return false;
    std::cerr << idx << std::endl;
    
    for (std::size_t i = 0; i < m_assigned_type.size(); ++ i)
      if (m_assigned_type[i] == (std::size_t)(-1))
          continue;
      else if (m_assigned_type[i] > idx)
        m_assigned_type[i] --;
      else if (m_assigned_type[i] == idx)
        m_assigned_type[i] = (std::size_t)(-1);

    for (std::size_t i = 0; i < m_training_type.size(); ++ i)
      if (m_assigned_type[i] == (std::size_t)(-1))
        continue;
      else if (m_training_type[i] > idx)
        m_training_type[i] --;
      else if (m_training_type[i] == idx)
        m_training_type[i] = (std::size_t)(-1);

    return true;
  }

  /// \cond SKIP_IN_MANUAL
  std::size_t number_of_classification_types () const
  {
    return m_types.size();
  }

  Type_handle get_classification_type (std::size_t idx)
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

    \param attribute %Handle of the attribute to add.
   */
  void add_attribute (Attribute_handle attribute)
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

  Attribute_handle get_attribute(std::size_t idx)
  {
    return m_attributes[idx];
  }
  /// \endcond
  
  /// @}

  /// \name Output
  /// @{

  /*!
    \brief Gets the classification type of an indexed item.

    \note If classification was not performed (using `run()`,
    `run_with_local_smoothing()` or `run_with_graphcut()`), this
    function always returns the default `Type_handle`.

    \param index index of the input item

    \return handle to the classification type
  */
  Type_handle classification_type_of (std::size_t index) const
  {
    if (m_assigned_type.size() <= index
        || m_assigned_type[index] == (std::size_t)(-1))
      {
      return Type_handle();
      }
    return m_types[m_assigned_type[index]];
  }

  /// \cond SKIP_IN_MANUAL
  bool classification_prepared() const
  {
    return !(m_assigned_type.empty());
  }
  void set_classification_type_of (std::size_t index, Type_handle class_type)
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
    \brief Gets the confidence of the classification type of an indexed item.

    \note If classification was not performed (using `run()`,
    `run_with_local_smoothing()` or `run_with_graphcut()`), this
    function always returns 0.

    \param index index of the input item
    \return confidence ranging from 0 (not confident at all) to 1 (very confident).
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

    All the `Classification::Type` and `Classification::Attribute`
    necessary for the classification should have been registered in
    the object before running this function.

    Each classification type must be given a small set of user-defined
    inliers to provide the training algorithm with a ground truth.

    \param nb_tests number of tests to perform. Higher values may
    provide the user with better results at the cost of higher
    computation time. Using a value of at least 10 times the number of
    attributes is advised.

    \return minimum ratio (over all classification types) of provided
    ground truth items correctly classified using the best
    configuration found.
  */
  
  double train (std::size_t nb_tests = 300)
  {
    if (m_training_type.empty())
      return 0.;

    std::vector<std::vector<std::size_t> > training_sets (m_types.size());
    for (std::size_t i = 0; i < m_training_type.size(); ++ i)
      if (m_training_type[i] != (std::size_t)(-1))
        training_sets[m_training_type[i]].push_back (i);

    for (std::size_t i = 0; i < training_sets.size(); ++ i)
      if (training_sets[i].empty())
        std::cerr << "WARNING: \"" << m_types[i]->name() << "\" doesn't have a training set." << std::endl;

    std::vector<double> best_weights (m_attributes.size(), 1.);

    struct Attribute_training
    {
      bool skipped;
      double wmin;
      double wmax;
      double factor;
    };
    std::vector<Attribute_training> att_train;
    std::size_t nb_trials = 100;
    double wmin = 1e-5, wmax = 1e5;
    double factor = std::pow (wmax/wmin, 1. / (double)nb_trials);
    std::size_t att_used = 0;
    for (std::size_t j = 0; j < m_attributes.size(); ++ j)
      {
        Attribute_handle att = m_attributes[j];
        best_weights[j] = att->weight();

        std::size_t nb_useful = 0;
        double min = (std::numeric_limits<double>::max)();
        double max = -(std::numeric_limits<double>::max)();

        att->set_weight(wmin);
        for (std::size_t i = 0; i < 100; ++ i)
          {
            estimate_attribute_effect(training_sets, att);
            if (attribute_useful(att))
              {
                CGAL_CLASSTRAINING_CERR << "#";
                nb_useful ++;
                min = (std::min) (min, att->weight());
                max = (std::max) (max, att->weight());
              }
            else
              CGAL_CLASSTRAINING_CERR << "-";
            att->set_weight(factor * att->weight());
          }
        CGAL_CLASSTRAINING_CERR << std::endl;
        CGAL_CLASSTRAINING_CERR << att->name() << " useful in "
                  << nb_useful << "% of the cases, in interval [ "
                  << min << " ; " << max << " ]" << std::endl;
        att_train.push_back (Attribute_training());
        att_train.back().skipped = false;
        att_train.back().wmin = min / factor;
        att_train.back().wmax = max * factor;
        if (nb_useful < 2)
          {
            att_train.back().skipped = true;
            att->set_weight(0.);
            best_weights[j] = att->weight();
          }
        else if (best_weights[j] == 1.)
          {
            att->set_weight(0.5 * (att_train.back().wmin + att_train.back().wmax));
            best_weights[j] = att->weight();
            ++ att_used;
          }
        else
          {
            att->set_weight(best_weights[j]);
            ++ att_used;
          }
        estimate_attribute_effect(training_sets, att);
      }

    std::size_t nb_trials_per_attribute = 1 + (std::size_t)(nb_tests / (double)(att_used));
    std::cerr << "Trials = " << nb_tests << ", attributes = " << att_used
              << ", trials per att = " << nb_trials_per_attribute << std::endl;
    for (std::size_t i = 0; i < att_train.size(); ++ i)
      if (!(att_train[i].skipped))
        att_train[i].factor = std::pow (att_train[i].wmax / att_train[i].wmin,
                                        1. / (double)nb_trials_per_attribute);
    
    
    prepare_classification();
    
    double best_score = training_compute_worst_score(training_sets, 0.);
    double best_confidence = training_compute_worst_confidence(training_sets, 0.);
    
    std::cerr << "TRAINING GLOBALLY: Best score evolution: " << std::endl;

    std::cerr << 100. * best_score << "% (found at initialization)" << std::endl;

    std::size_t current_att_changed = 0;
    for (std::size_t i = 0; i < att_used; ++ i)
      {
        while (att_train[current_att_changed].skipped)
          {
            ++ current_att_changed;
            if (current_att_changed == m_attributes.size())
              current_att_changed = 0;
          }

        std::size_t nb_used = 0;
        for (std::size_t j = 0; j < m_attributes.size(); ++ j)
          {
            if (j == current_att_changed)
              continue;
            
            m_attributes[j]->set_weight(best_weights[j]);
            estimate_attribute_effect(training_sets, m_attributes[j]);
            if (attribute_useful(m_attributes[j]))
              nb_used ++;
          }
        Attribute_handle current_att = m_attributes[current_att_changed];
        const Attribute_training& tr = att_train[current_att_changed];
        
        current_att->set_weight(tr.wmin);
        for (std::size_t j = 0; j < nb_trials_per_attribute; ++ j)
          {
            estimate_attribute_effect(training_sets, current_att);

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
                          << (i * nb_trials_per_attribute) + j << "/" << nb_tests << ", "
                          << nb_used + (attribute_useful(current_att) ? 1 : 0)
                          << "/" << m_attributes.size() << " attribute(s) used)" << std::endl;
                for (std::size_t k = 0; k < m_attributes.size(); ++ k)
                  {
                    Attribute_handle att = m_attributes[k];
                    best_weights[k] = att->weight();
                  }
              }
            
            current_att->set_weight(current_att->weight() * tr.factor);
          }

        ++ current_att_changed;
      }

    for (std::size_t i = 0; i < best_weights.size(); ++ i)
      {
        Attribute_handle att = m_attributes[i];
        att->set_weight(best_weights[i]);
      }

    estimate_attributes_effects(training_sets);
    
    std::cerr << std::endl << "Best score found is at least " << 100. * best_score
              << "% of correct classification" << std::endl;

    std::size_t nb_removed = 0;
    for (std::size_t i = 0; i < best_weights.size(); ++ i)
      {
        Attribute_handle att = m_attributes[i];
        CGAL_CLASSTRAINING_CERR << "ATTRIBUTE " << att->name() << ": " << best_weights[i] << std::endl;
        att->set_weight(best_weights[i]);

        Classification::Attribute::Effect side = m_types[0]->attribute_effect(att);
        bool to_remove = true;
        for (std::size_t j = 0; j < m_types.size(); ++ j)
          {
            Type_handle ctype = m_types[j];
            if (ctype->attribute_effect(att) == Classification::Attribute::FAVORING)
              CGAL_CLASSTRAINING_CERR << " * Favored for ";
            else if (ctype->attribute_effect(att) == Classification::Attribute::PENALIZING)
              CGAL_CLASSTRAINING_CERR << " * Penalized for ";
            else
              CGAL_CLASSTRAINING_CERR << " * Neutral for ";
            if (ctype->attribute_effect(att) != side)
              to_remove = false;
            CGAL_CLASSTRAINING_CERR << ctype->name() << std::endl;
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
    \brief Adds the input item specified by index `idx` as an inlier
    of `class_type` for the training algorithm.

    \param class_type handle to the classification type.

    \param idx index of the input item.
  */
  bool set_inlier (Type_handle class_type, std::size_t idx)
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
    \brief Adds input items specified by a range of indices as
    inliers of `class_type` for the training algorithm.

    \param class_type handle to the classification type.
    \param indices Set of incides to add as inliers.

    \tparam IndexRange range of `std::size_t`, model of `ConstRange`.
  */
  template <class IndexRange>
  bool add_training_set (Type_handle class_type,
                         IndexRange indices)
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

    for (typename IndexRange::const_iterator it = indices.begin();
         it != indices.end(); ++ it)
      m_training_type[*it] = type_idx;

    return true;
  }
  
  /// @}

  
  /// \cond SKIP_IN_MANUAL
  Type_handle training_type_of (std::size_t index) const
  {
    if (m_training_type.size() <= index
        || m_training_type[index] == (std::size_t)(-1))
      return Type_handle();
    return m_types[m_training_type[index]];
  }

  void prepare_classification ()
  {
    // Reset data structure
    std::vector<std::size_t>(m_input.size(), (std::size_t)(-1)).swap (m_assigned_type);
    std::vector<double>(m_input.size()).swap (m_confidence);

    m_effect_table = std::vector<std::vector<Attribute_effect> >
      (m_types.size(), std::vector<Attribute_effect> (m_attributes.size(),
                                                                 Classification::Attribute::NEUTRAL));
    
    for (std::size_t i = 0; i < m_effect_table.size (); ++ i)
      for (std::size_t j = 0; j < m_effect_table[i].size (); ++ j)
        m_effect_table[i][j] = m_types[i]->attribute_effect (m_attributes[j]);

  }

  /// \endcond


private:

  double classification_value (const std::size_t& class_type,
                               const std::size_t& pt_index) const
  {
    double out = 0.;
    for (std::size_t i = 0; i < m_effect_table[class_type].size(); ++ i)
      {
        if (m_attributes[i]->weight() == 0.)
          continue;
        if (m_effect_table[class_type][i] == Classification::Attribute::FAVORING)
          out += m_attributes[i]->favored (pt_index);
        else if (m_effect_table[class_type][i] == Classification::Attribute::PENALIZING)
          out += m_attributes[i]->penalized (pt_index);
        else if (m_effect_table[class_type][i] == Classification::Attribute::NEUTRAL)
          out += m_attributes[i]->ignored (pt_index);
      }
    return out;
  }


  void estimate_attributes_effects
  (const std::vector<std::vector<std::size_t> >& training_sets)
  {
    for (std::size_t i = 0; i < m_attributes.size(); ++ i)
      estimate_attribute_effect (training_sets, m_attributes[i]);
  }

  void estimate_attribute_effect
  (const std::vector<std::vector<std::size_t> >& training_sets,
   Attribute_handle att)
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
        Type_handle ctype = m_types[j];
            
        for (std::size_t k = 0; k < training_sets[j].size(); ++ k)
          {
            double val = att->normalized(training_sets[j][k]);
            sd[j] += (val - mean[j]) * (val - mean[j]);
          }
        sd[j] = std::sqrt (sd[j] / training_sets[j].size());
        if (mean[j] - sd[j] > 0.5)
          ctype->set_attribute_effect (att, Classification::Attribute::FAVORING);
        else if (mean[j] + sd[j] < 0.5)
          ctype->set_attribute_effect (att, Classification::Attribute::PENALIZING);
        else
          ctype->set_attribute_effect (att, Classification::Attribute::NEUTRAL);
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

  bool attribute_useful (Attribute_handle att)
  {
    Classification::Attribute::Effect side = m_types[0]->attribute_effect(att);
    for (std::size_t k = 1; k < m_types.size(); ++ k)
      if (m_types[k]->attribute_effect(att) != side)
        return true;
    return false;
  }

};





} // namespace CGAL

#endif // CGAL_CLASSIFIER_H

