// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
// Copyright (c) 2017 GeometryFactory Sarl (France).
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
#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Classification/Label.h>

#include <CGAL/Memory_sizer.h>
#include <CGAL/Timer.h>

#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>
#include <tbb/mutex.h>
#endif // CGAL_LINKED_WITH_TBB

//#define CGAL_CLASSIFICATION_VERBOSE
#if defined(CGAL_CLASSIFICATION_VERBOSE)
#define CGAL_CLASSIFICATION_CERR std::cerr
#else
#define CGAL_CLASSIFICATION_CERR std::ostream(0)
#endif

namespace CGAL {

/*!
\ingroup PkgClassification

\brief Classifies a data set based on a set of features and a set of
labels.

This class implements the core of the classification algorithm
\cgalCite{cgal:lm-clscm-12} (section 2). It uses a data set as input
and assigns each input item to a label among a set of user defined
labels. To achieve this classification, a set of local geometric
features are used, such as planarity, elevation or vertical
dispersion. In addition, the user must define a set of labels such as
building, ground or vegetation.

Each pair of feature and label must be assigned an
[Feature::Effect](@ref CGAL::Classification::Feature::Effect) (for
example, vegetation has a low planarity and a high vertical
dispersion) and each feature must be assigned a weight. These
parameters can be set up by hand or by automatic training, provided a
small user-defined set of inlier is given for each classification
label.

\tparam ItemRange model of `ConstRange`. Its iterator type is
`RandomAccessIterator`.

\tparam ItemMap model of `ReadablePropertyMap` whose key
type is the value type of the iterator of `ItemRange` and value type is
the type of the items that are classified.

\tparam ConcurrencyTag enables sequential versus parallel
algorithm. Possible values are `Parallel_tag` (default value is CGAL
is linked with TBB) or `Sequential_tag` (default value otherwise).
*/
template <typename ItemRange,
          typename ItemMap,
#if defined(DOXYGEN_RUNNING)
          typename ConcurrencyTag>
#elif defined(CGAL_LINKED_WITH_TBB)
          typename ConcurrencyTag = CGAL::Parallel_tag>
#else
          typename ConcurrencyTag = CGAL::Sequential_tag>
#endif

class Classifier
{

  
public:
  typedef typename Classification::Label_handle   Label_handle;
  typedef typename Classification::Feature_handle Feature_handle;
  typedef Classification::Feature::Effect         Feature_effect;
  
  /// \cond SKIP_IN_MANUAL
  typedef typename ItemMap::value_type Item;

#ifdef CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
  typedef internal::Alpha_expansion_graph_cut_boost             Alpha_expansion;
#else
  typedef internal::Alpha_expansion_graph_cut_boykov_kolmogorov Alpha_expansion;
#endif
  
protected:

#ifdef CGAL_LINKED_WITH_TBB
  class Run
  {
    Classifier& m_classifier;
    const std::vector<std::vector<Feature_effect> >& m_effect_table;
    std::vector<std::size_t>& m_assigned_label;
    std::vector<double>& m_confidence;
    
    
  public:

    Run (Classifier& classifier,
         const std::vector<std::vector<Feature_effect> >& effect_table,
         std::vector<std::size_t>& assigned_label,
         std::vector<double>& confidence)
      : m_classifier (classifier), m_effect_table (effect_table),
        m_assigned_label (assigned_label), m_confidence (confidence)
    { }
          
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t s = r.begin(); s != r.end(); ++ s)
        {
          std::size_t nb_class_best=0; 

          double val_class_best = (std::numeric_limits<double>::max)();
          std::vector<double> values;
      
          for(std::size_t k = 0; k < m_effect_table.size(); ++ k)
            {
              double value = m_classifier.classification_value (k, s);
              values.push_back (value);
          
              if(val_class_best > value)
                {
                  val_class_best = value;
                  nb_class_best=k;
                }
            }

          m_assigned_label[s] = nb_class_best;

          std::sort (values.begin(), values.end());
          m_confidence[s] = values[1] - values[0];
        }
    }

  };

  class Run_with_local_smoothing_preprocessing
  {
    Classifier& m_classifier;
    std::vector<std::vector<double> >& m_values;
    const std::vector<std::vector<Feature_effect> >& m_effect_table;
    
  public:

    Run_with_local_smoothing_preprocessing (Classifier& classifier,
                                            std::vector<std::vector<double> >& values,
                                            const std::vector<std::vector<Feature_effect> >& effect_table)
      : m_classifier (classifier), m_values(values), m_effect_table (effect_table)
    { }
          
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t s = r.begin(); s != r.end(); ++ s)
        for(std::size_t k = 0; k < m_effect_table.size(); ++ k)
          m_values[k][s] = m_classifier.classification_value (k, s);
    }
  };
  
  template <typename NeighborQuery>
  class Run_with_local_smoothing
  {
    Classifier& m_classifier;
    const ItemRange& m_input;
    ItemMap m_item_map;
    const std::vector<std::vector<double> >& m_values;
    const NeighborQuery& m_neighbor_query;
    std::vector<std::size_t>& m_assigned_label;
    std::vector<double>& m_confidence;
    
  public:

    Run_with_local_smoothing (Classifier& classifier,
                              const ItemRange& input,
                              ItemMap item_map,
                              const std::vector<std::vector<double> >& values,
                              const NeighborQuery& neighbor_query,
                              std::vector<std::size_t>& assigned_label,
                              std::vector<double>& confidence)
    : m_classifier (classifier), m_input (input), m_item_map (item_map),
      m_values(values),
      m_neighbor_query (neighbor_query),
      m_assigned_label (assigned_label), m_confidence (confidence)
    { }
          
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t s = r.begin(); s != r.end(); ++ s)
        {
          std::vector<std::size_t> neighbors;
          m_neighbor_query (get (m_item_map, *(m_input.begin()+s)), std::back_inserter (neighbors));

          std::vector<double> mean (m_values.size(), 0.);
          for (std::size_t n = 0; n < neighbors.size(); ++ n)
            for (std::size_t j = 0; j < m_values.size(); ++ j)
              mean[j] += m_values[j][neighbors[n]];

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

          m_assigned_label[s] = nb_class_best;

          std::sort (mean.begin(), mean.end());
          m_confidence[s] = mean[1] - mean[0];      
        }
    }

  };

  tbb::mutex m_mutex;
  void mutex_lock() { m_mutex.lock(); }
  void mutex_unlock() { m_mutex.unlock(); }
#else // CGAL_LINKED_WITH_TBB
  void mutex_lock() { }
  void mutex_unlock() { }
#endif // CGAL_LINKED_WITH_TBB
  
  const ItemRange& m_input;
  ItemMap m_item_map;

  std::vector<std::size_t> m_assigned_label;
  std::vector<double> m_confidence;

  std::vector<Label_handle> m_labels; 
  std::vector<Feature_handle> m_features; 

  std::vector<std::vector<Feature_effect> > m_effect_table;
  /// \endcond

public:


  /// \name Constructor
  /// @{
  
  /*! 
    \brief Initializes a classification object.

    \param input input range.

    \param item_map property map to access the input items.
  */
  Classifier (const ItemRange& input,
              ItemMap item_map)
    : m_input (input), m_item_map (item_map)
  {
  }

  /// @}

  /// \cond SKIP_IN_MANUAL
  virtual ~Classifier() { }
  /// \endcond
  
  /// \name Features
  /// @{

  /*!
    \brief Adds a feature.

    \tparam Feature type of the feature, inherited from
    `Classification::Feature_base`.

    \tparam T types of the parameters of the feature's constructor
    (with the exception of the first parameter that is always of type
    `ItemRange&` and that is automatically passed by the classifier to
    the feature's constructor).

    \param t parameters of the feature's constructor (with the
    exception of the first parameter that is always the input item
    range and that is automatically passed by the classifier to the
    feature's constructor).

    \return a handle to the newly added feature.
   */
#if (!defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE)) || DOXYGEN_RUNNING
  template <typename Feature, typename ... T>
  Feature_handle add_feature (T&& ... t)
  {
    Feature_handle fh (new Feature(m_input, std::forward<T>(t)...));
    mutex_lock();
    m_features.push_back (fh);
    mutex_unlock();
    return fh;
  }
#else
  template <typename Feature>
  Feature_handle add_feature ()
  {
    Feature_handle fh (new Feature(m_input));
    mutex_lock();
    m_features.push_back (fh);
    mutex_unlock();
    return fh;
  }
  template <typename Feature, typename T1>
  Feature_handle add_feature (T1& t1)
  {
    Feature_handle fh (new Feature(m_input, t1));
    mutex_lock();
    m_features.push_back (fh);
    mutex_unlock();
    return fh;
  }
  template <typename Feature, typename T1, typename T2>
  Feature_handle add_feature (T1& t1, T2& t2)
  {
    Feature_handle fh (new Feature(m_input, t1, t2));
    mutex_lock();
    m_features.push_back (fh);
    mutex_unlock();
    return fh;
  }
  template <typename Feature, typename T1, typename T2, typename T3>
  Feature_handle add_feature (T1& t1, T2& t2, T3& t3)
  {
    Feature_handle fh (new Feature(m_input, t1, t2, t3));
    mutex_lock();
    m_features.push_back (fh);
    mutex_unlock();
    return fh;
  }
  template <typename Feature, typename T1, typename T2, typename T3, typename T4>
  Feature_handle add_feature (T1& t1, T2& t2, T3& t3, T4& t4)
  {
    Feature_handle fh (new Feature(m_input, t1, t2, t3, t4));
    mutex_lock();
    m_features.push_back (fh);
    mutex_unlock();
    return fh;
  }
  template <typename Feature, typename T1, typename T2, typename T3, typename T4, typename T5>
  Feature_handle add_feature (T1& t1, T2& t2, T3& t3, T4& t4, T5& t5)
  {
    Feature_handle fh (new Feature(m_input, t1, t2, t3, t4, t5));
    mutex_lock();
    m_features.push_back (fh);
    mutex_unlock();
    return fh;
  }
#endif

  /*!
    \brief Removes a feature.

    \param feature the handle to feature type that must be removed.

    \return `true` if the feature was correctly removed, `false` if
    its handle was not found.
   */ 
  bool remove_feature (Feature_handle feature)
  {
    for (std::size_t i = 0; i < m_features.size(); ++ i)
      if (m_features[i] == feature)
        {
          m_features.erase (m_features.begin() + i);
          return true;
        }
    return false;
  }

  /*!
    \brief Returns how many features are defined.
  */  
  std::size_t number_of_features() const
  {
    return m_features.size();
  }


  /*!
    \brief Returns the \f$i^{th}\f$ feature.
  */  
  Feature_handle feature(std::size_t i)
  {
    return m_features[i];
  }

  /*!
    \brief Removes all features.
   */
  void clear_features ()
  {
    m_features.clear();
  }

  /// \endcond
  
  /// @}


  /// \name Labels
  /// @{
  
  /*!
    \brief Adds a label.

    \param name name of the label.

    \return a handle to the newly added label.
   */
  Label_handle add_label (const char* name)
  {
    Label_handle out (new Classification::Label (name));
    m_labels.push_back (out);
    return out;
  }

  /*!
    \brief Removes a label.

    \param label the handle to the label that must be removed.

    \return `true` if the label was correctly removed,
    `false` if its handle was not found.
   */ 
 bool remove_label (Label_handle label)
  {
    std::size_t idx = (std::size_t)(-1);
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (m_labels[i] == label)
        {
          m_labels.erase (m_labels.begin() + i);
          idx = i;
          break;
        }
    if (idx == (std::size_t)(-1))
      return false;
    std::cerr << idx << std::endl;
    
    for (std::size_t i = 0; i < m_assigned_label.size(); ++ i)
      if (m_assigned_label[i] == (std::size_t)(-1))
          continue;
      else if (m_assigned_label[i] > idx)
        m_assigned_label[i] --;
      else if (m_assigned_label[i] == idx)
        m_assigned_label[i] = (std::size_t)(-1);

    return true;
  }

  /*!
    \brief Returns how many labels are defined.
  */  
  std::size_t number_of_labels () const
  {
    return m_labels.size();
  }
  
  /*!
    \brief Returns the \f$i^{th}\f$ label.
  */  
  Label_handle label (std::size_t i) const
  {
    return m_labels[i];
  }


  /*!
    \brief Removes all labels.
   */
  void clear_labels ()
  {
    m_labels.clear();
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
    
#ifndef CGAL_LINKED_WITH_TBB
    CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                               "Parallel_tag is enabled but TBB is unavailable.");
#else
    if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
      {
        Run f (*this, m_effect_table, m_assigned_label, m_confidence);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_input.size ()), f);
      }
    else
#endif
      {

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

            m_assigned_label[s] = nb_class_best;

            std::sort (values.begin(), values.end());
            m_confidence[s] = values[1] - values[0];
          }
      }
  }


  /*! 
    \brief Runs the classification algorithm with a local smoothing.

    The computed classification energy is smoothed on a user defined
    local neighborhood of items. This method is a compromise between
    efficiency and reliability.

    \tparam NeighborQuery model of `NeighborQuery`.
    \param neighbor_query used to access neighborhoods of items.
  */
  template <typename NeighborQuery>
  void run_with_local_smoothing (const NeighborQuery& neighbor_query)
  {
    prepare_classification ();
    
    // data term initialisation
    CGAL_CLASSIFICATION_CERR << "Labeling... ";

    std::vector<std::vector<double> > values
      (m_labels.size(),
       std::vector<double> (m_input.size(), -1.));

#ifndef CGAL_LINKED_WITH_TBB
    CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                               "Parallel_tag is enabled but TBB is unavailable.");
#else
    if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
      {
        Run_with_local_smoothing_preprocessing f1
          (*this, values, m_effect_table);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_input.size ()), f1);
        
        Run_with_local_smoothing<NeighborQuery> f2
          (*this, m_input, m_item_map, values, neighbor_query,m_assigned_label, m_confidence);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_input.size ()), f2);
      }
    else
#endif
      {
        for (std::size_t s=0; s < m_input.size(); ++ s)
          {
            std::vector<std::size_t> neighbors;
            neighbor_query (get (m_item_map, *(m_input.begin()+s)), std::back_inserter (neighbors));

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

            m_assigned_label[s] = nb_class_best;

            std::sort (mean.begin(), mean.end());
            m_confidence[s] = mean[1] - mean[0];      
          }
      }
  }


  /*! 
    \brief Runs the classification algorithm with a global
    regularization based on a graphcut.

    The computed classification energy is globally regularized through
    an alpha-expansion algorithm. This method is slow but provides
    the user with good quality results.

    \tparam NeighborQuery model of `NeighborQuery`.
    \param neighbor_query used to access neighborhoods of items.
    \param weight weight of the regularization with respect to the
    classification energy. Higher values produce more regularized
    output but may result in a loss of details.

  */
  template <typename NeighborQuery>
  void run_with_graphcut (const NeighborQuery& neighbor_query,
                          const double weight)
  {
    prepare_classification ();
    
    // data term initialisation
#ifdef CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
    CGAL_CLASSIFICATION_CERR << "Labeling using Boost with regularization weight " << weight << "... ";
#else
    CGAL_CLASSIFICATION_CERR << "Labeling using Boyvok Kolmogorov with regularization weight " << weight << "... ";
#endif

    std::vector<std::pair<std::size_t, std::size_t> > edges;
    std::vector<double> edge_weights;
    std::vector<std::vector<double> > probability_matrix
      (m_effect_table.size(), std::vector<double>(m_input.size(), 0.));
    std::vector<std::size_t>(m_input.size()).swap(m_assigned_label);

    std::cerr << "Size of probability matrix = " << m_effect_table.size() * m_input.size() << std::endl;
    std::cerr << "Size of assigned labels = " << m_assigned_label.size() << std::endl;
    
    for (std::size_t s = 0; s < m_input.size(); ++ s)
      {
        std::vector<std::size_t> neighbors;

        neighbor_query (get(m_item_map, *(m_input.begin()+s)), std::back_inserter (neighbors));

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
        m_assigned_label[s] = nb_class_best;
      }
    
    Alpha_expansion graphcut;
    graphcut(edges, edge_weights, probability_matrix, m_assigned_label);
    std::cerr << ((double)(CGAL::Memory_sizer().virtual_size()) / 1073741824.) << " GB allocated" << std::endl;
  }
  
  /// @}
  


  /// \name Output
  /// @{

  /*!
    \brief Returns the value of the energy of `label` at the item at
    position `index`.
  */
  double energy_of (Label_handle label, std::size_t index) const
  {
    double out = 0.;
    for (std::size_t i = 0; i < m_features.size(); ++ i)
      {
        if (m_features[i]->weight() == 0.)
          continue;

        Feature_effect eff = label->feature_effect (m_features[i]);

        if (eff == Classification::Feature::FAVORING)
          out += m_features[i]->favored (index);
        else if (eff == Classification::Feature::PENALIZING)
          out += m_features[i]->penalized (index);
        else if (eff == Classification::Feature::NEUTRAL)
          out += m_features[i]->ignored (index);
      }
    return out;
  }

  /*!
    \brief Returns the label of the item at position
    `index`.

    \note If classification was not performed (using `run()`,
    `run_with_local_smoothing()` or `run_with_graphcut()`), this
    function always returns the default `Label_handle`.
  */
  Label_handle label_of (std::size_t index) const
  {
    if (m_assigned_label.size() <= index
        || m_assigned_label[index] == (std::size_t)(-1))
      {
      return Label_handle();
      }
    return m_labels[m_assigned_label[index]];
  }

  /// \cond SKIP_IN_MANUAL
  void set_label_of (std::size_t index, Label_handle label)
  {
    if (index >= m_assigned_label.size())
      m_assigned_label.resize (index + 1, (std::size_t)(-1));
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (m_labels[i] == label)
        {
          m_assigned_label[index] = i;
          return;
        }
    m_assigned_label[index] = (std::size_t)(-1);
  }
  /// \endcond

  /*!

    \brief Returns the confidence of the label of the
    item at position `index`.

    \note If classification was not performed (using `run()`,
    `run_with_local_smoothing()` or `run_with_graphcut()`), this
    function always returns 0.

    \return confidence ranging from 0 (not confident at all) to 1
    (very confident).
  */
  double confidence_of (std::size_t index) const
  {
    if (m_confidence.size() <= index)
      return 0.;
    return m_confidence[index];
  }

  /// @}


  

protected:

  /// \cond SKIP_IN_MANUAL
  void prepare_classification ()
  {
    // Reset data structure
    std::vector<std::size_t>(m_input.size(), (std::size_t)(-1)).swap (m_assigned_label);
    std::vector<double>(m_input.size()).swap (m_confidence);

    m_effect_table = std::vector<std::vector<Feature_effect> >
      (m_labels.size(), std::vector<Feature_effect> (m_features.size(),
                                                                 Classification::Feature::NEUTRAL));
    
    for (std::size_t i = 0; i < m_effect_table.size (); ++ i)
      for (std::size_t j = 0; j < m_effect_table[i].size (); ++ j)
        m_effect_table[i][j] = m_labels[i]->feature_effect (m_features[j]);

  }
  
  double classification_value (const std::size_t& label,
                               const std::size_t& pt_index) const
  {
    double out = 0.;
    for (std::size_t i = 0; i < m_effect_table[label].size(); ++ i)
      {
        if (m_features[i]->weight() == 0.)
          continue;
        if (m_effect_table[label][i] == Classification::Feature::FAVORING)
          out += m_features[i]->favored (pt_index);
        else if (m_effect_table[label][i] == Classification::Feature::PENALIZING)
          out += m_features[i]->penalized (pt_index);
        else if (m_effect_table[label][i] == Classification::Feature::NEUTRAL)
          out += m_features[i]->ignored (pt_index);
      }
    return out;
  }
  /// \endcond
};





} // namespace CGAL

#endif // CGAL_CLASSIFIER_H

