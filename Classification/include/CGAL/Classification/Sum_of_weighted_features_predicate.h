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

#ifndef CLASSIFICATION_SUM_OF_WEIGHTED_FEATURES_PREDICATE_H
#define CLASSIFICATION_SUM_OF_WEIGHTED_FEATURES_PREDICATE_H

#include <CGAL/Classification/Feature_set.h>
#include <CGAL/Classification/Label_set.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

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

//#define CGAL_CLASSTRAINING_VERBOSE
#if defined(CGAL_CLASSTRAINING_VERBOSE)
#define CGAL_CLASSTRAINING_CERR std::cerr
#else
#define CGAL_CLASSTRAINING_CERR std::ostream(0)
#endif

namespace CGAL {

namespace Classification {

class Sum_of_weighted_features_predicate
{
public:

  enum Effect /// Defines the effect of an feature on a type.
    {
      FAVORING = 0, ///< High values of the feature favor this type
      NEUTRAL = 1, ///< The feature has no effect on this type
      PENALIZING = 2 ///< Low values of the feature favor this type
    };
  
private:

#ifdef CGAL_LINKED_WITH_TBB
  class Compute_worst_score_and_confidence
  {
    std::vector<std::size_t>& m_training_set;
    const Sum_of_weighted_features_predicate& m_predicate;
    std::size_t m_label;
    float& m_confidence;
    std::size_t& m_nb_okay;
    tbb::mutex& m_mutex;
    
  public:

    Compute_worst_score_and_confidence (std::vector<std::size_t>& training_set,
                                        const Sum_of_weighted_features_predicate& predicate,
                                        std::size_t label,
                                        float& confidence,
                                        std::size_t& nb_okay,
                                        tbb::mutex& mutex)
      : m_training_set (training_set)
      , m_predicate (predicate)
      , m_label (label)
      , m_confidence (confidence)
      , m_nb_okay (nb_okay)
      , m_mutex (mutex)
    { }
    
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t k = r.begin(); k != r.end(); ++ k)
        {
          std::vector<std::pair<float, std::size_t> > values;
          
          std::vector<float> v;
          m_predicate.probabilities (m_training_set[k], v);
  
          for(std::size_t l = 0; l < v.size(); ++ l)
            values.push_back (std::make_pair (v[l], l));
          
          std::sort (values.begin(), values.end());

          if (values[0].second == m_label)
            {
              m_mutex.lock();
              m_confidence += values[1].first - values[0].first;
              ++ m_nb_okay;
              m_mutex.unlock();
            }
        }
    }

  };
#endif // CGAL_LINKED_WITH_TBB


  Label_set& m_labels;
  Feature_set& m_features;
  std::vector<float> m_weights;
  std::vector<std::vector<Effect> > m_effect_table;
  mutable std::map<Label_handle, std::size_t> m_map_labels;
  mutable std::map<Feature_handle, std::size_t> m_map_features;

public:
  Sum_of_weighted_features_predicate (Label_set& labels,
                                      Feature_set& features)
    : m_labels (labels), m_features (features),
      m_weights (features.size(), 1.),
      m_effect_table (labels.size(), std::vector<Effect>
                      (features.size(),
                        NEUTRAL))
  {
    for (std::size_t i = 0; i < labels.size(); ++ i)
      m_map_labels[labels[i]] = i;
    for (std::size_t i = 0; i < features.size(); ++ i)
      m_map_features[features[i]] = i;
  }

  /*!
    \brief Sets the weight of the feature (`weight` must be positive).
  */
  void set_weight (Feature_handle feature, float weight)
  {
    m_weights[m_map_features[feature]] = weight;
  }
  /// \cond SKIP_IN_MANUAL
  void set_weight (std::size_t feature, float weight)
  {
    m_weights[feature] = weight;
  }
  /// \endcond

  /*!
    \brief Returns the weight of the feature.
  */
  float weight (Feature_handle feature) const
  {
    return m_weights[m_map_features[feature]];
  }
  /// \cond SKIP_IN_MANUAL
  float weight (std::size_t feature) const
  {
    return m_weights[feature];
  }
  /// \endcond

  /*! 
    \brief Sets the `effect` of `feature` on `label`.
  */ 
  void set_effect (Label_handle label, Feature_handle feature,
                   Effect effect)
  {
    m_effect_table[m_map_labels[label]][m_map_features[feature]] = effect;
  }
  /// \cond SKIP_IN_MANUAL
  void set_effect (std::size_t label, std::size_t feature,
                   Effect effect)
  {
    m_effect_table[label][feature] = effect;
  }
  /// \endcond

  /*! 
    \brief Returns the `effect` of `feature` on `label`.
  */ 
  Effect effect (Label_handle label, Feature_handle feature) const
  {
    return m_effect_table[m_map_labels[label]][m_map_features[feature]];
  }
  /// \cond SKIP_IN_MANUAL
  Effect effect (std::size_t label, std::size_t feature) const
  {
    return m_effect_table[label][feature];
  }
  /// \endcond
  
  void probabilities (std::size_t item_index,
                      std::vector<float>& out) const
  {
    out.resize (m_labels.size());
    for (std::size_t l = 0; l < m_labels.size(); ++ l)
      {
        out[l] = 0.;
        for (std::size_t f = 0; f < m_features.size(); ++ f)
          if (weight(f) != 0.)
            out[l] += value (l, f, item_index);
      }
  }


  /*!
    \brief Saves the current configuration in the stream `output`.

    This allows to easily save and recover a specific classification
    configuration, that is to say:

    - The size of the smallest scale
    - The features and their respective weights
    - The labels and the effects of the features on them

    The output file is written in an XML format that is readable by
    the `load_configuration()` method.
  */
  void save_configuration (std::ostream& output)
  {
    boost::property_tree::ptree tree;

    for (std::size_t i = 0; i < m_features.size(); ++ i)
      {
        if (weight(m_features[i]) == 0)
          continue;
        boost::property_tree::ptree ptr;
        
        ptr.put("name", m_features[i]->name());
        ptr.put("weight", weight(m_features[i]));
        tree.add_child("classification.features.feature", ptr);
      }


    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      {
        boost::property_tree::ptree ptr;
        ptr.put("name", m_labels[i]->name());
        for (std::size_t j = 0; j < m_features.size(); ++ j)
          {
            if (weight(j) == 0)
              continue;
            boost::property_tree::ptree ptr2;
            ptr2.put("name", m_features[j]->name());
            Effect e = effect(i, j);
            if (e == PENALIZING)
              ptr2.put("effect", "penalized");
            else if (e == NEUTRAL)
              ptr2.put("effect", "neutral");
            else if (e == FAVORING)
              ptr2.put("effect", "favored");
            ptr.add_child("feature", ptr2);
          }
        tree.add_child("classification.labels.label", ptr);
      }

    // Write property tree to XML file
    boost::property_tree::xml_writer_settings<std::string> settings(' ', 3);
    boost::property_tree::write_xml(output, tree, settings);
  }
  
  /*!
    \brief Loads a configuration from the stream `input`.

    All data structures, features and labels specified in the input
    stream `input` are instantiated if possible (in particular,
    property maps needed should be provided), similarly to what is
    done in `generate_features()`.

    The input file should be in the XML format written by the
    `save_configuration()` method.

    \tparam VectorMap model of `ReadablePropertyMap` whose key type is
    the value type of the iterator of `PointRange` and value type is
    `Geom_traits::Vector_3`.
    \tparam ColorMap model of `ReadablePropertyMap`  whose key type is
    the value type of the iterator of `PointRange` and value type is
    `CGAL::Classification::RGB_Color`.
    \tparam EchoMap model of `ReadablePropertyMap` whose key type is
    the value type of the iterator of `PointRange` and value type is
    `std::size_t`.
    \param input input stream.
    \param normal_map property map to access the normal vectors of the input points (if any).
    \param color_map property map to access the colors of the input points (if any).
    \param echo_map property map to access the echo values of the input points (if any).
  */
  bool load_configuration (std::istream& input, bool verbose = false)
  {
    bool out = true;
    std::map<std::string, std::size_t> map_n2l;
    std::map<std::string, std::size_t> map_n2f;
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      map_n2l.insert (std::make_pair (m_labels[i]->name(), i));
    for (std::size_t i = 0; i < m_features.size(); ++ i)
      map_n2f.insert (std::make_pair (m_features[i]->name(), i));

    boost::property_tree::ptree tree;
    boost::property_tree::read_xml(input, tree);

    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, tree.get_child("classification.features"))
      {
        std::string name = v.second.get<std::string>("name");
        typename std::map<std::string, std::size_t>::iterator
          found = map_n2f.find (name);
        if (found != map_n2f.end())
          m_weights[found->second] = v.second.get<float>("weight");
        else
          {
            if (verbose)
              std::cerr << "Warning: feature \"" << name << "\" in configuration file not found" << std::endl;
            out = false;
          }
      }

    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, tree.get_child("classification.labels"))
      {
        std::string label_name = v.second.get<std::string>("name");
        typename std::map<std::string, std::size_t>::iterator
          found = map_n2l.find (label_name);
        std::size_t l = 0;
        if (found != map_n2l.end())
          l = found->second;
        else
          {
            if (verbose)
              std::cerr << "Warning: label \"" << label_name << "\" in configuration file not found" << std::endl;
            out = false;
            continue;
          }
        
        BOOST_FOREACH(boost::property_tree::ptree::value_type &v2, v.second)
          {
            if (v2.first == "name")
              continue;
            
            std::string feature_name = v2.second.get<std::string>("name");
            
            typename std::map<std::string, std::size_t>::iterator
              found2 = map_n2f.find (feature_name);
            std::size_t f = 0;
            if (found2 != map_n2f.end())
              f = found2->second;
            else if (verbose)
              {
                if (verbose)
                  std::cerr << "Warning: feature \"" << feature_name << "\" in configuration file not found" << std::endl;
                out = false;
                continue;
              }
            std::string e = v2.second.get<std::string>("effect");
            if (e == "penalized")
              set_effect (l, f, PENALIZING);
            else if (e == "neutral")
              set_effect (l, f, NEUTRAL);
            else
              set_effect (l, f, FAVORING);
          }
      }
    return out;
  }
  /// \name Training
  /// @{

  /*!
    \brief Runs the training algorithm.

    All the `Classification::Label` and `Classification::Feature`
    necessary for classification should have been added before running
    this function. After training, the user can call `run()`,
    `run_with_local_smoothing()` or `run_with_graphcut()` to compute
    the classification using the estimated parameters.

    \param nb_tests number of tests to perform. Higher values may
    provide the user with better results at the cost of a higher
    computation time. Using a value of at least 10 times the number of
    features is advised.

    \return minimum ratio (over all labels) of provided
    ground truth items correctly classified using the best
    configuration found.
  */
  template <typename ConcurrencyTag>  
  float train (const std::vector<std::size_t>& ground_truth,
                std::size_t nb_tests = 300)
  {
    std::vector<std::vector<std::size_t> > training_sets (m_labels.size());
    std::size_t nb_tot = 0;
    for (std::size_t i = 0; i < ground_truth.size(); ++ i)
      if (ground_truth[i] != std::size_t(-1))
        {
          training_sets[ground_truth[i]].push_back (i);
          ++ nb_tot;
        }

    CGAL_CLASSIFICATION_CERR << "Training using " << nb_tot << " inliers" << std::endl;
    
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (training_sets.size() <= i || training_sets[i].empty())
        std::cerr << "WARNING: \"" << m_labels[i]->name() << "\" doesn't have a training set." << std::endl;

    std::vector<float> best_weights (m_features.size(), 1.);

    struct Feature_training
    {
      std::size_t i;
      float wmin;
      float wmax;
      float factor;

      bool operator<(const Feature_training& other) const
      {
        return (wmin / wmax) < (other.wmin / other.wmax);
      }
    };
    std::vector<Feature_training> feature_train;
    std::size_t nb_trials = 100;
    float wmin = 1e-5, wmax = 1e5;
    float factor = std::pow (wmax/wmin, 1. / (float)nb_trials);
    
    for (std::size_t j = 0; j < m_features.size(); ++ j)
      {
        Feature_handle feature = m_features[j];
        best_weights[j] = weight(j);

        std::size_t nb_useful = 0;
        float min = (std::numeric_limits<float>::max)();
        float max = -(std::numeric_limits<float>::max)();

        set_weight(j, wmin);
        for (std::size_t i = 0; i < 100; ++ i)
          {
            estimate_feature_effect(j, training_sets);
            if (feature_useful(j))
              {
                CGAL_CLASSTRAINING_CERR << "#";
                nb_useful ++;
                min = (std::min) (min, weight(j));
                max = (std::max) (max, weight(j));
              }
            else
              CGAL_CLASSTRAINING_CERR << "-";
            set_weight(j, factor * weight(j));
          }
        CGAL_CLASSTRAINING_CERR << std::endl;
        CGAL_CLASSTRAINING_CERR << feature->name() << " useful in "
                                << nb_useful << "% of the cases, in interval [ "
                                << min << " ; " << max << " ]" << std::endl;
        if (nb_useful < 2)
          {
            set_weight(j, 0.);
            best_weights[j] = weight(j);
            continue;
          }

        feature_train.push_back (Feature_training());
        feature_train.back().i = j;
        feature_train.back().wmin = min / factor;
        feature_train.back().wmax = max * factor;

        if (best_weights[j] == 1.)
          {
            set_weight(j, 0.5 * (feature_train.back().wmin + feature_train.back().wmax));
            best_weights[j] = weight(j);
          }
        else
          set_weight(j, best_weights[j]);
        estimate_feature_effect(j, training_sets);
      }

    std::size_t nb_trials_per_feature = 1 + (std::size_t)(nb_tests / (float)(feature_train.size()));
    CGAL_CLASSIFICATION_CERR << "Trials = " << nb_tests << ", features = " << feature_train.size()
              << ", trials per feature = " << nb_trials_per_feature << std::endl;
    for (std::size_t i = 0; i < feature_train.size(); ++ i)
      feature_train[i].factor
        = std::pow (feature_train[i].wmax / feature_train[i].wmin,
                    1. / (float)nb_trials_per_feature);
    
    
    float best_score = 0.;
    float best_confidence = 0.;
    boost::tie (best_confidence, best_score)
      = compute_worst_confidence_and_score<ConcurrencyTag> (0., 0., training_sets);
    
    CGAL_CLASSIFICATION_CERR << "TRAINING GLOBALLY: Best score evolution: " << std::endl;

    CGAL_CLASSIFICATION_CERR << 100. * best_score << "% (found at initialization)" << std::endl;

    std::sort (feature_train.begin(), feature_train.end());
    for (std::size_t i = 0; i < feature_train.size(); ++ i)
      {
        const Feature_training& tr = feature_train[i];
        std::size_t current_feature_changed = tr.i;
        Feature_handle current_feature = m_features[current_feature_changed];
        
        std::size_t nb_used = 0;
        for (std::size_t j = 0; j < m_features.size(); ++ j)
          {
            if (j == current_feature_changed)
              continue;
            
            set_weight(j, best_weights[j]);
            estimate_feature_effect(j, training_sets);
            if (feature_useful(j))
              nb_used ++;
            else
              set_weight(j, 0.);
          }
        
        set_weight(current_feature_changed, tr.wmin);
        for (std::size_t j = 0; j < nb_trials_per_feature; ++ j)
          {
            estimate_feature_effect(current_feature_changed, training_sets);

            float worst_confidence = 0., worst_score = 0.;
            boost::tie (worst_confidence, worst_score)
              = compute_worst_confidence_and_score<ConcurrencyTag> (best_confidence, best_score, training_sets);

            if (worst_score > best_score
                && worst_confidence > best_confidence)
              {
                best_score = worst_score;
                best_confidence = worst_confidence;
                CGAL_CLASSIFICATION_CERR << 100. * best_score << "% (found at iteration "
                          << (i * nb_trials_per_feature) + j << "/" << nb_tests << ", "
                          << nb_used + (feature_useful(current_feature_changed) ? 1 : 0)
                          << "/" << m_features.size() << " feature(s) used)" << std::endl;
                for (std::size_t k = 0; k < m_features.size(); ++ k)
                  best_weights[k] = weight(k);
              }
            set_weight(current_feature_changed, weight(current_feature_changed) * tr.factor);
          }
      }

    for (std::size_t i = 0; i < best_weights.size(); ++ i)
      set_weight(i, best_weights[i]);

    estimate_features_effects(training_sets);
    
    CGAL_CLASSIFICATION_CERR << std::endl << "Best score found is at least " << 100. * best_score
              << "% of correct classification" << std::endl;

    std::size_t nb_removed = 0;
    for (std::size_t i = 0; i < best_weights.size(); ++ i)
      {
        Feature_handle feature = m_features[i];
        CGAL_CLASSTRAINING_CERR << "FEATURE " << feature->name() << ": " << best_weights[i] << std::endl;
        set_weight(i, best_weights[i]);

        Effect side = effect(0, i);
        bool to_remove = true;
        for (std::size_t j = 0; j < m_labels.size(); ++ j)
          {
            Label_handle clabel = m_labels[j];
            if (effect(j,i) == FAVORING)
              CGAL_CLASSTRAINING_CERR << " * Favored for ";
            else if (effect(j,i) == PENALIZING)
              CGAL_CLASSTRAINING_CERR << " * Penalized for ";
            else
              CGAL_CLASSTRAINING_CERR << " * Neutral for ";
            if (effect(j,i) != side)
              to_remove = false;
            CGAL_CLASSTRAINING_CERR << clabel->name() << std::endl;
          }
        if (to_remove)
          {
            CGAL_CLASSTRAINING_CERR << "   -> Useless! Should be removed" << std::endl;
            ++ nb_removed;
          }
      }
    CGAL_CLASSIFICATION_CERR << nb_removed
              << " feature(s) out of " << m_features.size() << " are useless" << std::endl;

    return best_score;
  }

  /// @}

private:

  float value (std::size_t label, std::size_t feature, std::size_t index) const
  {
    if (m_effect_table[label][feature] == FAVORING)
      return favored (feature, index);
    else if (m_effect_table[label][feature] == PENALIZING)
      return penalized (feature, index);
    else
      return ignored (feature, index);
  }
  
  float normalized (std::size_t feature, std::size_t index) const
  {
    return (std::max) (0.f, (std::min) (1.f, m_features[feature]->value(index) / m_weights[feature]));
  }
  float favored (std::size_t feature, std::size_t index) const
  {
    return (1. - normalized (feature, index));
  }
  float penalized (std::size_t feature, std::size_t index) const
  {
    return normalized (feature, index);
  }
  float ignored (std::size_t, std::size_t) const
  {
    return 0.5;
  }

  void estimate_features_effects(std::vector<std::vector<std::size_t> >& training_sets)
  {
    for (std::size_t i = 0; i < m_features.size(); ++ i)
      estimate_feature_effect (i, training_sets);
  }

  void estimate_feature_effect (std::size_t feature,
                                std::vector<std::vector<std::size_t> >& training_sets)
  {
    std::vector<float> mean (m_labels.size(), 0.);
                                  
    for (std::size_t j = 0; j < m_labels.size(); ++ j)
      {
        for (std::size_t k = 0; k < training_sets[j].size(); ++ k)
          {
            float val = normalized(feature, training_sets[j][k]);
            mean[j] += val;
          }
        mean[j] /= training_sets[j].size();
      }

    std::vector<float> sd (m_labels.size(), 0.);
        
    for (std::size_t j = 0; j < m_labels.size(); ++ j)
      {
        Label_handle clabel = m_labels[j];
            
        for (std::size_t k = 0; k < training_sets[j].size(); ++ k)
          {
            float val = normalized(feature, training_sets[j][k]);
            sd[j] += (val - mean[j]) * (val - mean[j]);
          }
        sd[j] = std::sqrt (sd[j] / training_sets[j].size());
        if (mean[j] - sd[j] > 0.75)
          set_effect (j, feature, FAVORING);
        else if (mean[j] + sd[j] < 0.25)
          set_effect (j, feature, PENALIZING);
        else
          set_effect (j, feature, NEUTRAL);
      }
  }

  template <typename ConcurrencyTag>
  std::pair<float, float> compute_worst_confidence_and_score (float lower_conf, float lower_score,
                                                                std::vector<std::vector<std::size_t> >& training_sets)
  {
    float worst_confidence = (std::numeric_limits<float>::max)();
    float worst_score = (std::numeric_limits<float>::max)();
    
    for (std::size_t j = 0; j < m_labels.size(); ++ j)
      {
        float confidence = 0.;
        std::size_t nb_okay = 0;

#ifndef CGAL_LINKED_WITH_TBB
        CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                                   "Parallel_tag is enabled but TBB is unavailable.");
#else
        if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
          {
            tbb::mutex mutex;
            Compute_worst_score_and_confidence f(training_sets[j], *this, j, confidence, nb_okay, mutex);
            tbb::parallel_for(tbb::blocked_range<size_t>(0, training_sets[j].size ()), f);
          }
        else
#endif
          {
            for (std::size_t k = 0; k < training_sets[j].size(); ++ k)
              {
                std::vector<std::pair<float, std::size_t> > values;

                std::vector<float> v;
                probabilities (training_sets[j][k], v);
                
                for(std::size_t l = 0; l < m_labels.size(); ++ l)
                  values.push_back (std::make_pair (v[l], l));
                
                std::sort (values.begin(), values.end());

                if (values[0].second == j)
                  {
                    confidence += values[1].first - values[0].first;
                    ++ nb_okay;
                  }
              }
          }
        
        float score = nb_okay / (float)(training_sets[j].size());
        confidence /= (float)(training_sets[j].size() * m_features.size());

        if (confidence < worst_confidence)
          worst_confidence = confidence;
        if (score < worst_score)
          worst_score = score;
        
        if (worst_confidence < lower_conf || worst_score < lower_score)
          return std::make_pair (worst_confidence, worst_score);
      }
    return std::make_pair (worst_confidence, worst_score);
  }

  bool feature_useful (std::size_t feature)
  {
    Effect side = effect(0, feature);
    for (std::size_t k = 1; k < m_labels.size(); ++ k)
      if (effect(k, feature) != side)
        return true;
    return false;
  }
  
};

}

}

#endif //  CLASSIFICATION_SUM_OF_WEIGHTED_FEATURES_PREDICATE_H
