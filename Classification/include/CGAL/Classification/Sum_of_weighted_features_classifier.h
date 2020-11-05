// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot, Florent Lafarge

#ifndef CLASSIFICATION_SUM_OF_WEIGHTED_FEATURES_CLASSIFIER_H
#define CLASSIFICATION_SUM_OF_WEIGHTED_FEATURES_CLASSIFIER_H

#include <CGAL/license/Classification.h>

#include <CGAL/Classification/Feature_set.h>
#include <CGAL/Classification/Label_set.h>
#include <CGAL/Classification/internal/verbosity.h>
#include <CGAL/tags.h>
#include <CGAL/algorithm.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <map>
#include <iostream>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>
#include <mutex>
#endif // CGAL_LINKED_WITH_TBB

#define CLASSIFICATION_TRAINING_QUICK_ESTIMATION

namespace CGAL {

namespace Classification {

/*!
  \ingroup PkgClassificationClassifiers

  \brief %Classifier based on the sum of weighted features with
  user-defined effects on labels.

  \cgalModels `CGAL::Classification::Classifier`
*/
class Sum_of_weighted_features_classifier
{
public:

  enum Effect /// Defines the effect of a feature on a type.
  {
    FAVORING = 0, ///< High values of the feature favor this type
    NEUTRAL = 1, ///< The feature has no effect on this type
    PENALIZING = 2 ///< Low values of the feature favor this type
  };

private:

#ifdef CGAL_LINKED_WITH_TBB
  class Compute_iou
  {
    std::vector<std::size_t>& m_training_set;
    const Sum_of_weighted_features_classifier& m_classifier;
    std::size_t m_label;
    std::vector<std::size_t>& m_true_positives;
    std::vector<std::size_t>& m_false_positives;
    std::vector<std::size_t>& m_false_negatives;
    std::vector<std::mutex>& m_tp_mutex;
    std::vector<std::mutex>& m_fp_mutex;
    std::vector<std::mutex>& m_fn_mutex;


  public:

    Compute_iou (std::vector<std::size_t>& training_set,
                 const Sum_of_weighted_features_classifier& classifier,
                 std::size_t label,
                 std::vector<std::size_t>& true_positives,
                 std::vector<std::size_t>& false_positives,
                 std::vector<std::size_t>& false_negatives,
                 std::vector<std::mutex>& tp_mutex,
                 std::vector<std::mutex>& fp_mutex,
                 std::vector<std::mutex>& fn_mutex)
      : m_training_set (training_set)
      , m_classifier (classifier)
      , m_label (label)
      , m_true_positives (true_positives)
      , m_false_positives (false_positives)
      , m_false_negatives (false_negatives)
      , m_tp_mutex (tp_mutex)
      , m_fp_mutex (fp_mutex)
      , m_fn_mutex (fn_mutex)
    { }

    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t k = r.begin(); k != r.end(); ++ k)
      {
        std::size_t res = 0;

        std::vector<float> v;
        m_classifier (m_training_set[k], v);

        float max = 0.f;
        for(std::size_t l = 0; l < v.size(); ++ l)
          if (v[l] > max)
          {
            max = v[l];
            res = l;
          }

        if (m_label == res)
        {
          m_tp_mutex[m_label].lock();
          ++ m_true_positives[m_label];
          m_tp_mutex[m_label].unlock();
          continue;
        }
        m_fp_mutex[res].lock();
        ++ m_false_positives[res];
        m_fp_mutex[res].unlock();

        m_fn_mutex[m_label].lock();
        ++ m_false_negatives[m_label];
        m_fn_mutex[m_label].unlock();
      }
    }

  };
#endif // CGAL_LINKED_WITH_TBB

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

  const Label_set& m_labels;
  const Feature_set& m_features;
  std::vector<float> m_weights;
  std::vector<std::vector<Effect> > m_effect_table;
  mutable std::map<Label_handle, std::size_t> m_map_labels;
  mutable std::map<Feature_handle, std::size_t> m_map_features;

public:

  /// \name Constructor
  /// @{

/*!

  \brief instantiates the classifier using the sets of `labels` and `features`.

  \note If the label set of the feature set are modified after
  instantiating this object (addition of removal of a label and/or of
  a feature), another classifier object should be instantiated as the
  internal data structures of this one are invalidated.
*/
  Sum_of_weighted_features_classifier (const Label_set& labels,
                                       const Feature_set& features)
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

  /// @}

  /// \name Weights and Effects
  /// @{

  /*!
    \brief sets the weight of `feature` (`weight` must be positive).
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
    \brief returns the weight of `feature`.
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
    \brief sets the `effect` of `feature` on `label`.
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
    \brief returns the `effect` of `feature` on `label`.
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

  /// @}

  /// \cond SKIP_IN_MANUAL
  void operator() (std::size_t item_index,
                   std::vector<float>& out) const
  {
    out.resize (m_labels.size());
    for (std::size_t l = 0; l < m_labels.size(); ++ l)
    {
      out[l] = 0.;
      for (std::size_t f = 0; f < m_features.size(); ++ f)
        if (weight(f) != 0.)
          out[l] += value (l, f, item_index);
      out[l] = std::exp (-out[l]);
    }
  }
  /// \endcond

  /// \name Training
  /// @{

  /*!
    \brief runs the training algorithm.

    From the set of provided ground truth, this algorithm estimates
    the sets of weights and effects that produce the most accurate
    result with respect to this ground truth. Old weights and effects
    are discarded.

    \pre At least one ground truth item should be assigned to each
    label.

    \param ground_truth vector of label indices. It should contain for
    each input item, in the same order as the input set, the index of
    the corresponding label in the `Label_set` provided in the
    constructor. Input items that do not have a ground truth
    information should be given the value `-1`.

    \param nb_tests number of tests to perform. Higher values may
    provide the user with better results at the cost of a higher
    computation time. Using a value of at least 10 times the number of
    features is advised.

    \return mean intersection-over-union over each label between the
    provided ground truth and the best classification found by the
    training set.
  */
  template <typename ConcurrencyTag, typename LabelIndexRange>
  float train (const LabelIndexRange& ground_truth,
               unsigned int nb_tests = 300)
  {
    CGAL_precondition (m_labels.is_valid_ground_truth (ground_truth));

    std::vector<std::vector<std::size_t> > training_sets (m_labels.size());
    std::size_t nb_tot = 0;
    std::size_t i = 0;
    for (const auto& gt_value : ground_truth)
    {
      if (int(gt_value) != -1)
      {
        training_sets[std::size_t(gt_value)].push_back (i);
        ++ nb_tot;
      }
      ++ i;
    }

#ifdef CLASSIFICATION_TRAINING_QUICK_ESTIMATION
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      CGAL::cpp98::random_shuffle (training_sets[i].begin(), training_sets[i].end());
#endif

    CGAL_CLASSIFICATION_CERR << "Training using " << nb_tot << " inliers" << std::endl;

    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (training_sets.size() <= i || training_sets[i].empty())
        std::cerr << "WARNING: \"" << m_labels[i]->name() << "\" doesn't have a training set." << std::endl;

    std::vector<float> best_weights (m_features.size(), 1.);

    std::vector<Feature_training> feature_train;
    std::size_t nb_trials = 100;
    float wmin = 1e-5f, wmax = 1e5f;
    float factor = std::pow (wmax/wmin, 1.f / float(nb_trials));

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
        set_weight(j, 0.5f * (feature_train.back().wmin + feature_train.back().wmax));
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
                    1.f / float(nb_trials_per_feature));


    float best_score = 0.;
    best_score = compute_mean_iou<ConcurrencyTag>(training_sets);

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

        float worst_score = 0.;
        worst_score = compute_mean_iou<ConcurrencyTag>(training_sets);
        if (worst_score > best_score)
        {
          best_score = worst_score;
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

  /// \cond SKIP_IN_MANUAL
  template <typename ConcurrencyTag, typename LabelIndexRange>
  float train_random (const LabelIndexRange& ground_truth,
                      unsigned int nb_tests = 300)
  {
    std::vector<std::vector<std::size_t> > training_sets (m_labels.size());
    std::size_t nb_tot = 0;
    for (std::size_t i = 0; i < ground_truth.size(); ++ i)
      if (ground_truth[i] != -1)
      {
        training_sets[std::size_t(ground_truth[i])].push_back (i);
        ++ nb_tot;
      }

    CGAL_CLASSIFICATION_CERR << "Training using " << nb_tot << " inliers" << std::endl;

    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (training_sets.size() <= i || training_sets[i].empty())
        std::cerr << "WARNING: \"" << m_labels[i]->name() << "\" doesn't have a training set." << std::endl;

    std::vector<float> best_weights (m_features.size(), 1.);

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

    CGAL_CLASSIFICATION_CERR << "Trials = " << nb_tests << ", features = " << feature_train.size() << std::endl;


    float best_score = compute_mean_iou<ConcurrencyTag>(training_sets);

    CGAL_CLASSIFICATION_CERR << "TRAINING GLOBALLY: Best score evolution: " << std::endl;

    CGAL_CLASSIFICATION_CERR << 100. * best_score << "% (found at initialization)" << std::endl;

    for (std::size_t i = 0; i < std::size_t(nb_tests); ++ i)
    {
      std::size_t nb_used = 0;
      std::size_t j = rand() % feature_train.size();
      set_weight (feature_train[j].i,
                  feature_train[j].wmin + ((feature_train[j].wmax - feature_train[j].wmin)
                                           * (rand() / float(RAND_MAX))));
      estimate_feature_effect(feature_train[j].i, training_sets);

      float worst_score = compute_mean_iou<ConcurrencyTag>(training_sets);

      if (worst_score > best_score)
      {
        best_score = worst_score;
        CGAL_CLASSIFICATION_CERR << 100. * best_score << "% (found at iteration "
                                 << i << "/" << nb_tests << ", "
                                 << nb_used
                                 << "/" << m_features.size() << " feature(s) used)" << std::endl;
        for (std::size_t k = 0; k < m_features.size(); ++ k)
          best_weights[k] = weight(k);
      }
      set_weight (feature_train[j].i,
                  best_weights[feature_train[j].i]);
      estimate_feature_effect(feature_train[j].i, training_sets);
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
  /// \endcond

  /// \name Input/Output
  /// @{

  /*!
    \brief saves the current configuration in the stream `output`.

    This allows to easily save and recover a specific classification
    configuration, that is to say:

    - The weight of each feature
    - The effects of each feature on each label

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
    boost::property_tree::write_xml(output, tree,
#if BOOST_VERSION >= 105600
                                    boost::property_tree::xml_writer_make_settings<std::string>(' ', 3));
#else
                                    boost::property_tree::xml_writer_make_settings<char>(' ', 3));
#endif
  }

  /*!
    \brief loads a configuration from the stream `input`. A
    configuration is a set of weights and effects.

    The input file should be in the XML format written by the
    `save_configuration()` method. Labels and features are described
    in the XML file by their name and the corresponding `Label` and
    `Feature_base` object should therefore be given the same names as
    the ones they had when saving the configuration.

    \note If a feature (or label) found in the input file is not found
    in the `Feature_set` (`Label_set`) provided by the user in the
    constructor, and if `verbose` is set up to `true`, a warning is
    displayed.

    \note If a feature (or label) provided by the user in the
    constructor is not described in the input file, the corresponding
    weights and effects are kept to their default values (1 for the
    weight and `NEUTRAL` for the effect).

    \param input input stream.
    \param verbose displays warning if set to `true`. The method is
    silent otherwise.

    \return `true` if all weights and effects found in the
    configuration file were applicable to the feature set and label
    set of this classifier, `false` otherwise.
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

    for(boost::property_tree::ptree::value_type& v : tree.get_child("classification.features"))
    {
      std::string name = v.second.get<std::string>("name");
      std::map<std::string, std::size_t>::iterator
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

    for(boost::property_tree::ptree::value_type& v : tree.get_child("classification.labels"))
    {
      std::string label_name = v.second.get<std::string>("name");
      std::map<std::string, std::size_t>::iterator
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

      for(boost::property_tree::ptree::value_type& v2 : v.second)
      {
        if (v2.first == "name")
          continue;

        std::string feature_name = v2.second.get<std::string>("name");

        std::map<std::string, std::size_t>::iterator
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
    return (1.f - normalized (feature, index));
  }
  float penalized (std::size_t feature, std::size_t index) const
  {
    return normalized (feature, index);
  }
  float ignored (std::size_t, std::size_t) const
  {
    return 0.5f;
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
#ifdef CLASSIFICATION_TRAINING_QUICK_ESTIMATION
      std::size_t training_set_size = (std::min) (std::size_t(0.1 * training_sets[j].size()),
                                                  std::size_t(10000));
#else
      std::size_t training_set_size = training_sets[j].size();
#endif

      for (std::size_t k = 0; k < training_set_size; ++ k)
      {
        float val = normalized(feature, training_sets[j][k]);
        mean[j] += val;
      }
      mean[j] /= training_set_size;
    }

    std::vector<float> sd (m_labels.size(), 0.);

    for (std::size_t j = 0; j < m_labels.size(); ++ j)
    {
      Label_handle clabel = m_labels[j];

#ifdef CLASSIFICATION_TRAINING_QUICK_ESTIMATION
      std::size_t training_set_size = (std::min) (std::size_t(0.1 * training_sets[j].size()),
                                                  std::size_t(10000));
#else
      std::size_t training_set_size = training_sets[j].size();
#endif
      for (std::size_t k = 0; k < training_set_size; ++ k)
      {
        float val = normalized(feature, training_sets[j][k]);
        sd[j] += (val - mean[j]) * (val - mean[j]);
      }
      sd[j] = std::sqrt (sd[j] / training_set_size);
      if (mean[j] - sd[j] > (2./3.))
        set_effect (j, feature, FAVORING);
      else if (mean[j] + sd[j] < (1./3.))
        set_effect (j, feature, PENALIZING);
      else
        set_effect (j, feature, NEUTRAL);
    }
  }

  template <typename ConcurrencyTag>
  float compute_mean_iou (std::vector<std::vector<std::size_t> >& training_sets)
  {
    std::vector<std::size_t> true_positives (m_labels.size());
    std::vector<std::size_t> false_positives (m_labels.size());
    std::vector<std::size_t> false_negatives (m_labels.size());

    for (std::size_t j = 0; j < training_sets.size(); ++ j)
    {
      std::size_t gt = j;

#ifndef CGAL_LINKED_WITH_TBB
      CGAL_static_assertion_msg (!(std::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                                 "Parallel_tag is enabled but TBB is unavailable.");
#else
      if (std::is_convertible<ConcurrencyTag,Parallel_tag>::value)
      {
        std::vector<std::mutex> tp_mutex (m_labels.size());
        std::vector<std::mutex> fp_mutex (m_labels.size());
        std::vector<std::mutex> fn_mutex (m_labels.size());
        Compute_iou f(training_sets[j], *this, j,
                      true_positives, false_positives, false_negatives,
                      tp_mutex, fp_mutex, fn_mutex);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, training_sets[j].size ()), f);
      }
      else
#endif
        for (std::size_t k = 0; k < training_sets[j].size(); ++ k)
        {
          std::size_t res = 0;

          std::vector<float> v;
          (*this) (training_sets[j][k], v);

          float max = 0.f;
          for(std::size_t l = 0; l < m_labels.size(); ++ l)
            if (v[l] > max)
            {
              max = v[l];
              res = l;
            }

          if (gt == res)
          {
            ++ true_positives[gt];
            continue;
          }
          ++ false_positives[res];
          ++ false_negatives[gt];
        }
    }

    float out = 0.;

    for (std::size_t j = 0; j < m_labels.size(); ++ j)
    {
      float iou = true_positives[j] / float(true_positives[j] + false_positives[j] + false_negatives[j]);
      out += iou;
    }

    return out / m_labels.size();
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

#endif //  CLASSIFICATION_SUM_OF_WEIGHTED_FEATURES_CLASSIFIER_H
