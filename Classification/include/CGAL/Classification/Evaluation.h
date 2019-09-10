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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_EVALUATION_H
#define CGAL_CLASSIFICATION_EVALUATION_H

#include <CGAL/license/Classification.h>

#include <CGAL/Classification/Label.h>
#include <CGAL/Classification/Label_set.h>
#include <map>
#include <cmath> // for std::isnan

namespace CGAL {

namespace Classification {

/*!
  \ingroup PkgClassificationDataStructures

  \brief Class to compute several measurements to evaluate the quality
  of a classification output.
*/
class Evaluation
{
  mutable std::map<Label_handle, std::size_t> m_map_labels;

  std::vector<float> m_precision;
  std::vector<float> m_recall;
  std::vector<float> m_iou; // intersection over union
  float m_accuracy;
  float m_mean_iou;
  float m_mean_f1;

public:

  /// \name Constructor
  /// @{

  
/*!

  \brief Instantiates an evaluation object and computes all
  measurements.

  \param labels labels used.

  \param ground_truth vector of label indices: it should contain the
  index of the corresponding label in the `Label_set` provided in the
  constructor. Input items that do not have a ground truth information
  should be given the value `-1`.

  \param result similar to `ground_truth` but contained the result of
  a classification.

*/
  template <typename LabelIndexRange>
  Evaluation (const Label_set& labels,
              const LabelIndexRange& ground_truth,
              const LabelIndexRange& result)
    : m_precision (labels.size()),
      m_recall (labels.size()),
      m_iou (labels.size())
  {
    for (std::size_t i = 0; i < labels.size(); ++ i)
      m_map_labels[labels[i]] = i;

    std::vector<std::size_t> true_positives (labels.size());
    std::vector<std::size_t> false_positives (labels.size());
    std::vector<std::size_t> false_negatives (labels.size());

    std::size_t sum_true_positives = 0;
    std::size_t total = 0;
    
    for (std::size_t j = 0; j < ground_truth.size(); ++ j)
    {
      int gt = static_cast<int>(ground_truth[j]);
      int res = static_cast<int>(result[j]);
      if (gt == -1 || res == -1)
        continue;
      ++ total;
      if (gt == res)
      {
        ++ true_positives[gt];
        ++ sum_true_positives;
        continue;
      }
      ++ false_positives[res];
      ++ false_negatives[gt];
    }

    m_mean_iou = 0.;
    m_mean_f1 = 0.;

    std::size_t correct_labels = 0;
    
    for (std::size_t j = 0; j < labels.size(); ++ j)
    {
      m_precision[j] = true_positives[j] / float(true_positives[j] + false_positives[j]);
      m_recall[j] = true_positives[j] / float(true_positives[j] + false_negatives[j]);
      m_iou[j] = true_positives[j] / float(true_positives[j] + false_positives[j] + false_negatives[j]);

      if (std::isnan(m_iou[j]))
        continue;

      ++ correct_labels;
      m_mean_iou += m_iou[j];
      m_mean_f1 += 2.f * (m_precision[j] * m_recall[j])
        / (m_precision[j] + m_recall[j]);
    }

    m_mean_iou /= correct_labels;
    m_mean_f1 /= correct_labels;
    m_accuracy = sum_true_positives / float(total);
  }

  /// @}

  /// \name Label Evaluation
  /// @{


  /*!

    \brief Returns the precision of the training for the given label.

    Precision is the number of true positives divided by the sum of
    the true positives and the false positives.

  */
  float precision (Label_handle label) const
  {
    return m_precision[m_map_labels[label]];
  }

  /*!

    \brief Returns the recall of the training for the given label.

    Recall is the number of true positives divided by the sum of
    the true positives and the false negatives.

  */
  float recall (Label_handle label) const
  {
    return m_recall[m_map_labels[label]];
  }

  /*!

    \brief Returns the \f$F_1\f$ score of the training for the given label.

    \f$F_1\f$ score is the harmonic mean of `precision()` and `recall()`:

    \f[
    F_1 = 2 \times \frac{precision \times recall}{precision + recall}
    \f]

  */
  float f1_score (Label_handle label) const
  {
    std::size_t label_idx = m_map_labels[label];
    return 2.f * (m_precision[label_idx] * m_recall[label_idx])
      / (m_precision[label_idx] + m_recall[label_idx]);
  }

  /*!
    \brief Returns the intersection over union of the training for the
    given label.

    Intersection over union is the number of true positives divided by
    the sum of the true positives, of the false positives and of the
    false negatives.
  */
  float intersection_over_union (Label_handle label) const
  {
    return m_iou[m_map_labels[label]];
  }

  /// @}
  
  /// \name Global Evaluation
  /// @{

  
  /*!
    \brief Returns the accuracy of the training.

    Accuracy is the total number of true positives divided by the
    total number of provided inliers.
  */
  float accuracy() const { return m_accuracy; }
  
  /*!
    \brief Returns the mean \f$F_1\f$ score of the training over all
    labels (see `f1_score()`).
  */
  float mean_f1_score() const { return m_mean_f1; }
  
  /*!
    \brief Returns the mean intersection over union of the training
    over all labels (see `intersection_over_union()`).
  */
  float mean_intersection_over_union() const { return m_mean_iou; }

  /// @}
  
};
  
  
} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_EVALUATION_H

