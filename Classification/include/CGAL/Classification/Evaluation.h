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
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_EVALUATION_H
#define CGAL_CLASSIFICATION_EVALUATION_H_

#include <CGAL/Classification/Label.h>
#include <CGAL/Classification/Label_set.h>

namespace CGAL {

namespace Classification {

class Evaluation
{
  mutable std::map<Label_handle, std::size_t> m_map_labels;

  std::vector<double> m_precision;
  std::vector<double> m_recall;
  std::vector<double> m_iou; // intersection over union
  double m_accuracy;
  double m_mean_iou;
  double m_mean_f1;

public:
  
  Evaluation (const Label_set& labels,
              const std::vector<std::size_t>& ground_truth,
              const std::vector<std::size_t>& result)
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
        std::size_t gt = ground_truth[j];
        std::size_t res = result[j];
        if (gt == std::size_t(-1) || res == std::size_t(-1))
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
    
    for (std::size_t j = 0; j < labels.size(); ++ j)
      {
        m_precision[j] = true_positives[j] / double(true_positives[j] + false_positives[j]);
        m_recall[j] = true_positives[j] / double(true_positives[j] + false_negatives[j]);
        m_iou[j] = true_positives[j] / double(true_positives[j] + false_positives[j] + false_negatives[j]);

        m_mean_iou += m_iou[j];
        m_mean_f1 += 2. * (m_precision[j] * m_recall[j])
          / (m_precision[j] + m_recall[j]);
      }

    m_mean_iou /= labels.size();
    m_mean_f1 /= labels.size();
    m_accuracy = sum_true_positives / double(total);
  }


  /*!

    \brief Returns the precision of the training for the given label.

    Precision is the number of true positives divided by the sum of
    the true positives and the false positives.

  */
  double precision (Label_handle label) const
  {
    return m_precision[m_map_labels[label]];
  }

  /*!

    \brief Returns the recall of the training for the given label.

    Recall is the number of true positives divided by the sum of
    the true positives and the false negatives.

  */
  double recall (Label_handle label) const
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
  double f1_score (Label_handle label) const
  {
    std::size_t label_idx = m_map_labels[label];
    return 2. * (m_precision[label_idx] * m_recall[label_idx])
      / (m_precision[label_idx] + m_recall[label_idx]);
  }

  /*!
    \brief Returns the intersection over union of the training for the
    given label.

    Intersection over union is the number of true positives divided by
    the sum of the true positives, of the false positives and of the
    false negatives.
  */
  double intersection_over_union (Label_handle label) const
  {
    return m_iou[m_map_labels[label]];
  }

  /*!
    \brief Returns the accuracy of the training.

    Accuracy is the total number of true positives divided by the
    total number of provided inliers.
  */
  double accuracy() const { return m_accuracy; }
  
  /*!
    \brief Returns the mean \f$F_1\f$ score of the training over all
    labels (see `f1_score()`).
  */
  double mean_f1_score() const { return m_mean_f1; }
  
  /*!
    \brief Returns the mean intersection over union of the training
    over all labels (see `intersection_over_union()`).
  */
  double mean_intersection_over_union() const { return m_mean_iou; }
  
};
  
  
} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_EVALUATION_H_
