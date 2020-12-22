// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_EVALUATION_H
#define CGAL_CLASSIFICATION_EVALUATION_H

#include <CGAL/license/Classification.h>

#include <CGAL/Classification/Label.h>
#include <CGAL/Classification/Label_set.h>

#include <CGAL/Iterator_range.h>

#include <boost/iterator/zip_iterator.hpp>

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
  const Label_set& m_labels;
  std::vector<std::vector<std::size_t> > m_confusion; // confusion matrix

public:

  /// \name Constructors
  /// @{

  /*!
    \brief instantiates an empty evaluation object.

    \param labels labels used.
  */
  Evaluation (const Label_set& labels)
    : m_labels (labels)
  {
    init();
  }

  /*!

    \brief instantiates an evaluation object and computes all
    measurements.

    \param labels labels used.

    \param ground_truth vector of label indices: it should contain the
    index of the corresponding label in the `Label_set` provided in the
    constructor. Input items that do not have a ground truth information
    should be given the value `-1`.

    \param result similar to `ground_truth` but contained the result of
    a classification.

  */
  template <typename GroundTruthIndexRange, typename ResultIndexRange>
  Evaluation (const Label_set& labels,
              const GroundTruthIndexRange& ground_truth,
              const ResultIndexRange& result)
    : m_labels (labels)
  {
    init();
    append(ground_truth, result);
  }

  /// \cond SKIP_IN_MANUAL
  void init()
  {
    m_confusion.resize (m_labels.size());
    for (std::size_t i = 0; i < m_confusion.size(); ++ i)
      m_confusion[i].resize (m_labels.size(), 0);
  }

  bool label_has_ground_truth (std::size_t label_idx) const
  {
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (m_confusion[i][label_idx] != 0)
        return true;
    return false;
  }
  /// \endcond

  /// @}

  /// \name Modification
  /// @{

  /*!
    \brief appends more items to the evaluation object.

    \param ground_truth vector of label indices: it should contain the
    index of the corresponding label in the `Label_set` provided in the
    constructor. Input items that do not have a ground truth information
    should be given the value `-1`.

    \param result similar to `ground_truth` but contained the result of
    a classification.

  */
  template <typename GroundTruthIndexRange, typename ResultIndexRange>
  void append (const GroundTruthIndexRange& ground_truth,
               const ResultIndexRange& result)
  {
    CGAL_precondition (m_labels.is_valid_ground_truth (ground_truth));
    CGAL_precondition (m_labels.is_valid_ground_truth (result));

    for (const auto& p : CGAL::make_range
           (boost::make_zip_iterator(boost::make_tuple(ground_truth.begin(), result.begin())),
            boost::make_zip_iterator(boost::make_tuple(ground_truth.end(), result.end()))))
    {
      int gt = static_cast<int>(boost::get<0>(p));
      int res = static_cast<int>(boost::get<1>(p));
      if (gt == -1 || res == -1)
        continue;

      ++ m_confusion[std::size_t(res)][std::size_t(gt)];
    }
  }

  /// @}

  /// \name Label Evaluation
  /// @{

  /*!
    \brief returns the number of items whose ground truth is
    `ground_truth` and which were classified as `result`.
  */
  std::size_t confusion (Label_handle ground_truth, Label_handle result)
  {
    std::size_t idx_gt = ground_truth->index();
    std::size_t idx_r = result->index();
    return m_confusion[idx_gt][idx_r];
  }

  /*!

    \brief returns the precision of the training for the given label.

    Precision is the number of true positives divided by the sum of
    the true positives and the false positives.

  */
  float precision (Label_handle label) const
  {
    std::size_t idx = label->index();
    if (!label_has_ground_truth(idx))
      return std::numeric_limits<float>::quiet_NaN();

    std::size_t total = 0;
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      total += m_confusion[idx][i];

    if (total == 0)
      return 0.f;

    return m_confusion[idx][idx] / float(total);
  }

  /*!

    \brief returns the recall of the training for the given label.

    Recall is the number of true positives divided by the sum of
    the true positives and the false negatives.

  */
  float recall (Label_handle label) const
  {
    std::size_t idx = label->index();
    if (!label_has_ground_truth(idx))
      return std::numeric_limits<float>::quiet_NaN();

    std::size_t total = 0;
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      total += m_confusion[i][idx];
    return m_confusion[idx][idx] / float(total);
  }

  /*!

    \brief returns the \f$F_1\f$ score of the training for the given label.

    \f$F_1\f$ score is the harmonic mean of `precision()` and `recall()`:

    \f[
    F_1 = 2 \times \frac{precision \times recall}{precision + recall}
    \f]

  */
  float f1_score (Label_handle label) const
  {
    float p = precision(label);
    float r = recall(label);

    if (p == 0.f && r == 0.f)
      return 0.f;

    return 2.f * p * r / (p + r);
  }

  /*!
    \brief returns the intersection over union of the training for the
    given label.

    Intersection over union is the number of true positives divided by
    the sum of the true positives, of the false positives and of the
    false negatives.
  */
  float intersection_over_union (Label_handle label) const
  {
    std::size_t idx = label->index();

    std::size_t total = 0;
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
    {
      total += m_confusion[i][idx];
      if (i != idx)
        total += m_confusion[idx][i];
    }

    return m_confusion[idx][idx] / float(total);
  }

  /// @}

  /// \name Global Evaluation
  /// @{

  /*!
    \brief returns the number of misclassified items.
  */
  std::size_t number_of_misclassified_items() const
  {
    std::size_t total = 0;
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      for (std::size_t j = 0; j < m_labels.size(); ++ j)
        if (i != j)
          total += m_confusion[i][j];
    return total;
  }

  /*!
    \brief returns the total number of items used for evaluation.
  */
  std::size_t number_of_items() const
  {
    std::size_t total = 0;
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      for (std::size_t j = 0; j < m_labels.size(); ++ j)
        total += m_confusion[i][j];
    return total;
  }

  /*!
    \brief returns the accuracy of the training.

    Accuracy is the total number of true positives divided by the
    total number of provided inliers.
  */
  float accuracy() const
  {
    std::size_t true_positives = 0;
    std::size_t total = 0;
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
    {
      true_positives += m_confusion[i][i];
      for (std::size_t j = 0; j < m_labels.size(); ++ j)
        total += m_confusion[i][j];
    }
    return true_positives / float(total);
  }

  /*!
    \brief returns the mean \f$F_1\f$ score of the training over all
    labels (see `f1_score()`).
  */
  float mean_f1_score() const
  {
    float mean = 0;
    std::size_t nb = 0;
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (label_has_ground_truth(i))
      {
        mean += f1_score(m_labels[i]);
        ++ nb;
      }
    return mean / nb;
  }

  /*!
    \brief returns the mean intersection over union of the training
    over all labels (see `intersection_over_union()`).
  */
  float mean_intersection_over_union() const
  {
    float mean = 0;
    std::size_t nb = 0;
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
    {
      float iou = intersection_over_union(m_labels[i]);
      if (!std::isnan(iou))
      {
        mean += iou;
        ++ nb;
      }
    }
    return mean / nb;
  }

  /// @}

  /// \name Output Formatting Functions
  /// @{

  /*!
    \brief outputs the evaluation in a simple ASCII format to the stream `os`.
  */
  friend std::ostream& operator<< (std::ostream& os, const Evaluation& evaluation)
  {
    os << "Evaluation of classification:" << std::endl;
    os << " * Global results:" << std::endl;
    os << "   - " << evaluation.number_of_misclassified_items()
       << " misclassified item(s) out of " << evaluation.number_of_items() << std::endl
       << "   - Accuracy = " << evaluation.accuracy() << std::endl
       << "   - Mean F1 score = " << evaluation.mean_f1_score() << std::endl
       << "   - Mean IoU = " << evaluation.mean_intersection_over_union() << std::endl;
    os << " * Detailed results:" << std::endl;
    for (std::size_t i = 0; i < evaluation.m_labels.size(); ++ i)
    {
      os << "   - \"" << evaluation.m_labels[i]->name() << "\": ";
      if (evaluation.label_has_ground_truth(i))
        os << "Precision = " << evaluation.precision(evaluation.m_labels[i]) << " ; "
           << "Recall = " << evaluation.recall(evaluation.m_labels[i]) << " ; "
           << "F1 score = " << evaluation.f1_score(evaluation.m_labels[i]) << " ; "
           << "IoU = " << evaluation.intersection_over_union(evaluation.m_labels[i]) << std::endl;
      else
        os << "(no ground truth)" << std::endl;
    }
    return os;
  }

  /*!
    \brief outputs the evaluation as an HTML page to the stream `os`.
  */
  static std::ostream& output_to_html (std::ostream& os, const Evaluation& evaluation)
  {
    os <<  "<!DOCTYPE html>" << std::endl
       << "<html>" << std::endl
       << "<head>" << std::endl
       << "<style type=\"text/css\">" << std::endl
       << "  body{margin:40px auto; max-width:900px; line-height:1.5; color:#333}" << std::endl
       << "  h1,h2{line-height:1.2}" << std::endl
       << "  table{width:100%}" << std::endl
       << "  table,th,td{border: 1px solid black; border-collapse: collapse; }" << std::endl
       << "  th,td{padding: 5px;}" << std::endl
       << "</style>" << std::endl
       << "<title>Evaluation of CGAL Classification results</title>" << std::endl
       << "</head>" << std::endl
       << "<body>" << std::endl
       << "<h1>Evaluation of CGAL Classification results</h1>" << std::endl;

    os <<  "<h2>Global Results</h2>" << std::endl
       << "<ul>" << std::endl
       << "  <li>" << evaluation.number_of_misclassified_items()
       << " misclassified item(s) out of " << evaluation.number_of_items() << "</li>" << std::endl
       << "  <li>Accuracy = " << evaluation.accuracy() << "</li>" << std::endl
       << "  <li>Mean F1 score = " << evaluation.mean_f1_score() << "</li>" << std::endl
       << "  <li>Mean IoU = " << evaluation.mean_intersection_over_union() << "</li>" << std::endl
       << "</ul>" << std::endl;

    const Label_set& labels = evaluation.m_labels;

    os <<  "<h2>Detailed Results</h2>" << std::endl
       << "<table>" << std::endl
       << "  <tr>" << std::endl
       << "    <th><strong>Label</strong></th>" << std::endl
       << "    <th><strong>Precision</strong></th>" << std::endl
       << "    <th><strong>Recall</strong></th>" << std::endl
       << "    <th><strong>F1 score</strong></th>" << std::endl
       << "    <th><strong>IoU</strong></th>" << std::endl
       << "  </tr>" << std::endl;
    for (std::size_t i = 0; i < labels.size(); ++ i)
      if (evaluation.label_has_ground_truth(i))
        os <<  "  <tr>" << std::endl
           << "    <td>" << labels[i]->name() << "</td>" << std::endl
           << "    <td>" << evaluation.precision(labels[i]) << "</td>" << std::endl
           << "    <td>" << evaluation.recall(labels[i]) << "</td>" << std::endl
           << "    <td>" << evaluation.f1_score(labels[i]) << "</td>" << std::endl
           << "    <td>" << evaluation.intersection_over_union(labels[i]) << "</td>" << std::endl
           << "  </tr>" << std::endl;
      else
        os <<  "  <tr>" << std::endl
           << "    <td>" << labels[i]->name() << "</td>" << std::endl
           << "    <td><em>(no ground truth)</em></td>" << std::endl
           << "    <td></td>" << std::endl
           << "    <td></td>" << std::endl
           << "    <td></td>" << std::endl
           << "  </tr>" << std::endl;

    os <<  "</table>" << std::endl;

    os <<  "<h2>Confusion Matrix</h2>" << std::endl
       << "<table>" << std::endl
       << "  <tr>" << std::endl
       << "    <th></th>" << std::endl;
    for (std::size_t i = 0; i < labels.size(); ++ i)
      os <<  "    <th><strong>" << labels[i]->name() << "</strong></th>" << std::endl;
    os << "    <th><strong>PREDICTIONS</strong></th>" << std::endl;
    os <<  "  </tr>" << std::endl;

    std::vector<std::size_t> sums (labels.size(), 0);
    for (std::size_t i = 0; i < labels.size(); ++ i)
    {
      os <<  "  <tr>" << std::endl
         << "    <td><strong>" << labels[i]->name() << "</strong></td>" << std::endl;
      std::size_t sum = 0;
      for (std::size_t j = 0; j < labels.size(); ++ j)
      {
        if (i == j)
          os <<  "    <td><strong>" << evaluation.m_confusion[i][j] << "</strong></td>" << std::endl;
        else
          os <<  "    <td>" << evaluation.m_confusion[i][j] << "</td>" << std::endl;
        sum += evaluation.m_confusion[i][j];
        sums[j] += evaluation.m_confusion[i][j];
      }
      os << "    <td><strong>" << sum << "</strong></td>" << std::endl;
      os <<  "  </tr>" << std::endl;
    }

    os <<  "  <tr>" << std::endl
       << "    <td><strong>GROUND TRUTH</strong></td>" << std::endl;
    std::size_t total = 0;
    for (std::size_t j = 0; j < labels.size(); ++ j)
    {
      os << "    <td><strong>" << sums[j] << "</strong></td>" << std::endl;
      total += sums[j];
    }
    os << "    <td><strong>" << total << "</strong></td>" << std::endl
       <<  "  </tr>" << std::endl
       <<  "</table>" << std::endl
       << "<p><em>This page was generated by the <a \
href=\"https://doc.cgal.org/latest/Classification/index.html\">CGAL \
Classification package</a>.</em></p>" << std::endl
       <<  "</body>" << std::endl
       << "</html>" << std::endl;

    return os;
  }

  /// @}

};


} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_EVALUATION_H
