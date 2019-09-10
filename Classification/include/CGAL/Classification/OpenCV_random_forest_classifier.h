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

#ifndef CGAL_CLASSIFICATION_OPENCV_RANDOM_FOREST_CLASSIFIER_H
#define CGAL_CLASSIFICATION_OPENCV_RANDOM_FOREST_CLASSIFIER_H

#include <CGAL/license/Classification.h>

#include <CGAL/Classification/Feature_set.h>
#include <CGAL/Classification/Label_set.h>
#if (CV_MAJOR_VERSION < 3)
#include <cv.h>
#include <ml.h>
#else
#include <opencv/cv.h>
#include <opencv/ml.h>
#endif

namespace CGAL {

namespace Classification {

/*!
  \ingroup PkgClassificationClassifiers

  \brief %Classifier based on the OpenCV version of random forest algorithm.

  \note This class requires the \ref thirdpartyOpenCV library.

  \cgalModels `CGAL::Classification::Classifier`
*/
class OpenCV_random_forest_classifier
{
  const Label_set& m_labels;
  const Feature_set& m_features;
  int m_max_depth;
  int m_min_sample_count;
  int m_max_categories;
  int m_max_number_of_trees_in_the_forest;
  float m_forest_accuracy;
  
#if (CV_MAJOR_VERSION < 3)
  CvRTrees* rtree;
#else
  cv::Ptr<cv::ml::RTrees> rtree;
#endif

public:
  
  /// \name Constructor
  /// @{
  
/*!
  \brief Instantiate the classifier using the sets of `labels` and `features`.

  Parameters documentation is copy-pasted from [the official documentation of OpenCV](http://docs.opencv.org/2.4/modules/ml/doc/random_trees.html). For more details on this method, please refer to it.

  \param labels label set used.
  \param features feature set used.
  \param max_depth the depth of the tree. A low value will likely underfit and conversely a high value will likely overfit. The optimal value can be obtained using cross validation or other suitable methods.
  \param min_sample_count minimum samples required at a leaf node for it to be split. A reasonable value is a small percentage of the total data e.g. 1%.
  \param max_categories Cluster possible values of a categorical variable into \f$ K \leq max\_categories \f$ clusters to find a suboptimal split. If a discrete variable, on which the training procedure tries to make a split, takes more than max_categories values, the precise best subset estimation may take a very long time because the algorithm is exponential. Instead, many decision trees engines (including ML) try to find sub-optimal split in this case by clustering all the samples into max_categories clusters that is some categories are merged together. The clustering is applied only in \f$ n>2-class \f$ classification problems for categorical variables with \f$ N > max\_categories \f$ possible values. In case of regression and 2-class classification the optimal split can be found efficiently without employing clustering, thus the parameter is not used in these cases.
  \param max_number_of_trees_in_the_forest The maximum number of trees in the forest (surprise, surprise). Typically the more trees you have the better the accuracy. However, the improvement in accuracy generally diminishes and asymptotes pass a certain number of trees. Also to keep in mind, the number of tree increases the prediction time linearly.
  \param forest_accuracy Sufficient accuracy (OOB error).
*/
  OpenCV_random_forest_classifier (const Label_set& labels,
                                   const Feature_set& features,
                                   int max_depth = 20,
                                   int min_sample_count = 5,
                                   int max_categories = 15,
                                   int max_number_of_trees_in_the_forest = 100,
                                   float forest_accuracy = 0.01f)
    : m_labels (labels), m_features (features),
      m_max_depth (max_depth), m_min_sample_count (min_sample_count),
      m_max_categories (max_categories),
      m_max_number_of_trees_in_the_forest (max_number_of_trees_in_the_forest),
      m_forest_accuracy (forest_accuracy)
#if (CV_MAJOR_VERSION < 3)
    , rtree (NULL)
#endif
  {  }

  /// \cond SKIP_IN_MANUAL
  ~OpenCV_random_forest_classifier ()
  {
#if (CV_MAJOR_VERSION < 3)
    if (rtree != NULL)
      delete rtree;
#endif
  }
  /// \endcond
  
  /// @}

  /// \name Parameters
  /// @{
  
  void set_max_depth (int max_depth) { m_max_depth = max_depth; }
  void set_min_sample_count (int min_sample_count) { m_min_sample_count = min_sample_count; }
  void set_max_categories (int max_categories) { m_max_categories = max_categories; }
  void set_max_number_of_trees_in_the_forest (int max_number_of_trees_in_the_forest)
  { m_max_number_of_trees_in_the_forest = max_number_of_trees_in_the_forest; }
  void set_forest_accuracy (float forest_accuracy) { m_forest_accuracy = forest_accuracy; }
  
  
  /// @}

  /// \name Training
  /// @{
  
  /*!
    \brief Runs the training algorithm.

    From the set of provided ground truth, this algorithm estimates
    sets up the random trees that produce the most accurate result
    with respect to this ground truth.

    \pre At least one ground truth item should be assigned to each
    label.

    \param ground_truth vector of label indices. It should contain for
    each input item, in the same order as the input set, the index of
    the corresponding label in the `Label_set` provided in the
    constructor. Input items that do not have a ground truth
    information should be given the value `-1`.
  */
  template <typename LabelIndexRange>
  void train (const LabelIndexRange& ground_truth)
  {
#if (CV_MAJOR_VERSION < 3)
    if (rtree != NULL)
      delete rtree;
#endif

#ifdef CGAL_CLASSIFICATION_VERBOSE

    std::cerr << "Training random forest (OpenCV "
              << CV_MAJOR_VERSION << "."
              << CV_MINOR_VERSION << ")" << std::endl;
#endif
    
    std::size_t nb_samples = 0;
    for (std::size_t i = 0; i < ground_truth.size(); ++ i)
      if (ground_truth[i] != -1)
        ++ nb_samples;

    cv::Mat training_features (int(nb_samples), int(m_features.size()), CV_32FC1);
    cv::Mat training_labels (int(nb_samples), 1, CV_32FC1);

    for (std::size_t i = 0, index = 0; i < ground_truth.size(); ++ i)
      if (ground_truth[i] != -1)
      {
        for (std::size_t f = 0; f < m_features.size(); ++ f)
          training_features.at<float>(int(index), int(f)) = m_features[f]->value(i);
        training_labels.at<float>(int(index), 0) = static_cast<float>(ground_truth[i]);
        ++ index;
      }

#if (CV_MAJOR_VERSION < 3)
    float* priors = new float[m_labels.size()];
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      priors[i] = 1.;

    CvRTParams params (m_max_depth, m_min_sample_count,
                       0, false, m_max_categories, priors, false, 0,
                       m_max_number_of_trees_in_the_forest,
                       m_forest_accuracy,
                       CV_TERMCRIT_ITER | CV_TERMCRIT_EPS
      );

    cv::Mat var_type (m_features.size() + 1, 1, CV_8U);
    var_type.setTo (cv::Scalar(CV_VAR_NUMERICAL));

    rtree = new CvRTrees;
    rtree->train (training_features, CV_ROW_SAMPLE, training_labels,
                  cv::Mat(), cv::Mat(), var_type, cv::Mat(), params);

    delete[] priors;
#else
    rtree = cv::ml::RTrees::create();
    rtree->setMaxDepth (m_max_depth);
    rtree->setMinSampleCount (m_min_sample_count);
    rtree->setMaxCategories (m_max_categories);
    rtree->setCalculateVarImportance (false);
    rtree->setRegressionAccuracy (m_forest_accuracy);
    rtree->setUseSurrogates(false);
    rtree->setPriors(cv::Mat());
    rtree->setCalculateVarImportance(false);
    
    cv::TermCriteria criteria (cv::TermCriteria::EPS + cv::TermCriteria::COUNT, m_max_number_of_trees_in_the_forest, 0.01f);
    rtree->setTermCriteria (criteria);

    cv::Ptr<cv::ml::TrainData> tdata = cv::ml::TrainData::create
      (training_features, cv::ml::ROW_SAMPLE, training_labels);

    rtree->train (tdata);

#endif

  }

  /// @}

  /// \cond SKIP_IN_MANUAL
  void operator() (std::size_t item_index, std::vector<float>& out) const
  {
    out.resize (m_labels.size(), 0.);
    
    cv::Mat feature (1, int(m_features.size()), CV_32FC1);
    for (std::size_t f = 0; f < m_features.size(); ++ f)
      feature.at<float>(0, int(f)) = m_features[f]->value(item_index);

//compute the result of each tree
#if (CV_MAJOR_VERSION < 3)
    std::size_t nb_trees = std::size_t(rtree->get_tree_count());
    for (std::size_t i = 0; i < nb_trees; i++)
    {
      std::size_t l = rtree->get_tree(int(i))->predict(feature, cv::Mat())->value;
      out[l] += 1.;
    }

    for (std::size_t i = 0; i < out.size(); ++ i)
      out[i] = - std::log (out[i] / nb_trees);
#else

    std::vector<float> result (1, 0);
    
    rtree->predict (feature, result);
    
    for (std::size_t i = 0; i < out.size(); ++ i)
      if (i == std::size_t(result[0]))
        out[i] = 0.;
      else
        out[i] = 1.;

#endif
  }
  /// \endcond

  /// \name Input/Output
  /// @{


  /*!
    \brief Saves the current configuration in the file named `filename`.

    This allows to easily save and recover a specific classification
    configuration.

    The output file is written in an XML format that is readable by
    the `load_configuration()` method.
  */
  void save_configuration (const char* filename)
  {
    rtree->save(filename);
  }

  /*!
    \brief Loads a configuration from the file named `filename`.

    The input file should be in the XML format written by the
    `save_configuration()` method. The feature set of the classifier
    should contain the exact same features in the exact same order as
    the ones present when the file was generated using
    `save_configuration()`.
  */
  void load_configuration (const char* filename)
  {
#if (CV_MAJOR_VERSION < 3)
    if (rtree != NULL)
      delete rtree;
    rtree = new CvRTrees;
    rtree->load(filename);
#else
    rtree = cv::ml::StatModel::load<cv::ml::RTrees> (filename);
#endif
  }


};

}

}

#endif // CGAL_CLASSIFICATION_OPENCV_RANDOM_FOREST_CLASSIFIER_H
