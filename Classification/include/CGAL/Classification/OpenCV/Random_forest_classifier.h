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

#ifndef CGAL_CLASSIFICATION_OPENCV_RANDOM_FOREST_CLASSIFIER_H
#define CGAL_CLASSIFICATION_OPENCV_RANDOM_FOREST_CLASSIFIER_H

#include <CGAL/license/Classification.h>

#include <CGAL/Classification/Feature_set.h>
#include <CGAL/Classification/Label_set.h>

#include <opencv2/opencv.hpp>

//In opencv version 2.X the first digit is named EPOCH,
//until version 3.0 where EPOCH disappears and it becomes MAJOR. Hence this
//weird condition
#ifdef CV_VERSION_EPOCH
  #if  CV_VERSION_MAJOR == 4 && CV_VERSION_MINOR>= 11
    #include <opencv2/ml.hpp>
  #else
    #include <opencv2/ml/ml.hpp>
  #endif
#else
  #include <opencv2/ml.hpp>
#endif



namespace CGAL {

namespace Classification {

namespace OpenCV {

/*!
  \ingroup PkgClassificationClassifiersOpenCV

  \brief %Classifier based on the OpenCV version of the random forest algorithm.

  \note This class requires the \ref thirdpartyOpenCV library.

  \cgalModels `CGAL::Classification::Classifier`
*/
class Random_forest_classifier
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
  \brief instantiates the classifier using the sets of `labels` and `features`.

  Parameters documentation is copy-pasted from [the official documentation of OpenCV](https://docs.opencv.org/2.4/modules/ml/doc/random_trees.html). For more details on this method, please refer to it.

  \param labels label set used.
  \param features feature set used.
  \param max_depth the depth of the tree. A low value will likely underfit and conversely a high value will likely overfit. The optimal value can be obtained using cross validation or other suitable methods.
  \param min_sample_count minimum samples required at a leaf node for it to be split. A reasonable value is a small percentage of the total data e.g. 1%.
  \param max_categories Cluster possible values of a categorical variable into \f$ K \leq max\_categories \f$ clusters to find a suboptimal split. If a discrete variable, on which the training procedure tries to make a split, takes more than max_categories values, the precise best subset estimation may take a very long time because the algorithm is exponential. Instead, many decision trees engines (including ML) try to find sub-optimal split in this case by clustering all the samples into max_categories clusters that is some categories are merged together. The clustering is applied only in \f$ n>2-class \f$ classification problems for categorical variables with \f$ N > max\_categories \f$ possible values. In case of regression and 2-class classification the optimal split can be found efficiently without employing clustering, thus the parameter is not used in these cases.
  \param max_number_of_trees_in_the_forest The maximum number of trees in the forest (surprise, surprise). Typically the more trees you have the better the accuracy. However, the improvement in accuracy generally diminishes and asymptotes pass a certain number of trees. Also to keep in mind, the number of tree increases the prediction time linearly.
  \param forest_accuracy Sufficient accuracy (OOB error).
*/
  Random_forest_classifier (const Label_set& labels,
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
    , rtree (nullptr)
#endif
  {  }

  /// \cond SKIP_IN_MANUAL
  ~Random_forest_classifier ()
  {
#if (CV_MAJOR_VERSION < 3)
    if (rtree != nullptr)
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
    \brief runs the training algorithm.

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
    CGAL_precondition (m_labels.is_valid_ground_truth (ground_truth));

#if (CV_MAJOR_VERSION < 3)
    if (rtree != nullptr)
      delete rtree;
#endif

#ifdef CGAL_CLASSIFICATION_VERBOSE

    std::cerr << "Training random forest (OpenCV "
              << CV_MAJOR_VERSION << "."
              << CV_MINOR_VERSION << ")" << std::endl;
#endif

    std::size_t nb_samples = 0;
    for (const auto& gt_value : ground_truth)
      if (int(gt_value) != -1)
        ++ nb_samples;

    cv::Mat training_features (int(nb_samples), int(m_features.size()), CV_32FC1);
    cv::Mat training_labels (int(nb_samples), 1, CV_32FC1);

    std::size_t i = 0, index = 0;
    for (const auto& gt_value : ground_truth)
    {
      if (int(gt_value) != -1)
      {
        for (std::size_t f = 0; f < m_features.size(); ++ f)
          training_features.at<float>(int(index), int(f)) = m_features[f]->value(i);
        training_labels.at<float>(int(index), 0) = static_cast<float>(gt_value);
        ++ index;
      }
      ++ i;
    }


#if (CV_MAJOR_VERSION < 3)
    float* priors = new float[m_labels.size()];
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      priors[i] = 1.;

    CvRTParams params;

    if (m_forest_accuracy == 0.f)
      params = CvRTParams
        (m_max_depth, m_min_sample_count,
         0, false, m_max_categories, priors, false, 0,
         m_max_number_of_trees_in_the_forest,
         m_forest_accuracy, CV_TERMCRIT_ITER);
    else
      params = CvRTParams
        (m_max_depth, m_min_sample_count,
         0, false, m_max_categories, priors, false, 0,
         m_max_number_of_trees_in_the_forest,
         m_forest_accuracy, CV_TERMCRIT_EPS | CV_TERMCRIT_ITER);

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

    cv::TermCriteria criteria;
    if (m_forest_accuracy == 0.f)
      criteria = cv::TermCriteria (cv::TermCriteria::COUNT, m_max_number_of_trees_in_the_forest, m_forest_accuracy);
    else
      criteria = cv::TermCriteria (cv::TermCriteria::EPS + cv::TermCriteria::COUNT, m_max_number_of_trees_in_the_forest, m_forest_accuracy);

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
      out[i] = out[i] / nb_trees;
#else

    std::vector<float> result (1, 0);

    rtree->predict (feature, result);

    for (std::size_t i = 0; i < out.size(); ++ i)
      if (i == std::size_t(result[0]))
        out[i] = 1.f;
      else
        out[i] = 0.f;

#endif
  }
  /// \endcond

  /// \name Input/Output
  /// @{


  /*!
    \brief saves the current configuration in the file named `filename`.

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
    \brief loads a configuration from the file named `filename`.

    The input file should be in the XML format written by the
    `save_configuration()` method. The feature set of the classifier
    should contain the exact same features in the exact same order as
    the ones present when the file was generated using
    `save_configuration()`.
  */
  void load_configuration (const char* filename)
  {
#if (CV_MAJOR_VERSION < 3)
    if (rtree != nullptr)
      delete rtree;
    rtree = new CvRTrees;
    rtree->load(filename);
#else
    rtree = cv::ml::StatModel::load<cv::ml::RTrees> (filename);
#endif
  }

  /// @}

};

}

/// \cond SKIP_IN_MANUAL
// Backward compatibility
typedef OpenCV::Random_forest_classifier OpenCV_random_forest_classifier;
/// \endcond

}

}

#endif // CGAL_CLASSIFICATION_OPENCV_RANDOM_FOREST_CLASSIFIER_H
