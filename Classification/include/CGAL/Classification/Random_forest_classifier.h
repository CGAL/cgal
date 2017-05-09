#ifndef CGAL_CLASSIFICATION_RANDOM_FOREST_CLASSIFIER_H
#define CGAL_CLASSIFICATION_RANDOM_FOREST_CLASSIFIER_H

#include <CGAL/Classification/Feature_set.h>
#include <CGAL/Classification/Label_set.h>

#include <opencv2/opencv.hpp>

namespace CGAL {

namespace Classification {

/*!
  \ingroup PkgClassificationClassifiers

  \brief Classifier based on a random forest algorithm.

  \note This class requires the \ref thirdpartyOpenCV library.

  \cgalModels `CGAL::Classification::Classifier`
*/
class Random_forest_classifier
{
  Label_set& m_labels;
  Feature_set& m_features;

#if (CV_MAJOR_VERSION < 3)
  CvRTrees* rtree;
#else
  cv::Ptr<cv::ml::RTrees> rtree;
#endif

public:
  
/*!
  \brief Instantiate the classifier using the sets of `labels` and `features`.
*/
  Random_forest_classifier (Label_set& labels,
                            Feature_set& features)
    : m_labels (labels), m_features (features)
#if (CV_MAJOR_VERSION < 3)
    , rtree (NULL)
#endif
  {  }

  /// \cond SKIP_IN_MANUAL
  ~Random_forest_classifier ()
  {
#if (CV_MAJOR_VERSION < 3)
    if (rtree != NULL)
      delete rtree;
#endif
  }
  /// \endcond
  
  /*!
    \brief Runs the training algorithm.

    From the set of provided ground truth, this algorithm estimates
    sets up the random trees that produce the most accurate result
    with respect to this ground truth.

    For more details on the parameters of this algorithm, please refer
    to [the official documentation of OpenCV](http://docs.opencv.org/2.4/modules/ml/doc/random_trees.html).

    \note Each label should be assigned at least one ground truth
    item.

    \param ground_truth vector of label indices. It should contain for
    each input item, in the same order as the input set, the index of
    the corresponding label in the `Label_set` provided in the
    constructor. Input items that do not have a ground truth
    information should be given the value `std::size_t(-1)`.
    \param max_depth see [OpenCV](http://docs.opencv.org/2.4/modules/ml/doc/random_trees.html#cvrtparams-cvrtparams).
    \param min_sample_count see [OpenCV](http://docs.opencv.org/2.4/modules/ml/doc/random_trees.html#cvrtparams-cvrtparams).
    \param max_categories see [OpenCV](http://docs.opencv.org/2.4/modules/ml/doc/random_trees.html#cvrtparams-cvrtparams).
    \param max_number_of_trees_in_the_forest see [OpenCV](http://docs.opencv.org/2.4/modules/ml/doc/random_trees.html#cvrtparams-cvrtparams).
    \param forest_accuracy see [OpenCV](http://docs.opencv.org/2.4/modules/ml/doc/random_trees.html#cvrtparams-cvrtparams).
  */
  void train (const std::vector<std::size_t>& ground_truth,
              int max_depth = 20,
              int min_sample_count = 5,
              int max_categories = 15,
              int max_number_of_trees_in_the_forest = 100,
              float forest_accuracy = 0.01f)
              
  {
#if (CV_MAJOR_VERSION < 3)
    if (rtree != NULL)
      delete rtree;
#endif
    
    std::size_t nb_samples = 0;
    for (std::size_t i = 0; i < ground_truth.size(); ++ i)
      if (ground_truth[i] != std::size_t(-1))
        ++ nb_samples;


    cv::Mat training_features (nb_samples, m_features.size(), CV_32FC1);
    cv::Mat training_labels (nb_samples, 1, CV_32FC1);

    for (std::size_t i = 0, index = 0; i < ground_truth.size(); ++ i)
      if (ground_truth[i] != std::size_t(-1))
      {
        for (std::size_t f = 0; f < m_features.size(); ++ f)
          training_features.at<float>(index, f) = m_features[f]->value(i);
        training_labels.at<float>(index, 0) = ground_truth[i];
        ++ index;
      }

#if (CV_MAJOR_VERSION < 3)
    float* priors = new float[m_labels.size()];
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      priors[i] = 1.;

    CvRTParams params (max_depth, min_sample_count,
                       0, false, max_categories, priors, false, 0,
                       max_number_of_trees_in_the_forest,
                       forest_accuracy,
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
    rtree->setMaxDepth (max_depth);
    rtree->setMinSampleCount (min_sample_count);
    rtree->setMaxCategories (max_categories);
    rtree->setCalculateVarImportance (false);
    rtree->setRegressionAccuracy (forest_accuracy);
    rtree->setUseSurrogates(false);
    rtree->setPriors(cv::Mat());
    rtree->setCalculateVarImportance(false);
    
    cv::TermCriteria criteria (cv::TermCriteria::EPS + cv::TermCriteria::COUNT, max_number_of_trees_in_the_forest, 0.01f);
    rtree->setTermCriteria (criteria);

    cv::Ptr<cv::ml::TrainData> tdata = cv::ml::TrainData::create
      (training_features, cv::ml::ROW_SAMPLE, training_labels);

    rtree->train (tdata);

#endif

  }

  /// \cond SKIP_IN_MANUAL
  void operator() (std::size_t item_index, std::vector<float>& out) const
  {
    out.resize (m_labels.size(), 0.);
    
    cv::Mat feature (1, m_features.size(), CV_32FC1);
    for (std::size_t f = 0; f < m_features.size(); ++ f)
      feature.at<float>(0, f) = m_features[f]->value(item_index);

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

  void save_configuration (const char* filename)
  {
    rtree->save(filename);
  }

  void load_configuration (const char* filename)
  {
#if (CV_MAJOR_VERSION < 3)
    if (rtree != NULL)
      delete rtree;
    rtree = new CvRTrees;
#else
    rtree = cv::ml::StatModel::load<cv::ml::RTrees> (filename);
#endif
  }
  /// \endcond

};

}

}

#endif // CGAL_CLASSIFICATION_RANDOM_FOREST_CLASSIFIER_H
