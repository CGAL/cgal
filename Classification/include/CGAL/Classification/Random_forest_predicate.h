#ifndef CGAL_CLASSIFICATION_RANDOM_FOREST_PREDICATE_H
#define CGAL_CLASSIFICATION_RANDOM_FOREST_PREDICATE_H

#include <CGAL/Classification/Feature_set.h>
#include <CGAL/Classification/Label_set.h>

#include <cv.h>       // opencv general include file
#include <ml.h>		  // opencv machine learning include file

namespace CGAL {

namespace Classification {

/*!
  \ingroup PkgClassificationPredicates

  \brief %Classification predicate based on a random forest algorithm.

  \note This class requires the \ref thirdpartyOpenCV library.

  \cgalModels `CGAL::Classification::Predicate`
*/
class Random_forest_predicate
{
  Label_set& m_labels;
  Feature_set& m_features;
  CvRTrees* rtree;

public:
  
/*!
  \brief Instantiate the predicate using the sets of `labels` and `features`.
*/
  Random_forest_predicate (Label_set& labels,
                           Feature_set& features)
    : m_labels (labels), m_features (features), rtree (NULL)
  {  }

  /// \cond SKIP_IN_MANUAL
  ~Random_forest_predicate ()
  {
    if (rtree != NULL)
      delete rtree;
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
    if (rtree != NULL)
      delete rtree;
    
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
  }

  /// \cond SKIP_IN_MANUAL
  void operator() (std::size_t item_index, std::vector<float>& out) const
  {
    out.resize (m_labels.size(), 1.);
    
    cv::Mat feature (1, m_features.size(), CV_32FC1);
    for (std::size_t f = 0; f < m_features.size(); ++ f)
      feature.at<float>(0, f) = m_features[f]->value(item_index);

    float result = rtree->predict (feature, cv::Mat());
    std::size_t label = std::size_t(result);
    if (label < out.size())
      out[label] = 0.;
  }
  /// \endcond

};

}

}

#endif // CGAL_CLASSIFICATION_RANDOM_FOREST_PREDICATE_H
