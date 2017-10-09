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

#ifndef CGAL_CLASSIFICATION_ETHZ_RANDOM_FOREST_CLASSIFIER_H
#define CGAL_CLASSIFICATION_ETHZ_RANDOM_FOREST_CLASSIFIER_H

#include <CGAL/license/Classification.h>

#include <CGAL/Classification/Feature_set.h>
#include <CGAL/Classification/Label_set.h>

#ifdef CGAL_CLASSIFICATION_VERBOSE
#define VERBOSE_TREE_PROGRESS true
#endif

#include <CGAL/Classification/internal/auxiliary/random-forest/node-gini.hpp>
#include <CGAL/Classification/internal/auxiliary/random-forest/forest.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>


namespace CGAL {

namespace Classification {

/*!
  \ingroup PkgClassificationClassifiers

  \brief %Classifier based on the ETHZ version of random forest algorithm.

  \cgalModels `CGAL::Classification::Classifier`
*/
class ETHZ_random_forest_classifier
{
  typedef liblearning::RandomForest::RandomForest
  < liblearning::RandomForest::NodeGini
    < liblearning::RandomForest::AxisAlignedSplitter> > Forest;
  
  const Label_set& m_labels;
  const Feature_set& m_features;
  Forest* m_rfc;

public:
  
/*!
  \brief Instantiate the classifier using the sets of `labels` and `features`.

*/
  ETHZ_random_forest_classifier (const Label_set& labels,
                                 const Feature_set& features)
    : m_labels (labels), m_features (features), m_rfc (NULL)
  { }
  
  /// \cond SKIP_IN_MANUAL
  ~ETHZ_random_forest_classifier ()
  {
    if (m_rfc != NULL)
      delete m_rfc;
  }
  /// \endcond
  
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
  void train (const LabelIndexRange& ground_truth,
              bool reset_trees = true,
              std::size_t num_trees = 25,
              std::size_t max_depth = 20)
  {
    liblearning::RandomForest::ForestParams params;
    params.n_trees   = num_trees;
    params.max_depth = max_depth;

    std::vector<int> gt;
    std::vector<float> ft;
    
    for (std::size_t i = 0; i < ground_truth.size(); ++ i)
      if (ground_truth[i] != std::size_t(-1))
      {
        for (std::size_t f = 0; f < m_features.size(); ++ f)
          ft.push_back(m_features[f]->value(i));
        gt.push_back(ground_truth[i]);
      }

    liblearning::DataView2D<int> label_vector (&(gt[0]), gt.size(), 1);    
    liblearning::DataView2D<float> feature_vector(&(ft[0]), gt.size(), ft.size() / gt.size());

    if (m_rfc != NULL && reset_trees)
    {
      delete m_rfc;
      m_rfc = NULL;
    }
    
    if (m_rfc == NULL)
      m_rfc = new Forest (params);

    liblearning::RandomForest::AxisAlignedRandomSplitGenerator generator;
    
    m_rfc->train(feature_vector, label_vector, liblearning::DataView2D<int>(), generator, 0, false, reset_trees);
  }

  /// \cond SKIP_IN_MANUAL
  void operator() (std::size_t item_index, std::vector<float>& out) const
  {
    out.resize (m_labels.size(), 0.);
    
    std::vector<float> ft;
    ft.reserve (m_features.size());
    for (std::size_t f = 0; f < m_features.size(); ++ f)
      ft.push_back (m_features[f]->value(item_index));

    std::vector<float> prob (m_labels.size());

    m_rfc->evaluate (ft.data(), prob.data());
    
    for (std::size_t i = 0; i < out.size(); ++ i)
      out[i] = - std::log (prob[i]);
  }

  void save_configuration (const char* filename)
  {
    std::ofstream ofs(filename, std::ios_base::out | std::ios_base::binary);
    boost::iostreams::filtering_ostream outs;
    outs.push(boost::iostreams::gzip_compressor());
    outs.push(ofs);
    boost::archive::text_oarchive oas(outs);
    oas << BOOST_SERIALIZATION_NVP(*m_rfc);
  }

  void load_configuration (const char* filename,
                           std::size_t num_trees = 25,
                           std::size_t max_depth = 20)
  {
    liblearning::RandomForest::ForestParams params;
    params.n_trees   = num_trees;
    params.max_depth = max_depth;
    if (m_rfc != NULL)
      delete m_rfc;
    m_rfc = new Forest (params);
    
    std::ifstream ifs(filename, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_istream ins;
    ins.push(boost::iostreams::gzip_decompressor());
    ins.push(ifs);
    boost::archive::text_iarchive ias(ins);
    ias >> BOOST_SERIALIZATION_NVP(*m_rfc);
  }
  /// \endcond

};

}

}

#endif // CGAL_CLASSIFICATION_ETHZ_RANDOM_FOREST_CLASSIFIER_H
