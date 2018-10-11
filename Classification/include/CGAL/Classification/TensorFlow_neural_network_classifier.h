// Copyright (c) 2018 GeometryFactory Sarl (France).
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

#ifndef CGAL_CLASSIFICATION_TENSORFLOW_NEURAL_NETWORK_CLASSIFIER_H
#define CGAL_CLASSIFICATION_TENSORFLOW_NEURAL_NETWORK_CLASSIFIER_H

#include <CGAL/license/Classification.h>

#include <CGAL/Classification/Feature_set.h>
#include <CGAL/Classification/Label_set.h>

#include <tensorflow/cc/client/client_session.h>
#include <tensorflow/cc/ops/standard_ops.h>
#include <tensorflow/core/framework/tensor.h>
#include <tensorflow/core/framework/types.h>
#include <tensorflow/cc/framework/gradients.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <CGAL/boost/iterator/counting_iterator.hpp>

#define CGAL_NEURAL_NETWORKS_VERSION "0.1"

namespace CGAL {

namespace Classification {

namespace TF = tensorflow;
namespace TFops = tensorflow::ops;

/*!
  \ingroup PkgClassificationClassifiers

  \cgalModels `CGAL::Classification::Classifier`
*/
class TensorFlow_neural_network_classifier
{
#define USE_RELU
#if defined(USE_TANH)
  typedef TFops::Tanh Activation_function;
#elif defined(USE_RELU)
  typedef TFops::Relu Activation_function;
#elif defined(USE_RELU6)
  typedef TFops::Relu6 Activation_function;
#elif defined(USE_ELU)
  typedef TFops::Elu Activation_function;
#elif defined(USE_SELU)
  typedef TFops::Selu Activation_function;
#elif defined(USE_SIGMOID)
  typedef TFops::Sigmoid Activation_function;
#endif

  bool m_verbose;
  
  const Label_set& m_labels;
  const Feature_set& m_features;
  std::vector<float> m_feature_means;
  std::vector<float> m_feature_sd;

  TF::Scope* m_root;
  TFops::Placeholder* m_ph_ft;
  TFops::Placeholder* m_ph_gt;

//#define USE_IOU
#ifdef USE_IOU
  TFops::ReduceMean* m_loss;
#else
  TFops::ReduceMean* m_loss;
#endif
  
  std::vector<std::vector<float> > m_weights;
  std::vector<std::vector<float> > m_bias;
  float m_learning_rate;
  std::vector<TF::Output> m_layers;

#define USE_ADAM
#ifdef USE_ADAM
  std::vector<TFops::ApplyAdam> m_gradient_descent;
#else
  std::vector<TFops::ApplyGradientDescent> m_gradient_descent;
#endif
  
  TF::ClientSession* m_session;

  
public:

  /// \cond SKIP_IN_MANUAL
  const Feature_set& features() const { return m_features; }
  /// \endcond
  
  /// \name Constructor
  /// @{
  
  /*!
    \brief Instantiate the classifier using the sets of `labels` and `features`.

  */
  TensorFlow_neural_network_classifier (const Label_set& labels,
                                        const Feature_set& features)
    : m_verbose (true), m_labels (labels), m_features (features)
    , m_root (NULL), m_ph_ft (NULL), m_ph_gt (NULL), m_loss(NULL), m_session (NULL)
  { }
  
  /// \cond SKIP_IN_MANUAL
  ~TensorFlow_neural_network_classifier ()
  {
    clear();
  }

  bool initialized() const { return (m_root != NULL); }

  bool& verbose() { return m_verbose; }

  void clear()
  {
    clear (m_ph_gt);
    clear (m_ph_ft);
    clear (m_loss);
    clear (m_session);
    clear (m_root);
    m_weights.clear();
    m_bias.clear();
    m_gradient_descent.clear();
    m_layers.clear();
  }

  template <typename T>
  void clear (T* t)
  {
    if (t != NULL)
      delete t;
    t = NULL;
  }

  void compute_normalization_coefficients (const std::vector<std::size_t>& indices)
  {
//#define DO_NOT_NORMALIZE_FEATURES
#ifdef DO_NOT_NORMALIZE_FEATURES
    m_feature_means.clear();
    m_feature_means.resize (m_features.size(), 0.f);
    m_feature_sd.clear();
    m_feature_sd.resize (m_features.size(), 1.f);
#else
    m_feature_means.clear();
    m_feature_means.resize (m_features.size(), 0.f);
    
    for (std::size_t i = 0; i < indices.size(); ++ i)
      for (std::size_t f = 0; f < m_features.size(); ++ f)
        m_feature_means[f] += m_features[f]->value(indices[i]);
    for (std::size_t f = 0; f < m_features.size(); ++ f)
      m_feature_means[f] /= indices.size();

    m_feature_sd.clear();
    m_feature_sd.resize (m_features.size(), 0.f);

    for (std::size_t i = 0; i < indices.size(); ++ i)
      for (std::size_t f = 0; f < m_features.size(); ++ f)
        m_feature_sd[f] +=
          (m_features[f]->value(indices[i]) - m_feature_means[f]) *
          (m_features[f]->value(indices[i]) - m_feature_means[f]);
    for (std::size_t f = 0; f < m_features.size(); ++ f)
    {
      m_feature_sd[f] = std::sqrt (m_feature_sd[f] / indices.size());
      if (m_feature_sd[f] == 0.f)
        m_feature_sd[f] = 1.f;

      if (std::isnan(m_feature_means[f]))
        m_feature_means[f] = 0.f;
      if (std::isnan(m_feature_sd[f]))
        m_feature_sd[f] = 1.f;

      // if (m_verbose)
      //   std::cerr << "#" << f << ": " << m_features[f]->name() << " = "
      //             << m_feature_means[f] << " +/- " << m_feature_sd[f] << std::endl;
    }
#endif
  }
  /// \endcond
  
  /// @}

  /// \name Training
  
  /// @{
  /*!
    \brief Runs the training algorithm.

  */
  template <typename LabelIndexRange>
  void train (const LabelIndexRange& ground_truth,
              bool restart_from_scratch = false,
              std::size_t number_of_iterations = 5000,
              float learning_rate = 0.1,
              std::size_t batch_size = 1000,
              const std::vector<std::size_t>& hidden_layers
              = std::vector<std::size_t>())
  {
    if (restart_from_scratch)
      clear();

#define EQUAL_DISTRIBUTION

#ifdef EQUAL_DISTRIBUTION
    std::vector<std::vector<std::size_t> > random_indices (m_labels.size());
#endif
    
    std::vector<std::size_t> indices;
    std::vector<int> raw_gt;
    for (std::size_t i = 0; i < ground_truth.size(); ++ i)
    {
      int gc = int(ground_truth[i]);
      if (gc != -1)
      {
        indices.push_back (i);
        raw_gt.push_back (gc);
#ifdef EQUAL_DISTRIBUTION
        random_indices[std::size_t(gc)].push_back (indices.size() - 1);
#endif
      }
    }

    if (!initialized())
    {
      m_learning_rate = learning_rate;
      build_architecture (hidden_layers);
      if (m_verbose) std::cerr << "II - NORMALIZING FEATURES" << std::endl;
      compute_normalization_coefficients (indices);
    }
    
    if (m_verbose) std::cerr << "III - TRAINING NEURAL NETWORK" << std::endl;
    
#ifndef EQUAL_DISTRIBUTION
    batch_size = std::min (std::size_t(batch_size), indices.size());
    std::vector<std::size_t> random_indices
      (boost::make_counting_iterator<std::size_t>(0),
       boost::make_counting_iterator<std::size_t>(indices.size()));
#endif
    
    std::vector<TF::Output> operations;
    for (std::size_t i = 0; i < m_gradient_descent.size(); ++ i)
      operations.push_back (m_gradient_descent[i]);
    operations.push_back (m_layers.back());
    
    std::vector<TF::Tensor> outputs;

    m_weights.clear();
    m_bias.clear();
    // std::ofstream log ("loss.log");
    // log.precision(18);

    for (std::size_t i = 0; i < number_of_iterations; ++ i)
    {
#ifdef EQUAL_DISTRIBUTION
      std::size_t total_batch_size = 0;
      for (std::size_t j = 0; j < random_indices.size(); ++ j)
      {
        if (random_indices[j].empty())
          continue;
        std::size_t local_batch_size = std::min (std::size_t(batch_size / m_labels.size()),
                                                 random_indices[j].size());
        random_unique (random_indices[j].begin(), random_indices[j].end(), local_batch_size);
        total_batch_size += local_batch_size;

        if (i == 0)
        {
          if (m_verbose) std::cerr << "   * Size of batch for " << m_labels[j]->name() << ": "
                                  << local_batch_size << std::endl;
        }
      }
      
      TF::Tensor ft
        (TF::DataTypeToEnum<float>::v(), 
         TF::TensorShape {(long long)(total_batch_size), (long long)(m_features.size())});
      TF::Tensor gt
        (TF::DataTypeToEnum<float>::v(), 
         TF::TensorShape {(long long)(total_batch_size), (long long)(m_labels.size())});

      float* ft_data = ft.flat<float>().data();
      float* gt_data = gt.flat<float>().data();

      // Fill input tensors
      std::size_t current_i = 0;
      for (std::size_t j = 0; j < random_indices.size(); ++ j)
      {
        if (random_indices[j].empty())
          continue;
        std::size_t local_batch_size = std::min (std::size_t(batch_size / m_labels.size()),
                                                 random_indices[j].size());
        for (std::size_t k = 0; k < local_batch_size; ++ k)
        {
          std::size_t idx = indices[random_indices[j][k]];
          int g = raw_gt[random_indices[j][k]];

          for (std::size_t f = 0; f < m_features.size(); ++ f)
            ft_data[current_i * m_features.size() + f]
              = (m_features[f]->value(idx) - m_feature_means[f]) / m_feature_sd[f];

          for (std::size_t l = 0; l < m_labels.size(); ++ l)
            if (std::size_t(g) == l)
              gt_data[current_i * m_labels.size() + l] = 1.f;
            else
              gt_data[current_i * m_labels.size() + l] = 0.f;

          ++ current_i;
        }
      }
#else
      random_unique (random_indices.begin(), random_indices.end(), batch_size);

      TF::Tensor ft
        (TF::DataTypeToEnum<float>::v(), 
         TF::TensorShape {(long long)(batch_size), (long long)(m_features.size())});
      TF::Tensor gt
        (TF::DataTypeToEnum<float>::v(), 
         TF::TensorShape {(long long)(batch_size), (long long)(m_labels.size())});

      float* ft_data = ft.flat<float>().data();
      float* gt_data = gt.flat<float>().data();

      // Fill input tensors
      for (std::size_t i = 0; i < batch_size; ++ i)
      {
        std::size_t idx = indices[random_indices[i]];
        int g = raw_gt[random_indices[i]];

        for (std::size_t f = 0; f < m_features.size(); ++ f)
          ft_data[i * m_features.size() + f]
            = (m_features[f]->value(idx) - m_feature_means[f]) / m_feature_sd[f];

        for (std::size_t l = 0; l < m_labels.size(); ++ l)
          if (std::size_t(g) == l)
            gt_data[i * m_labels.size() + l] = 1.f;
          else
            gt_data[i * m_labels.size() + l] = 0.f;
      }
#endif
      
      TF_CHECK_OK (m_session->Run ({{*m_ph_ft, ft}, {*m_ph_gt, gt}}, {*m_loss}, &outputs));

//      log << outputs[0].scalar<float>() << std::endl;
      if ((i+1) % (number_of_iterations / 20) == 0)
      {
        if (m_verbose) std::cerr << "   * Step " << i+1 << "/" << number_of_iterations << ": loss = "
                                << outputs[0].scalar<float>() << std::endl;
      }
      if (!std::isfinite(*outputs[0].scalar<float>().data()))
      {
        if (m_verbose)
          std::cerr << "Loss is " << outputs[0].scalar<float>() << ", aborting" << std::endl;
        return;
      }

      if (i == number_of_iterations - 1)
      {
        TF_CHECK_OK(m_session->Run({{*m_ph_ft, ft}, {*m_ph_gt, gt}}, operations, &outputs));
        for (std::size_t o = 0; o < outputs.size() - 1; ++ o)
        {
          std::size_t size = outputs[o].dim_size(0) * outputs[o].dim_size(1);
          if (o < outputs.size() / 2)
          {
            m_weights.push_back (std::vector<float>(size));
            float* weights_data = outputs[o].flat<float>().data();
            std::copy_n (weights_data, size, m_weights.back().begin());
          }
          else
          {
            m_bias.push_back (std::vector<float>(size));
            float* bias_data = outputs[o].flat<float>().data();
            std::copy_n (bias_data, size, m_bias.back().begin());
          }
        }
      }
      else
        TF_CHECK_OK(m_session->Run({{*m_ph_ft, ft}, {*m_ph_gt, gt}}, operations, nullptr));
    }

  }

  /// \cond SKIP_IN_MANUAL
  void operator() (std::size_t item_index, std::vector<float>& out) const
  {
    out.resize (m_labels.size(), 0.);
    
    TF::Tensor ft
      (TF::DataTypeToEnum<float>::v(), 
       TF::TensorShape {(long long)(1), (long long)(m_features.size())});

    float* ft_data = ft.flat<float>().data();

    // Fill input tensor
    for (std::size_t f = 0; f < m_features.size(); ++ f)
      ft_data[f] = (m_features[f]->value(item_index) - m_feature_means[f]) / m_feature_sd[f];
    
    std::vector<TF::Tensor> outputs;
    TF_CHECK_OK(m_session->Run({{*m_ph_ft, ft}}, {m_layers.back()}, &outputs));

    float* output_data = outputs[0].flat<float>().data();

    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      out[i] = output_data[i];
  }

  
  void operator() (const std::vector<std::size_t>& item_indices,
                   std::vector<std::vector<float> >& out) const
  {
    out.resize (item_indices.size(), std::vector<float>(m_labels.size(), 0.));
    
    TF::Tensor ft
      (TF::DataTypeToEnum<float>::v(), 
       TF::TensorShape {(long long)(item_indices.size()), (long long)(m_features.size())});

    float* ft_data = ft.flat<float>().data();

    // Fill input tensor
    for (std::size_t i = 0; i < item_indices.size(); ++ i)
      for (std::size_t f = 0; f < m_features.size(); ++ f)
        ft_data[i * m_features.size() + f]
          = (m_features[f]->value(item_indices[i]) - m_feature_means[f]) / m_feature_sd[f];
    
    std::vector<TF::Tensor> outputs;
    TF_CHECK_OK(m_session->Run({{*m_ph_ft, ft}}, {m_layers.back()}, &outputs));

    float* output_data = outputs[0].flat<float>().data();

    for (std::size_t i = 0; i < item_indices.size(); ++ i)
      for (std::size_t l = 0; l < m_labels.size(); ++ l)
        out[i][l] = output_data[i * m_labels.size() + l];
  }
  /// \endcond
  
  /// @}

  /// \name Input/Output
  /// @{
  
  /*!
    \brief Saves the current configuration in the stream `output`.

    This allows to easily save and recover a specific classification
    configuration.

    The output file is written in an GZIP container that is readable
    by the `load_configuration()` method.
  */
  void save_configuration (std::ostream& output)
  {
    boost::property_tree::ptree tree;

    {
      boost::property_tree::ptree ptr;
      ptr.put("classifier_version", CGAL_NEURAL_NETWORKS_VERSION);
      ptr.put("number_of_features", m_features.size());
      ptr.put("number_of_labels", m_labels.size());
      ptr.put("number_of_hidden_layers", m_layers.size() - 1);
      ptr.put("learning_rate", m_learning_rate);
      tree.add_child("classification.metadata", ptr);
    }
    
    for (std::size_t i = 0; i < m_features.size(); ++ i)
    {
      boost::property_tree::ptree ptr;
      ptr.put("name", m_features[i]->name());
      ptr.put("mean", m_feature_means[i]);
      ptr.put("standard_deviation", m_feature_sd[i]);
      tree.add_child("classification.features.feature", ptr);
    }

    for (std::size_t i = 0; i < m_labels.size(); ++ i)
    {
      boost::property_tree::ptree ptr;
      ptr.put("name", m_labels[i]->name());
      tree.add_child("classification.labels.label", ptr);
    }

    for (std::size_t i = 0; i < m_weights.size(); ++ i)
    {
      boost::property_tree::ptree ptr;
      ptr.put("height", (m_weights[i].size() / m_bias[i].size()));
      ptr.put("width", m_bias[i].size());

      std::ostringstream oss;
      oss.precision(std::numeric_limits<float>::digits10 + 2);
      for (std::size_t j = 0; j < m_weights[i].size(); ++ j)
      {
        oss << m_weights[i][j];
        if (j != m_weights[i].size() - 1)
          oss << " ";
      }
      ptr.put("weights", oss.str());
      
      oss = std::ostringstream();
      oss.precision(std::numeric_limits<float>::digits10 + 2);
      for (std::size_t j = 0; j < m_bias[i].size(); ++ j)
      {
        oss << m_bias[i][j];
        if (j != m_bias[i].size() - 1)
          oss << " ";
      }
      ptr.put("bias", oss.str());
      
      tree.add_child("classification.layers.layer", ptr);
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
    \brief Loads a configuration from the stream `input`.

    The input file should be a GZIP container written by the
    `save_configuration()` method. The feature set of the classifier
    should contain the exact same features in the exact same order as
    the ones present when the file was generated using
    `save_configuration()`.
  */
  bool load_configuration (std::istream& input, bool verbose = false)
  {
    clear();
    
    bool out = true;

    boost::property_tree::ptree tree;
    boost::property_tree::read_xml(input, tree);

    std::string version = tree.get<std::string>("classification.metadata.classifier_version");
    if (version != CGAL_NEURAL_NETWORKS_VERSION)
    {
      if (verbose)
        std::cerr << "Error: CGAL Neural Network version mismatch " << version << "/" << CGAL_NEURAL_NETWORKS_VERSION << std::endl;
      return false;
    }
    std::size_t nb_features = std::size_t(tree.get<int>("classification.metadata.number_of_features"));
    if (nb_features != m_features.size())
    {
      if (verbose)
        std::cerr << "Error: number of features mismatch " << nb_features << "/" << m_features.size() << std::endl;
      return false;
    }

    std::size_t nb_labels = std::size_t(tree.get<int>("classification.metadata.number_of_labels"));
    if (nb_labels != m_labels.size())
    {
      if (verbose)
        std::cerr << "Error: number of labels mismatch " << nb_labels << "/" << m_labels.size() << std::endl;
      return false;
    }

    std::size_t nb_layers = std::size_t(tree.get<int>("classification.metadata.number_of_hidden_layers")) + 1;
    m_learning_rate = tree.get<float>("classification.metadata.learning_rate");

    std::size_t idx = 0;
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, tree.get_child("classification.features"))
    {
      std::string name = v.second.get<std::string>("name");

      if (name != m_features[idx]->name())
      {
        if (verbose)
          std::cerr << "Warning: feature mismatch " << name << "/" << m_features[idx]->name() << std::endl;
        out = false;
      }

      m_feature_means.push_back(v.second.get<float>("mean"));
      m_feature_sd.push_back(v.second.get<float>("standard_deviation"));
      
      ++ idx;
    }

    if (idx != nb_features)
    {
      if (verbose)
        std::cerr << "Error: number of features mismatch " << nb_features << "/" << idx << std::endl;
      return false;
    }

    idx = 0;
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, tree.get_child("classification.labels"))
    {
      std::string name = v.second.get<std::string>("name");
      if (name != m_labels[idx]->name())
      {
        if (verbose)
          std::cerr << "Warning: label mismatch " << name << "/" << m_labels[idx]->name() << std::endl;
        out = false;
      }
      ++ idx;
    }
    
    if (idx != nb_labels)
    {
      if (verbose)
        std::cerr << "Error: number of labels mismatch " << nb_labels << "/" << idx << std::endl;
      return false;
    }

    idx = 0;
    m_weights.resize (nb_layers);
    m_bias.resize (nb_layers);

    std::vector<std::size_t> hidden_layers;
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, tree.get_child("classification.layers"))
    {
      if (idx >= nb_layers)
      {
        if (verbose)
          std::cerr << "Error: number of layers mismatch " << nb_layers << "/" << idx << std::endl;
        return false;
      }
      std::size_t height = std::size_t(v.second.get<int>("height"));
      std::size_t width = std::size_t(v.second.get<int>("width"));
      std::size_t size = height * width;
      m_weights[idx].reserve (size);
      m_bias[idx].reserve (width);

      std::string weights = v.second.get<std::string>("weights");
      std::istringstream iss (weights);
      float f;
      while (iss >> f)
        m_weights[idx].push_back (f);

      std::string bias = v.second.get<std::string>("bias");
      iss = std::istringstream (bias);
      while (iss >> f)
        m_bias[idx].push_back (f);

      if (idx != nb_layers - 1)
        hidden_layers.push_back (width);
      
      ++ idx;
    }
    
    if (idx != nb_layers)
    {
      if (verbose)
        std::cerr << "Error: number of layers mismatch " << nb_layers << "/" << idx << std::endl;
      return false;
    }

    build_architecture(hidden_layers);

    return out;
  }


private:

  template<class BidiIter >
  BidiIter random_unique(BidiIter begin, BidiIter end, size_t num_random) {
    size_t left = std::distance(begin, end);
    while (num_random--) {
      BidiIter r = begin;
      std::advance(r, rand()%left);
      std::swap(*begin, *r);
      ++begin;
      --left;
    }
    return begin;
  }

  void build_architecture (const std::vector<std::size_t>& hidden_layers)
  {
    m_root = new TF::Scope (TF::Scope::NewRootScope());

    
    if (m_verbose) std::cerr << "I - BUILDING NEURAL NETWORK ARCHITECTURE" << std::endl;
    
    // Get layers and sizes or init with default values
    std::vector<std::size_t> hl = hidden_layers;
    if (hl.empty())
    {
      hl.push_back (m_features.size());
      hl.push_back ((m_features.size() + m_labels.size()) / 2);
    }
    
    if (m_verbose) std::cerr << " 1) Initializing architecture:" << std::endl
                            << "   * Layer 0: " << m_features.size() << " neuron(s) (input features)" << std::endl;
      
    if (m_verbose)
      for (std::size_t i = 0; i < hl.size(); ++ i)
        std::cerr << "   * Layer " << i+1 << ": " << hl[i] << " neuron(s)" << std::endl;
    
    if (m_verbose) std::cerr << "   * Layer " << hl.size() + 1 << ": "
                            << m_labels.size() << " neuron(s) (output labels)" << std::endl;


    
    m_ph_ft = new TFops::Placeholder (*m_root, TF::DT_FLOAT);
    m_ph_gt = new TFops::Placeholder (*m_root, TF::DT_FLOAT);
    
    if (m_verbose) std::cerr << " 2) Creating weight matrices and bias" << std::endl;
    
    // create weight matrices and bias for each layer transition
    std::vector<TFops::Variable> weights;
    std::vector<TFops::Variable> bias;
    std::vector<TFops::Assign> assign_weights;
    std::vector<TFops::Assign> assign_bias;

    for (std::size_t i = 0; i <= hl.size(); ++ i)
    {
      long long size_from = 0;
      long long size_to = 0;
      if (i == 0)
      {
        size_from = m_features.size();
        size_to = hl[0];
      }
      else if (i == hl.size())
      {
        size_from = hl.back();
        size_to = m_labels.size();
      }
      else
      {
        size_from = hl[i-1];
        size_to = hl[i];
      }

      if (m_verbose) std::cerr << "   * Weight matrix " << i << " [" << size_from << ";" << size_to << "]" << std::endl;
      weights.push_back (TFops::Variable (*m_root, { size_from, size_to }, TF::DT_FLOAT));

      if (m_weights.empty())
      {
#ifdef USE_TANH // Weights initialized by Xavier method
        assign_weights.push_back (TFops::Assign (*m_root, weights.back(),
                                                 TFops::Mul (*m_root, TFops::Const(*m_root,
                                                                                   std::sqrt (1.f / (float)(size_from))),
                                                             TFops::RandomNormal (*m_root, { size_from, size_to}, TF::DT_FLOAT))));
#elif defined(USE_RELU) || defined(USE_RELU6) // Weights initialized by He method
        assign_weights.push_back (TFops::Assign (*m_root, weights.back(),
          TFops::Mul (*m_root, TFops::Const(*m_root,
                                            std::sqrt (2.f / (float)(size_from))),
                      TFops::RandomNormal (*m_root, { size_from, size_to}, TF::DT_FLOAT))));
#else // Default: weights truncated normal
        assign_weights.push_back (TFops::Assign (*m_root, weights.back(),
                                                 TFops::TruncatedNormal (*m_root, { size_from, size_to}, TF::DT_FLOAT)));
#endif
      }
      else
      {
        TF::Tensor init_weight
          (TF::DataTypeToEnum<float>::v(), 
           TF::TensorShape {(long long)(size_from), (long long)(size_to)});
        float* init_weight_data = init_weight.flat<float>().data();
        std::copy_n (m_weights[i].begin(), size_from * size_to, init_weight_data);

        assign_weights.push_back (TFops::Assign (*m_root, weights.back(), init_weight));
      }
      
      if (m_verbose) std::cerr << "   * Bias " << i << " [" << size_to << "]" << std::endl;
      bias.push_back (TFops::Variable (*m_root, { (long long)(1), size_to }, TF::DT_FLOAT));

      TF::Tensor init_bias
        (TF::DataTypeToEnum<float>::v(), 
         TF::TensorShape {(long long)(1), (long long)(size_to)});
      float* init_bias_data = init_bias.flat<float>().data();
      if (m_bias.empty())
        for (std::size_t s = 0; s < size_to; ++ s)
          init_bias_data[s] = 0.f;
      else
        std::copy_n (m_bias[i].begin(), size_to, init_bias_data);
      
      assign_bias.push_back (TFops::Assign (*m_root, bias.back(), init_bias));
    }

    if (!m_root->status().ok())
    {
      if (m_verbose) std::cerr << "Error: " << m_root->status().ToString() << std::endl;
      return;
    }
    
    if (m_verbose) std::cerr << " 3) Creating layers" << std::endl;
    
    for (std::size_t i = 0; i < weights.size(); ++ i)
    {
      if (i == 0)
        m_layers.push_back (Activation_function (*m_root, TFops::Add
                                                 (*m_root, TFops::MatMul (*m_root, *m_ph_ft, weights[0]),
                                                  bias[0])));
      else if (i == weights.size() - 1)
        m_layers.push_back (TFops::Softmax (*m_root, TFops::Add
                                            (*m_root, TFops::MatMul (*m_root, m_layers.back(), weights[i]),
                                             bias[i])));
      else
        m_layers.push_back (Activation_function (*m_root, TFops::Add
                                                 (*m_root, TFops::MatMul (*m_root, m_layers.back(), weights[i]),
                                                  bias[i])));
    }

    if (!m_root->status().ok())
    {
      if (m_verbose) std::cerr << "Error: " << m_root->status().ToString() << std::endl;
      return;
    }

    if (m_verbose) std::cerr << " 4) Setting up loss calculation" << std::endl;
    
    // loss calculation based on cross-entropy
      
//      TFops::Log log_ypred = TFops::Log (*m_root, m_layers.back());

#ifdef USE_IOU

    TFops::Mul ygt_times_ypred = TFops::Mul (*m_root, *m_ph_gt, m_layers.back());
    TFops::ReduceSum intersection = TFops::ReduceSum(*m_root, ygt_times_ypred, {1});

    TFops::Add ygt_plus_ypred = TFops::Add (*m_root, *m_ph_gt, m_layers.back());
    TFops::Mul ygt_times_ypred_2 = TFops::Mul (*m_root, *m_ph_gt, m_layers.back());
    TFops::Sub sum_minus_product = TFops::Sub (*m_root, ygt_plus_ypred, ygt_times_ypred_2);
    TFops::ReduceSum my_union = TFops::ReduceSum(*m_root, sum_minus_product, {1});

    TFops::Div IoU = TFops::Div(*m_root, intersection, my_union);
    TFops::Sub one_minus_IoU = TFops::Sub (*m_root, TFops::Const(*m_root, 1.f), IoU);
    m_loss = new TFops::ReduceMean (*m_root, one_minus_IoU, {0});

    // TFops::ReduceMean mean_IoU = TFops::ReduceMean (*m_root, IoU, {0});
    // m_loss = new TFops::Sub (*m_root, TFops::Const(*m_root, 1.f), mean_IoU);
    
#else
    TFops::Maximum truncated_ypred = TFops::Maximum (*m_root, TFops::Const(*m_root, 0.0001f), m_layers.back());
    TFops::Log log_ypred = TFops::Log (*m_root, truncated_ypred);
    TFops::Mul ygt_times_log_ypred = TFops::Mul (*m_root, *m_ph_gt, log_ypred);
    TFops::ReduceSum sum_ygt_times_log_ypred = TFops::ReduceSum(*m_root, ygt_times_log_ypred, {1});
    TFops::Mul minus_sum_ygt_times_log_ypred = TFops::Mul (*m_root, TFops::Const(*m_root, -1.f), sum_ygt_times_log_ypred);

    m_loss = new TFops::ReduceMean (*m_root, minus_sum_ygt_times_log_ypred, {0});
#endif
    
    if (!m_root->status().ok())
    {
      if (m_verbose) std::cerr << "Error: " << m_root->status().ToString() << std::endl;
      return;
    }

    if (m_verbose) std::cerr << " 5) Setting up gradient descent" << std::endl;
    
    std::vector<TF::Output> weights_and_bias;
    for (std::size_t i = 0; i < weights.size(); ++ i)
      weights_and_bias.push_back (TF::Output(weights[i]));
    for (std::size_t i = 0; i < bias.size(); ++ i)
      weights_and_bias.push_back (TF::Output(bias[i]));
    
    std::vector<TF::Output> gradients;

    TF_CHECK_OK(TF::AddSymbolicGradients(*m_root, {*m_loss},
                                         weights_and_bias,
                                         &gradients));
  
    if (!m_root->status().ok())
    {
      if (m_verbose) std::cerr << "Error: " << m_root->status().ToString() << std::endl;
      return;
    }
  
    std::vector<TF::Output> assigners;
  
#ifdef USE_ADAM
    for (std::size_t i = 0; i < weights.size(); ++ i)
    {
      long long size_from = 0;
      long long size_to = 0;
      if (i == 0)
      {
        size_from = m_features.size();
        size_to = hl[0];
      }
      else if (i == weights.size() - 1)
      {
        size_from = hl.back();
        size_to = m_labels.size();
      }
      else
      {
        size_from = hl[i-1];
        size_to = hl[i];
      }

      TFops::Variable m (*m_root, { size_from, size_to}, TF::DT_FLOAT);
      TF::Tensor init_m
        (TF::DataTypeToEnum<float>::v(), 
         TF::TensorShape {(long long)(size_from), (long long)(size_to)});
      float* init_m_data = init_m.flat<float>().data();
      for (std::size_t s = 0; s < size_from; ++ s)
        for (std::size_t ss = 0; ss < size_to; ++ ss)
          init_m_data[ss * size_from + s] = 0.f;
      assigners.push_back (TFops::Assign (*m_root, m, init_m));

      TFops::Variable v (*m_root, { size_from, size_to}, TF::DT_FLOAT);
      TF::Tensor init_v
        (TF::DataTypeToEnum<float>::v(), 
         TF::TensorShape {(long long)(size_from), (long long)(size_to)});
      float* init_v_data = init_v.flat<float>().data();
      for (std::size_t s = 0; s < size_from; ++ s)
        for (std::size_t ss = 0; ss < size_to; ++ ss)
          init_v_data[ss * size_from + s] = 0.f;
      assigners.push_back (TFops::Assign (*m_root, v, init_v));

      m_gradient_descent.push_back (TFops::ApplyAdam
                                    (*m_root,
                                     weights[i],
                                     m,
                                     v,
                                     TFops::Cast(*m_root, 0.9, TF::DT_FLOAT),
                                     TFops::Cast(*m_root, 0.999, TF::DT_FLOAT),
                                     TFops::Cast(*m_root, m_learning_rate, TF::DT_FLOAT),
                                     TFops::Cast(*m_root, 0.9, TF::DT_FLOAT),
                                     TFops::Cast(*m_root, 0.999, TF::DT_FLOAT),
                                     TFops::Cast(*m_root, 1e-8, TF::DT_FLOAT),
                                     {gradients[i]}));
    }
    
    for (std::size_t i = 0; i < bias.size(); ++ i)
    {
      long long size_from = 0;
      long long size_to = 0;
      if (i == 0)
        size_to = hl[0];
      else if (i == bias.size() - 1)
        size_to = m_labels.size();
      else
        size_to = hl[i];

      TFops::Variable m (*m_root, { 1, size_to}, TF::DT_FLOAT);
      TF::Tensor init_m
        (TF::DataTypeToEnum<float>::v(), 
         TF::TensorShape {(long long)(1), (long long)(size_to)});
      float* init_m_data = init_m.flat<float>().data();
      for (std::size_t s = 0; s < size_to; ++ s)
        init_m_data[s] = 0.f;
      assigners.push_back (TFops::Assign (*m_root, m, init_m));

      TFops::Variable v (*m_root, { 1, size_to}, TF::DT_FLOAT);
      TF::Tensor init_v
        (TF::DataTypeToEnum<float>::v(), 
         TF::TensorShape {(long long)(1), (long long)(size_to)});
      float* init_v_data = init_v.flat<float>().data();
      for (std::size_t s = 0; s < size_to; ++ s)
        init_v_data[s] = 0.f;
      assigners.push_back (TFops::Assign (*m_root, v, init_v));

      m_gradient_descent.push_back (TFops::ApplyAdam
                                    (*m_root,
                                     bias[i],
                                     m,
                                     v,
                                     TFops::Cast(*m_root, 0.9, TF::DT_FLOAT),
                                     TFops::Cast(*m_root, 0.999, TF::DT_FLOAT),
                                     TFops::Cast(*m_root, m_learning_rate, TF::DT_FLOAT),
                                     TFops::Cast(*m_root, 0.9, TF::DT_FLOAT),
                                     TFops::Cast(*m_root, 0.999, TF::DT_FLOAT),
                                     TFops::Cast(*m_root, 1e-8, TF::DT_FLOAT),
                                     {gradients[weights.size() + i]}));
    }
    
#else
    for (std::size_t i = 0; i < weights.size(); ++ i)
      m_gradient_descent.push_back (TFops::ApplyGradientDescent
                                    (*m_root, weights[i],
                                     TFops::Cast(*m_root, m_learning_rate, TF::DT_FLOAT), {gradients[i]}));
    for (std::size_t i = 0; i < bias.size(); ++ i)
      m_gradient_descent.push_back (TFops::ApplyGradientDescent
                                    (*m_root, bias[i],
                                     TFops::Cast(*m_root, m_learning_rate, TF::DT_FLOAT), {gradients[weights.size() + i]}));
#endif
    if (!m_root->status().ok())
    {
      if (m_verbose) std::cerr << "Error: " << m_root->status().ToString() << std::endl;
      return;
    }

    if (m_verbose) std::cerr << " 6) Starting session" << std::endl;
    
    TF::SessionOptions options = TF::SessionOptions();
    
    // options.config.mutable_gpu_options()->set_visible_device_list("0");
    // options.config.mutable_gpu_options()->set_per_process_gpu_memory_fraction(0.3);
    // options.config.mutable_gpu_options()->set_operation_timeout_in_ms(15000);

//    options.config.mutable_gpu_options()->set_visible_device_list("");
    options.config.mutable_gpu_options()->set_allow_growth(true);
    options.config.mutable_gpu_options()->set_per_process_gpu_memory_fraction(0.8);
    options.config.mutable_gpu_options()->set_force_gpu_compatible(true);

    m_session = new TF::ClientSession (*m_root, options);
    if (!m_root->status().ok())
    {
      if (m_verbose) std::cerr << "Error: " << m_root->status().ToString() << std::endl;
      return;
    }

    std::vector<TF::Tensor> outputs;

    for (std::size_t i = 0; i < assign_weights.size(); ++ i)
      assigners.push_back (TF::Output(assign_weights[i]));
    for (std::size_t i = 0; i < assign_bias.size(); ++ i)
      assigners.push_back (TF::Output(assign_bias[i]));
    TF_CHECK_OK (m_session->Run(assigners, nullptr));

    if (!m_root->status().ok())
    {
      if (m_verbose) std::cerr << "Error: " << m_root->status().ToString() << std::endl;
      return;
    }
  }

};


#if 1
// Specialization to use GPU parallelization
template <typename ConcurrencyTag,
          typename ItemRange,
          typename LabelIndexRange,
          typename ProbabilitiesRanges>
void classify (const ItemRange& input,
               const Label_set& labels,
               const TensorFlow_neural_network_classifier& classifier,
               LabelIndexRange& output,
               ProbabilitiesRanges& probabilities)
{
  std::cerr << "Classify with TensorFlow classifier" << std::endl;
  
  output.resize(input.size());
  probabilities.resize (labels.size());
  for (std::size_t i = 0; i < probabilities.size(); ++ i)
    probabilities[i].resize (input.size());

  const std::size_t mem_allocated = sizeof(float) * input.size() * (labels.size() + classifier.features().size());
  const std::size_t size_max = 1024 * 1024 * 1024;
  const std::size_t nb_subdivisions = (mem_allocated / size_max) + 1;
  std::cerr << nb_subdivisions << " subdivision(s) for GPU processing" << std::endl;

  std::size_t idx = 0;
  for (std::size_t n = 0; n < nb_subdivisions; ++ n)
  {
    std::vector<std::size_t> indices;
    indices.reserve (input.size() / nb_subdivisions);
    for (std::size_t i = 0; i < input.size() / nb_subdivisions && idx < input.size(); ++ i)
      indices.push_back(idx ++);
    
    std::vector<std::vector<float> > values;
    classifier (indices, values);
    for(std::size_t i = 0; i < indices.size(); ++ i)
    {
      std::size_t nb_class_best = 0;
      float val_class_best = 0.f;

      for (std::size_t j = 0; j < labels.size(); ++ j)
      {
        probabilities[j][indices[i]] = values[i][j];
        
        if(val_class_best < values[i][j])
        {
          val_class_best = values[i][j];
          nb_class_best = j;
        }
      }

      output[indices[i]] = nb_class_best;
    }
  }
}

// Specialization to use GPU parallelization
template <typename ConcurrencyTag,
          typename ItemRange,
          typename LabelIndexRange>
void classify (const ItemRange& input,
               const Label_set& labels,
               const TensorFlow_neural_network_classifier& classifier,
               LabelIndexRange& output)
{
  std::cerr << "Classify with TensorFlow classifier" << std::endl;
  
  output.resize(input.size());

  const std::size_t mem_allocated = sizeof(float) * input.size() * (labels.size() + classifier.features().size());
  const std::size_t size_max = 1024 * 1024 * 1024;
  const std::size_t nb_subdivisions = (mem_allocated / size_max) + 1;
  std::cerr << nb_subdivisions << " subdivision(s) for GPU processing" << std::endl;

  std::size_t idx = 0;
  for (std::size_t n = 0; n < nb_subdivisions; ++ n)
  {
    std::vector<std::size_t> indices;
    indices.reserve (input.size() / nb_subdivisions);
    for (std::size_t i = 0; i < input.size() / nb_subdivisions && idx < input.size(); ++ i)
      indices.push_back(idx ++);
    
    std::vector<std::vector<float> > values;
    classifier (indices, values);

    for(std::size_t i = 0; i < indices.size(); ++ i)
    {
      std::size_t nb_class_best = 0;
      float val_class_best = 0.f;

      for (std::size_t j = 0; j < labels.size(); ++ j)
      {
        if(val_class_best < values[i][j])
        {
          val_class_best = values[i][j];
          nb_class_best = j;
        }
      }

      output[indices[i]] = nb_class_best;
    }
  }
}
#endif  
}

}

#endif // CGAL_CLASSIFICATION_TENSORFLOW_NEURAL_NETWORK_CLASSIFIER_H
