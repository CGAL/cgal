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
  typedef TFops::Tanh Activation_function;
//  typedef TFops::Relu Activation_function;
  
  const Label_set& m_labels;
  const Feature_set& m_features;
  std::vector<float> m_feature_means;
  std::vector<float> m_feature_sd;

  TF::Scope* m_root;
  TFops::Placeholder* m_ph_ft;
  std::vector<TF::Output> m_layers;
  TF::ClientSession* m_session;
  
public:
  
  /// \name Constructor
  /// @{
  
  /*!
    \brief Instantiate the classifier using the sets of `labels` and `features`.

  */
  TensorFlow_neural_network_classifier (const Label_set& labels,
                                        const Feature_set& features)
    : m_labels (labels), m_features (features)
    , m_root (NULL), m_ph_ft (NULL), m_session (NULL)
  { }
  
  /// \cond SKIP_IN_MANUAL
  ~TensorFlow_neural_network_classifier ()
  {
  }

  void compute_normalization_coefficients (const std::vector<std::size_t>& indices)
  {
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
    }
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
              std::size_t number_of_iterations = 5000,
              float learning_rate = 0.5,
              const std::vector<std::size_t>& hidden_layers
              = std::vector<std::size_t>())
  {
    m_root = new TF::Scope (TF::Scope::NewRootScope());

    // Get input
    std::vector<std::size_t> indices;
    std::vector<int> raw_gt;
    for (std::size_t i = 0; i < ground_truth.size(); ++ i)
    {
      int gc = int(ground_truth[i]);
      if (gc != -1)
      {
        indices.push_back (i);
        raw_gt.push_back (gc);
      }
    }

    CGAL_CLASSTRAINING_CERR << "I - BUILDING NEURAL NETWORK ARCHITECTURE" << std::endl;
    
    // Get layers and sizes or init with default values
    std::vector<std::size_t> hl = hidden_layers;
    if (hl.empty())
    {
      hl.push_back (m_features.size());
      hl.push_back ((m_features.size() + m_labels.size()) / 2);
    }
    
    CGAL_CLASSTRAINING_CERR << " 1) Initializing architecture:" << std::endl
                            << "   * Layer 0: " << m_features.size() << " neuron(s) (input features)" << std::endl;
      
    if (!CGAL_CLASSTRAINING_SILENT)
      for (std::size_t i = 0; i < hl.size(); ++ i)
        std::cerr << "   * Layer " << i+1 << ": " << hl[i] << " neuron(s)" << std::endl;
    
    CGAL_CLASSTRAINING_CERR << "   * Layer " << hl.size() + 1 << ": "
                            << m_labels.size() << " neuron(s) (output labels)" << std::endl;
    
    m_ph_ft = new TFops::Placeholder (*m_root, TF::DT_FLOAT);
    TFops::Placeholder ph_gt (*m_root, TF::DT_FLOAT);
    
    

    CGAL_CLASSTRAINING_CERR << " 2) Creating weight matrices and bias" << std::endl;
    
    // create weight matrices and bias for each layer transition
    std::vector<TFops::Variable> weights;
    std::vector<TFops::Assign> assign_weights;
    std::vector<TFops::Variable> bias;
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

      CGAL_CLASSTRAINING_CERR << "   * Weight matrix " << i << " [" << size_from << ";" << size_to << "]" << std::endl;
      weights.push_back (TFops::Variable (*m_root, { size_from, size_to }, TF::DT_FLOAT));
      assign_weights.push_back (TFops::Assign (*m_root, weights.back(),
                                               TFops::RandomNormal (*m_root, { size_from, size_to}, TF::DT_FLOAT)));

      CGAL_CLASSTRAINING_CERR << "   * Bias " << i << " [" << size_to << "]" << std::endl;
      bias.push_back (TFops::Variable (*m_root, { (long long)(1), size_to }, TF::DT_FLOAT));

      TF::Tensor init_bias
        (TF::DataTypeToEnum<float>::v(), 
         TF::TensorShape {(long long)(1), (long long)(size_to)});
      float* init_bias_data = init_bias.flat<float>().data();
      for (std::size_t s = 0; s < size_to; ++ s)
        init_bias_data[s] = 0.f;

      assign_bias.push_back (TFops::Assign (*m_root, bias.back(), init_bias));
    }

    if (!m_root->status().ok())
    {
      CGAL_CLASSTRAINING_CERR << "Error: " << m_root->status().ToString() << std::endl;
      return;
    }
    
    CGAL_CLASSTRAINING_CERR << " 3) Creating layers" << std::endl;
    
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
      CGAL_CLASSTRAINING_CERR << "Error: " << m_root->status().ToString() << std::endl;
      return;
    }

    CGAL_CLASSTRAINING_CERR << " 4) Setting up loss calculation" << std::endl;
    
    // loss calculation
    TFops::ReduceMean loss
      = TFops::ReduceMean (*m_root,
                           TFops::Mul (*m_root, TFops::Const(*m_root, -1.f),
                                       TFops::ReduceSum(*m_root,
                                                        TFops::Mul (*m_root, ph_gt,
                                                                    TFops::Log (*m_root, m_layers.back())), {1})),
                           {0});
    
    if (!m_root->status().ok())
    {
      CGAL_CLASSTRAINING_CERR << "Error: " << m_root->status().ToString() << std::endl;
      return;
    }

    CGAL_CLASSTRAINING_CERR << " 5) Setting up gradient descent" << std::endl;
    
    std::vector<TF::Output> weights_and_bias;
    for (std::size_t i = 0; i < weights.size(); ++ i)
      weights_and_bias.push_back (TF::Output(weights[i]));
    for (std::size_t i = 0; i < bias.size(); ++ i)
      weights_and_bias.push_back (TF::Output(bias[i]));
    
    std::vector<TF::Output> gradients;

    TF_CHECK_OK(TF::AddSymbolicGradients(*m_root, {loss},
                                         weights_and_bias,
                                         &gradients));

    if (!m_root->status().ok())
    {
      CGAL_CLASSTRAINING_CERR << "Error: " << m_root->status().ToString() << std::endl;
      return;
    }

    std::vector<TFops::ApplyGradientDescent> gradient_descent;
    for (std::size_t i = 0; i < weights.size(); ++ i)
      gradient_descent.push_back (TFops::ApplyGradientDescent
                                  (*m_root, weights[i],
                                   TFops::Cast(*m_root, learning_rate, TF::DT_FLOAT), {gradients[i]}));
    for (std::size_t i = 0; i < bias.size(); ++ i)
      gradient_descent.push_back (TFops::ApplyGradientDescent
                                  (*m_root, bias[i],
                                   TFops::Cast(*m_root, learning_rate, TF::DT_FLOAT), {gradients[weights.size() + i]}));
    
    if (!m_root->status().ok())
    {
      CGAL_CLASSTRAINING_CERR << "Error: " << m_root->status().ToString() << std::endl;
      return;
    }

    CGAL_CLASSTRAINING_CERR << " 6) Starting session" << std::endl;
    
    m_session = new TF::ClientSession (*m_root);
    std::vector<TF::Tensor> outputs;

    std::vector<TF::Output> assign_weights_and_bias;
    for (std::size_t i = 0; i < assign_weights.size(); ++ i)
      assign_weights_and_bias.push_back (TF::Output(assign_weights[i]));
    for (std::size_t i = 0; i < assign_bias.size(); ++ i)
      assign_weights_and_bias.push_back (TF::Output(assign_bias[i]));
    TF_CHECK_OK (m_session->Run(assign_weights_and_bias, nullptr));

    if (!m_root->status().ok())
    {
      CGAL_CLASSTRAINING_CERR << "Error: " << m_root->status().ToString() << std::endl;
      return;
    }
    
    CGAL_CLASSTRAINING_CERR << " 7. Normalizing features" << std::endl;
    
    compute_normalization_coefficients (indices);
    
    CGAL_CLASSTRAINING_CERR << " 8. Constructing feature tensor and ground truth tensor" << std::endl;

    TF::Tensor ft
      (TF::DataTypeToEnum<float>::v(), 
       TF::TensorShape {(long long)(raw_gt.size()), (long long)(m_features.size())});
    TF::Tensor gt
      (TF::DataTypeToEnum<float>::v(), 
       TF::TensorShape {(long long)(raw_gt.size()), (long long)(m_labels.size())});

    float* ft_data = ft.flat<float>().data();
    float* gt_data = gt.flat<float>().data();
    
    // Fill input tensors
    for (std::size_t i = 0; i < indices.size(); ++ i)
    {
      std::size_t idx = indices[i];
      int g = raw_gt[i];

      for (std::size_t f = 0; f < m_features.size(); ++ f)
        ft_data[i * m_features.size() + f]
          = (m_features[f]->value(idx) - m_feature_means[f]) / m_feature_sd[f];

      for (std::size_t l = 0; l < m_labels.size(); ++ l)
        if (std::size_t(g) == l)
          gt_data[i * m_labels.size() + l] = 1.f;
        else
          gt_data[i * m_labels.size() + l] = 0.f;
    }

    CGAL_CLASSTRAINING_CERR << "II - TRAINING NEURAL NETWORK" << std::endl;
    
    std::vector<TF::Output> operations;
    for (std::size_t i = 0; i < gradient_descent.size(); ++ i)
      operations.push_back (gradient_descent[i]);
    operations.push_back (m_layers.back());
    
    for (std::size_t i = 0; i < number_of_iterations; ++ i)
    {
      TF_CHECK_OK (m_session->Run ({{*m_ph_ft, ft}, {ph_gt, gt}}, {loss}, &outputs));

      if ((i+1) % (number_of_iterations / 20) == 0)
      {
        CGAL_CLASSTRAINING_CERR << "   * Step " << i+1 << "/" << number_of_iterations << ": loss = "
                                << outputs[0].scalar<float>() << std::endl;
      }
      if (!std::isfinite(*outputs[0].scalar<float>().data()))
      {
        std::cerr << "Loss is " << outputs[0].scalar<float>() << ", aborting" << std::endl;
        return;
      }

      TF_CHECK_OK(m_session->Run({{*m_ph_ft, ft}, {ph_gt, gt}}, operations, nullptr));
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
  void save_configuration (std::ostream& )
  {
  }

  /*!
    \brief Loads a configuration from the stream `input`.

    The input file should be a GZIP container written by the
    `save_configuration()` method. The feature set of the classifier
    should contain the exact same features in the exact same order as
    the ones present when the file was generated using
    `save_configuration()`.
  */
  void load_configuration (std::istream& )
  {
  }

};

}

}

#endif // CGAL_CLASSIFICATION_TENSORFLOW_NEURAL_NETWORK_CLASSIFIER_H
