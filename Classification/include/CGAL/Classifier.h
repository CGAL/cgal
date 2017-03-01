// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Simon Giraudot, Florent Lafarge

#ifndef CGAL_CLASSIFIER_H
#define CGAL_CLASSIFIER_H

#include <cstdio>
#include <cassert>
#include <vector>
#include <list>
#include <set>
#include <string>
#include <queue>

#include <CGAL/bounding_box.h>
#include <CGAL/centroid.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <CGAL/Classification/Planimetric_grid.h>
#include <CGAL/Classification/Attribute_base.h>
#include <CGAL/Classification/Type.h>

#include <CGAL/Memory_sizer.h>
#include <CGAL/Timer.h>

#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>

//#define CGAL_CLASSIFICATION_VERBOSE
#if defined(CGAL_CLASSIFICATION_VERBOSE)
#define CGAL_CLASSIFICATION_CERR std::cerr
#else
#define CGAL_CLASSIFICATION_CERR std::ostream(0)
#endif

namespace CGAL {

/*!
\ingroup PkgClassification

\brief Classifies a data set based on a set of attributes and a set of
classification types.

This class implements the core of the classification algorithm
\cgalCite{cgal:lm-clscm-12}. It uses a data set as input and assigns
each input item to a classification type among a set of user defined
classification types. To achieve this classification, a set of local
geometric attributes are used, such as planarity, elevation or
vertical dispersion. In addition, the user must define a set of
classification types such as building, ground or vegetation.

Each pair of attribute and type must be assigned an
[Attribute::Effect](@ref CGAL::Classification::Attribute::Effect) (for
example, vegetation has a low planarity and a high vertical
dispersion) and each attribute must be assigned a weight. These
parameters can be set up by hand or by automatic training, provided a
small user-defined set of inlier is given for each classification
type.

\tparam ItemRange model of `ConstRange`. Its iterator type is
`RandomAccessIterator`.

\tparam ItemMap model of `ReadablePropertyMap` whose key
type is the value type of the iterator of `ItemRange` and value type is
the type of the items that are classified.
*/
template <typename ItemRange,
          typename ItemMap>
class Classifier
{

  
public:
  typedef typename Classification::Type_handle      Type_handle;
  typedef typename Classification::Attribute_handle Attribute_handle;
  
  /// \cond SKIP_IN_MANUAL
  typedef typename ItemMap::value_type Item;

#ifdef CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
  typedef internal::Alpha_expansion_graph_cut_boost             Alpha_expansion;
#else
  typedef internal::Alpha_expansion_graph_cut_boykov_kolmogorov Alpha_expansion;
#endif
  
protected:
  
  const ItemRange& m_input;
  ItemMap m_item_map;

  std::vector<std::size_t> m_assigned_type;
  std::vector<double> m_confidence;

  std::vector<Type_handle> m_types; 
  std::vector<Attribute_handle> m_attributes; 

  typedef Classification::Attribute::Effect Attribute_effect;
  std::vector<std::vector<Attribute_effect> > m_effect_table;

  /// \endcond

public:


  /// \name Constructor
  /// @{
  
  /*! 
    \brief Initializes a classification object.

    \param input input range.

    \param item_map property map to access the input items.
  */
  Classifier (const ItemRange& input,
              ItemMap item_map)
    : m_input (input), m_item_map (item_map)
  {
  }

  /// @}

  /// \cond SKIP_IN_MANUAL
  virtual ~Classifier() { }
  /// \endcond
  
  /// \name Attributes
  /// @{

  /*!
    \brief Adds an attribute.

    \tparam Attribute type of the attribute, inherited from
    `Classification::Attribute_base`.

    \tparam T types of the parameters of the attribute's constructor.

    \param t parameters of the attribute's constructor.

    \return a handle to the newly added attribute.
   */
#if (!defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE)) || DOXYGEN_RUNNING
  template <typename Attribute, typename ... T>
  Attribute_handle add_attribute (T&& ... t)
  {
    // CGAL::Timer timer;
    // timer.start();
    m_attributes.push_back (Attribute_handle (new Attribute(m_input, std::forward<T>(t)...)));
    // timer.stop();
    // CGAL_CLASSIFICATION_CERR << m_attributes.back()->name() << " took " << timer.time() << " second(s)" << std::endl;
    return m_attributes.back();
  }
#else
  template <typename Attribute>
  Attribute_handle add_attribute ()
  {
    m_attributes.push_back (Attribute_handle (new Attribute(m_input)));
    return m_attributes.back();
  }
  template <typename Attribute, typename T1>
  Attribute_handle add_attribute (T1& t1)
  {
    m_attributes.push_back (Attribute_handle (new Attribute(m_input, t1)));
    return m_attributes.back();
  }
  template <typename Attribute, typename T1, typename T2>
  Attribute_handle add_attribute (T1& t1, T2& t2)
  {
    m_attributes.push_back (Attribute_handle (new Attribute(m_input, t1, t2)));
    return m_attributes.back();
  }
  template <typename Attribute, typename T1, typename T2, typename T3>
  Attribute_handle add_attribute (T1& t1, T2& t2, T3& t3)
  {
    m_attributes.push_back (Attribute_handle (new Attribute(m_input, t1, t2, t3)));
    return m_attributes.back();
  }
  template <typename Attribute, typename T1, typename T2, typename T3, typename T4>
  Attribute_handle add_attribute (T1& t1, T2& t2, T3& t3, T4& t4)
  {
    m_attributes.push_back (Attribute_handle (new Attribute(m_input, t1, t2, t3, t4)));
    return m_attributes.back();
  }
  template <typename Attribute, typename T1, typename T2, typename T3, typename T4, typename T5>
  Attribute_handle add_attribute (T1& t1, T2& t2, T3& t3, T4& t4, T5& t5)
  {
    m_attributes.push_back (Attribute_handle (new Attribute(m_input, t1, t2, t3, t4, t5)));
    return m_attributes.back();
  }
#endif

  /*!
    \brief Removes an attribute.

    \param attribute the handle to attribute type that must be removed.

    \return `true` if the attribute was correctly removed, `false` if
    its handle was not found.
   */ 
  bool remove_attribute (Attribute_handle attribute)
  {
    for (std::size_t i = 0; i < m_attributes.size(); ++ i)
      if (m_attributes[i] == attribute)
        {
          m_attributes.erase (m_attributes.begin() + i);
          return true;
        }
    return false;
  }

  /*!
    \brief Returns how many attributes are defined.
  */  
  std::size_t number_of_attributes() const
  {
    return m_attributes.size();
  }


  /*!
    \brief Returns the i^{th} attribute.
  */  
  Attribute_handle attribute(std::size_t i)
  {
    return m_attributes[i];
  }

  /*!
    \brief Removes all attributes.
   */
  void clear_attributes ()
  {
    m_attributes.clear();
  }

  /// \endcond
  
  /// @}


  /// \name Classification Types
  /// @{
  
  /*!
    \brief Adds a classification type.

    \param name name of the classification type.

    \return a handle to the newly added classification type.
   */
  Type_handle add_classification_type (const char* name)
  {
    Type_handle out (new Classification::Type (name));
    m_types.push_back (out);
    return out;
  }

  /*!
    \brief Removes a classification type.

    \param type the handle to the classification type that must be removed.

    \return `true` if the classification type was correctly removed,
    `false` if its handle was not found.
   */ 
 bool remove_classification_type (Type_handle type)
  {
    std::size_t idx = (std::size_t)(-1);
    for (std::size_t i = 0; i < m_types.size(); ++ i)
      if (m_types[i] == type)
        {
          m_types.erase (m_types.begin() + i);
          idx = i;
          break;
        }
    if (idx == (std::size_t)(-1))
      return false;
    std::cerr << idx << std::endl;
    
    for (std::size_t i = 0; i < m_assigned_type.size(); ++ i)
      if (m_assigned_type[i] == (std::size_t)(-1))
          continue;
      else if (m_assigned_type[i] > idx)
        m_assigned_type[i] --;
      else if (m_assigned_type[i] == idx)
        m_assigned_type[i] = (std::size_t)(-1);

    return true;
  }

  /*!
    \brief Returns how many classification types are defined.
  */  
  std::size_t number_of_classification_types () const
  {
    return m_types.size();
  }
  
  /*!
    \brief Returns the i^{th} classification type.
  */  
  Type_handle classification_type (std::size_t i) const
  {
    return m_types[i];
  }


  /*!
    \brief Removes all classification types.
   */
  void clear_classification_types ()
  {
    m_types.clear();
  }

  /// @}

  /// \name Classification
  /// @{

  
  /*! 
    \brief Runs the classification algorithm without any regularization.

    There is no relationship between items, the classification energy
    is only minimized itemwise. This method is quick but produce
    suboptimal results.
  */
  void run()
  {
    prepare_classification ();
    
    // data term initialisation

    for (std::size_t s = 0; s < m_input.size(); s++)
      {

        std::size_t nb_class_best=0; 

        double val_class_best = (std::numeric_limits<double>::max)();
        std::vector<double> values;
      
        for(std::size_t k = 0; k < m_effect_table.size(); ++ k)
          {
            double value = classification_value (k, s);
            values.push_back (value);
          
            if(val_class_best > value)
              {
                val_class_best = value;
                nb_class_best=k;
              }
          }

        m_assigned_type[s] = nb_class_best;

        std::sort (values.begin(), values.end());
        m_confidence[s] = values[1] - values[0];
      }
  }


  /*! 
    \brief Runs the classification algorithm with a local smoothing.

    The computed classification energy is smoothed on a user defined
    local neighborhood of items. This method is a compromise between
    efficiency and reliability.

    \tparam NeighborQuery model of `NeighborQuery`.
    \param neighbor_query used to access neighborhoods of items.
  */
  template <typename NeighborQuery>
  void run_with_local_smoothing (const NeighborQuery& neighbor_query)
  {
    prepare_classification ();
    
    // data term initialisation
    CGAL_CLASSIFICATION_CERR << "Labeling... ";

    std::vector<std::vector<double> > values
      (m_types.size(),
       std::vector<double> (m_input.size(), -1.));

    for (std::size_t s=0; s < m_input.size(); ++ s)
      {
        std::vector<std::size_t> neighbors;
        neighbor_query (get (m_item_map, *(m_input.begin()+s)), std::back_inserter (neighbors));

        std::vector<double> mean (values.size(), 0.);
        for (std::size_t n = 0; n < neighbors.size(); ++ n)
          {
            if (values[0][neighbors[n]] < 0.)
              for(std::size_t k = 0; k < m_effect_table.size(); ++ k)
                {
                  values[k][neighbors[n]] = classification_value (k, neighbors[n]);
                  mean[k] += values[k][neighbors[n]];
                }
            else
              for (std::size_t j = 0; j < values.size(); ++ j)
                mean[j] += values[j][neighbors[n]];
          }

        std::size_t nb_class_best=0; 
        double val_class_best = (std::numeric_limits<double>::max)();
        for(std::size_t k = 0; k < mean.size(); ++ k)
          {
            mean[k] /= neighbors.size();
            if(val_class_best > mean[k])
              {
                val_class_best = mean[k];
                nb_class_best = k;
              }
          }

        m_assigned_type[s] = nb_class_best;

        std::sort (mean.begin(), mean.end());
        m_confidence[s] = mean[1] - mean[0];      
      }

  }


  /*! 
    \brief Runs the classification algorithm with a global
    regularization based on a graphcut.

    The computed classification energy is globally regularized through
    an alpha-expansion algorithm. This method is slow but provides
    the user with good quality results.

    \tparam NeighborQuery model of `NeighborQuery`.
    \param neighbor_query used to access neighborhoods of items.
    \param weight weight of the regularization with respect to the
    classification energy. Higher values produce more regularized
    output but may result in a loss of details.

  */
  template <typename NeighborQuery>
  void run_with_graphcut (const NeighborQuery& neighbor_query,
                          const double weight)
  {
    prepare_classification ();
    
    // data term initialisation
#ifdef CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
    CGAL_CLASSIFICATION_CERR << "Labeling using Boost with regularization weight " << weight << "... ";
#else
    CGAL_CLASSIFICATION_CERR << "Labeling using Boyvok Kolmogorov with regularization weight " << weight << "... ";
#endif

    std::vector<std::pair<std::size_t, std::size_t> > edges;
    std::vector<double> edge_weights;
    std::vector<std::vector<double> > probability_matrix
      (m_effect_table.size(), std::vector<double>(m_input.size(), 0.));
    std::vector<std::size_t>(m_input.size()).swap(m_assigned_type);

    std::cerr << "Size of probability matrix = " << m_effect_table.size() * m_input.size() << std::endl;
    std::cerr << "Size of assigned types = " << m_assigned_type.size() << std::endl;
    
    for (std::size_t s = 0; s < m_input.size(); ++ s)
      {
        std::vector<std::size_t> neighbors;

        neighbor_query (get(m_item_map, *(m_input.begin()+s)), std::back_inserter (neighbors));

        for (std::size_t i = 0; i < neighbors.size(); ++ i)
          if (s != neighbors[i])
            {
              edges.push_back (std::make_pair (s, neighbors[i]));
              edge_weights.push_back (weight);
            }
        
        std::size_t nb_class_best = 0;
        double val_class_best = (std::numeric_limits<double>::max)();
        for(std::size_t k = 0; k < m_effect_table.size(); ++ k)
          {
            double value = classification_value (k, s);
            probability_matrix[k][s] = value;
            
            if(val_class_best > value)
              {
                val_class_best = value;
                nb_class_best = k;
              }
          }
        m_assigned_type[s] = nb_class_best;
      }
    
    Alpha_expansion graphcut;
    graphcut(edges, edge_weights, probability_matrix, m_assigned_type);
    std::cerr << ((double)(CGAL::Memory_sizer().virtual_size()) / 1073741824.) << " GB allocated" << std::endl;
  }
  
  /// @}
  


  /// \name Output
  /// @{

  /*!
    \brief Returns the classification type of the item at position
    `index`.

    \note If classification was not performed (using `run()`,
    `run_with_local_smoothing()` or `run_with_graphcut()`), this
    function always returns the default `Type_handle`.
  */
  Type_handle classification_type_of (std::size_t index) const
  {
    if (m_assigned_type.size() <= index
        || m_assigned_type[index] == (std::size_t)(-1))
      {
      return Type_handle();
      }
    return m_types[m_assigned_type[index]];
  }

  /// \cond SKIP_IN_MANUAL
  bool classification_prepared() const
  {
    return !(m_assigned_type.empty());
  }
  void set_classification_type_of (std::size_t index, Type_handle class_type)
  {
    for (std::size_t i = 0; i < m_types.size(); ++ i)
      if (m_types[i] == class_type)
        {
          m_assigned_type[index] = i;
          return;
        }
    m_assigned_type[index] = (std::size_t)(-1);
  }
  /// \endcond

  /*!

    \brief Returns the confidence of the classification type of the
    item at position `index`.

    \note If classification was not performed (using `run()`,
    `run_with_local_smoothing()` or `run_with_graphcut()`), this
    function always returns 0.

    \return confidence ranging from 0 (not confident at all) to 1
    (very confident).
  */
  double confidence_of (std::size_t index) const
  {
    if (m_confidence.size() <= index)
      return 0.;
    return m_confidence[index];
  }

  /// @}


  double classification_value (Type_handle class_type, std::size_t pt_index)
  {
    double out = 0.;
    for (std::size_t i = 0; i < m_attributes.size(); ++ i)
      {
        if (m_attributes[i]->weight() == 0.)
          continue;

        Attribute_effect eff = class_type->attribute_effect (m_attributes[i]);

        if (eff == Classification::Attribute::FAVORING)
          out += m_attributes[i]->favored (pt_index);
        else if (eff == Classification::Attribute::PENALIZING)
          out += m_attributes[i]->penalized (pt_index);
        else if (eff == Classification::Attribute::NEUTRAL)
          out += m_attributes[i]->ignored (pt_index);
      }
    return out;
  }
  

protected:

  /// \cond SKIP_IN_MANUAL
  void prepare_classification ()
  {
    // Reset data structure
    std::vector<std::size_t>(m_input.size(), (std::size_t)(-1)).swap (m_assigned_type);
    std::vector<double>(m_input.size()).swap (m_confidence);

    m_effect_table = std::vector<std::vector<Attribute_effect> >
      (m_types.size(), std::vector<Attribute_effect> (m_attributes.size(),
                                                                 Classification::Attribute::NEUTRAL));
    
    for (std::size_t i = 0; i < m_effect_table.size (); ++ i)
      for (std::size_t j = 0; j < m_effect_table[i].size (); ++ j)
        m_effect_table[i][j] = m_types[i]->attribute_effect (m_attributes[j]);

  }
  
  double classification_value (const std::size_t& class_type,
                               const std::size_t& pt_index) const
  {
    double out = 0.;
    for (std::size_t i = 0; i < m_effect_table[class_type].size(); ++ i)
      {
        if (m_attributes[i]->weight() == 0.)
          continue;
        if (m_effect_table[class_type][i] == Classification::Attribute::FAVORING)
          out += m_attributes[i]->favored (pt_index);
        else if (m_effect_table[class_type][i] == Classification::Attribute::PENALIZING)
          out += m_attributes[i]->penalized (pt_index);
        else if (m_effect_table[class_type][i] == Classification::Attribute::NEUTRAL)
          out += m_attributes[i]->ignored (pt_index);
      }
    return out;
  }
  /// \endcond
};





} // namespace CGAL

#endif // CGAL_CLASSIFIER_H

