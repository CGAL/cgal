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

#ifndef CGAL_CLASSIFICATION_FEATURE_SET_H
#define CGAL_CLASSIFICATION_FEATURE_SET_H

#include <CGAL/license/Classification.h>

#include <CGAL/Classification/Feature_base.h>

#include <boost/make_shared.hpp>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/mutex.h>
#include <tbb/task_group.h>
#endif // CGAL_LINKED_WITH_TBB

#include <vector>
#include <utility>

namespace CGAL {

namespace Classification {

/*!
\ingroup PkgClassificationFeature

\brief Set of features (see `Feature_base`) used as input by
classification algorithms. This class handles both the instantiation,
the addition and the deletion of features.

*/
class Feature_set
{
  typedef std::vector<Feature_handle> Base;
  Base m_features;

  struct Compare_name
  {
    bool operator() (const Feature_handle& a, const Feature_handle& b) const
    {
      if (a->name() == b->name())
        return a < b;
      return a->name() < b->name();
    }
  };
  
#ifdef CGAL_LINKED_WITH_TBB
  tbb::task_group* m_tasks;
#endif // CGAL_LINKED_WITH_TBB
  
public:

  /// \name Constructor
  /// @{

  /*!
    \brief Creates an empty feature set.
  */
  Feature_set()
#ifdef CGAL_LINKED_WITH_TBB
    : m_tasks(NULL)
#endif
  { }

  /// @}
  
  /// \cond SKIP_IN_MANUAL
  virtual ~Feature_set()
  {
#ifdef CGAL_LINKED_WITH_TBB
    if (m_tasks != NULL)
      delete m_tasks;
    for (std::size_t i = 0; i < m_adders.size(); ++ i)
      delete m_adders[i];
#endif
  }
  /// \endcond

  /// \name Modifications
  /// @{
  
  /*!
    \brief Instantiates a new feature and adds it to the set.

    If several calls of `add()` are surrounded by
    `begin_parallel_additions()` and `end_parallel_additions()`, they
    are computed in parallel. They are still inserted in the specified
    order in the feature set (the first call of `add()` creates a
    feature at index 0, the second at index 1, etc.).

    \sa `begin_parallel_additions()`
    \sa `end_parallel_additions()`

    \tparam Feature type of the feature, inherited from
    `Feature_base`.

    \tparam T types of the parameters of the feature's constructor.

    \param t parameters of the feature's constructor.

    \return a handle to the newly added feature.
  */
  template <typename Feature, typename ... T>
  Feature_handle add (T&& ... t)
  {
#ifdef CGAL_LINKED_WITH_TBB
    if (m_tasks != NULL)
    {
      m_features.push_back (Feature_handle());
    
      Parallel_feature_adder<Feature, T...>* adder
        = new Parallel_feature_adder<Feature, T...>(m_features.back(), std::forward<T>(t)...);
      
      m_adders.push_back (adder);
      m_tasks->run (*adder);
    }
    else
#endif
    {
       m_features.push_back (Feature_handle (new Feature(std::forward<T>(t)...)));
    }
    return m_features.back();
  }

  /// \cond SKIP_IN_MANUAL
  template <typename Feature, typename ... T>
  Feature_handle add_with_scale_id (std::size_t i, T&& ... t)
  {
#ifdef CGAL_LINKED_WITH_TBB
    if (m_tasks != NULL)
    {
      m_features.push_back (Feature_handle());
    
      Parallel_feature_adder<Feature, T...>* adder
        = new Parallel_feature_adder<Feature, T...>(i, m_features.back(), std::forward<T>(t)...);
      
      m_adders.push_back (adder);
      m_tasks->run (*adder);
    }
    else
#endif
    {
       m_features.push_back (Feature_handle (new Feature(std::forward<T>(t)...)));
       m_features.back()->set_name (m_features.back()->name() + "_" + std::to_string(i));
    }
    return m_features.back();
  }
  /// \endcond

    
  /*!
    \brief Removes a feature.

    \param feature the handle to feature type that must be removed.

    \return `true` if the feature was correctly removed, `false` if
    its handle was not found.
  */ 
  bool remove (Feature_handle feature)
  {
    for (std::size_t i = 0; i < m_features.size(); ++ i)
      if (m_features[i] == feature)
      {
        m_features.erase (m_features.begin() + i);
        return true;
      }
    return false;
  }

  /*!
    \brief Removes all features.
  */
  void clear ()
  {
    m_features.clear();
  }

  /// @}

  /// \name Parallel Processing
  /// @{
  

#if defined(CGAL_LINKED_WITH_TBB) || defined(DOXYGEN_RUNNING)

  /*!
    \brief Initializes structures to compute features in parallel.

    If the user wants to add features in parallel, this function
    should be called before making several calls of `add()`. After the
    calls of `add()`, `end_parallel_additions()` should be called.

    \note This function requires \ref thirdpartyTBB.

    \warning As arguments of `add()` are passed by reference and that new
    threads are started if `begin_parallel_additions()` is used, it is
    highly recommended to always call `begin_parallel_additions()`,
    `add()` and `end_parallel_additions()` _within the same scope_, to
    avoid keeping references to temporary objects that might be
    deleted before the thread has terminated.

    \sa `end_parallel_additions()`
  */ 
  void begin_parallel_additions()
  {
    m_tasks = new tbb::task_group;
  }

  /*!

    \brief Waits for the end of parallel feature computation and
    clears dedicated data structures afterwards.

    If the user wants to add features in parallel, this function
    should be called after `begin_parallel_additions()` and several
    calls of `add()`.

    \note This function requires \ref thirdpartyTBB.

    \sa `begin_parallel_additions()`
  */ 
  void end_parallel_additions()
  {
    m_tasks->wait();
    delete m_tasks;
    m_tasks = NULL;
    
    for (std::size_t i = 0; i < m_adders.size(); ++ i)
      delete m_adders[i];
    m_adders.clear();
  }
#endif

  /// @}


  /// \name Access
  /// @{
  
  /*!
    \brief Returns how many features are defined.
  */  
  std::size_t size() const
  {
    return m_features.size();
  }


  /*!
    \brief Returns the \f$i^{th}\f$ feature.
  */  
  Feature_handle operator[](std::size_t i) const
  {
    return m_features[i];
  }

  /// @}

  /// \cond SKIP_IN_MANUAL
  void free_memory(std::size_t i)
  {
    m_features[i] = Feature_handle();
  }

  void sort_features_by_name()
  {
    std::sort (m_features.begin(), m_features.end(),
               Compare_name());               
  }
  /// \endcond

private:

  /// \cond SKIP_IN_MANUAL

  struct Abstract_parallel_feature_adder
  {
    virtual ~Abstract_parallel_feature_adder() { }
    virtual void operator()() const = 0;
  };
  
  template <typename Feature, typename ... T>
  struct Parallel_feature_adder : Abstract_parallel_feature_adder
  {
    std::size_t scale;
    mutable Feature_handle fh;
    boost::shared_ptr<std::tuple<T...> > args;
    
    Parallel_feature_adder (Feature_handle fh, T&& ... t)
      : scale (std::size_t(-1)), fh (fh)
    {
      args = boost::make_shared<std::tuple<T...> >(std::forward<T>(t)...);
    }
    
    Parallel_feature_adder (std::size_t scale, Feature_handle fh, T&& ... t)
      : scale(scale), fh (fh)
    {
      args = boost::make_shared<std::tuple<T...> >(std::forward<T>(t)...);
    }

    template<int ...>
    struct seq { };

    template<int N, int ...S>
    struct gens : gens<N-1, N-1, S...> { };

    template<int ...S>
    struct gens<0, S...> {
      typedef seq<S...> type;
    };

    template <typename Type>
    const Type& remove_ref_of_simple_type (const Type& t) const { return t; }

    
    template <typename Tuple, int ... S>
    void add_feature (Tuple& t, seq<S...>) const
    {
      fh.attach (new Feature (std::forward<T>(std::get<S>(t))...));
      if (scale != std::size_t(-1))
        fh->set_name (fh->name() + "_" + std::to_string(scale));
    }

    void operator()() const
    {
      add_feature(*args, typename gens<sizeof...(T)>::type());
    }

  };

  std::vector<Abstract_parallel_feature_adder*> m_adders;

  /// \endcond
};



} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURE_SET_H
