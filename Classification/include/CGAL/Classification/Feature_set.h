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

#ifndef CGAL_CLASSIFICATION_FEATURE_SET_H
#define CGAL_CLASSIFICATION_FEATURE_SET_H

#include <CGAL/Classification/Feature_base.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/mutex.h>
#endif // CGAL_LINKED_WITH_TBB

#include <vector>

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
  
#ifdef CGAL_LINKED_WITH_TBB
  tbb::mutex m_mutex;
  void mutex_lock() { m_mutex.lock(); }
  void mutex_unlock() { m_mutex.unlock(); }
#else // CGAL_LINKED_WITH_TBB
  void mutex_lock() { }
  void mutex_unlock() { }
#endif // CGAL_LINKED_WITH_TBB
  
public:

  Feature_set() { }
  
  /// \cond SKIP_IN_MANUAL
  virtual ~Feature_set() { }
  /// \endcond

  /*!
    \brief Instantiates a new feature and adds it to the set.

    \tparam Feature type of the feature, inherited from
    `Feature_base`.

    \tparam T types of the parameters of the feature's constructor.

    \param t parameters of the feature's constructor.

    \return a handle to the newly added feature.
  */
#if (!defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE)) || DOXYGEN_RUNNING
  template <typename Feature, typename ... T>
  Feature_handle add (T&& ... t)
  {
    Feature_handle fh (new Feature(std::forward<T>(t)...));
    mutex_lock();
    m_features.push_back (fh);
    mutex_unlock();
    return fh;
  }
#else
  template <typename Feature, typename T1>
  Feature_handle add (T1& t1)
  {
    Feature_handle fh (new Feature(t1));
    mutex_lock();
    m_features.push_back (fh);
    mutex_unlock();
    return fh;
  }
  template <typename Feature, typename T1, typename T2>
  Feature_handle add (T1& t1, T2& t2)
  {
    Feature_handle fh (new Feature(t1, t2));
    mutex_lock();
    m_features.push_back (fh);
    mutex_unlock();
    return fh;
  }
  template <typename Feature, typename T1, typename T2, typename T3>
  Feature_handle add (T1& t1, T2& t2, T3& t3)
  {
    Feature_handle fh (new Feature(t1, t2, t3));
    mutex_lock();
    m_features.push_back (fh);
    mutex_unlock();
    return fh;
  }
  template <typename Feature, typename T1, typename T2, typename T3, typename T4>
  Feature_handle add (T1& t1, T2& t2, T3& t3, T4& t4)
  {
    Feature_handle fh (new Feature(t1, t2, t3, t4));
    mutex_lock();
    m_features.push_back (fh);
    mutex_unlock();
    return fh;
  }
  template <typename Feature, typename T1, typename T2, typename T3, typename T4, typename T5>
  Feature_handle add (T1& t1, T2& t2, T3& t3, T4& t4, T5& t5)
  {
    Feature_handle fh (new Feature(t1, t2, t3, t4, t5));
    mutex_lock();
    m_features.push_back (fh);
    mutex_unlock();
    return fh;
  }
#endif

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

  /*!
    \brief Removes all features.
  */
  void clear ()
  {
    m_features.clear();
  }

  /// \cond SKIP_IN_MANUAL
  void free_memory(std::size_t i)
  {
    m_features[i] = Feature_handle();
  }
  /// \endcond
  
};



} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURE_SET_H
