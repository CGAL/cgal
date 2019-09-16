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

#ifndef CGAL_CLASSIFICATION_FEATURE_BASE_H
#define CGAL_CLASSIFICATION_FEATURE_BASE_H

#include <CGAL/license/Classification.h>

#include <boost/shared_ptr.hpp>

#include <vector>

namespace CGAL {

namespace Classification {
  
/*!
  \ingroup PkgClassificationFeature

  \brief Abstract class describing a classification feature that
  associates a scalar value to each item of the classification input.
*/


class Feature_base
{
  std::string m_name;
  
public:

  /// \cond SKIP_IN_MANUAL
  Feature_base() : m_name ("abstract_feature") { }
  virtual ~Feature_base() { }
  /// \endcond

  /*!
    \brief Returns the name of the feature (initialized to
    `abstract_feature` for `Feature_base`).
  */
  const std::string& name() const { return m_name; }

  /*!
    \brief Changes the name of the feature.
  */
  void set_name (const std::string& name) { m_name = name; }
  
  /*!
    \brief Returns the value taken by the feature for at the item for
    the item at position `index`. This method must be implemented by
    inherited classes.
  */
  virtual float value (std::size_t index) = 0;

};


#ifdef DOXYGEN_RUNNING
/*!
  \ingroup PkgClassificationFeature

  \brief %Handle to a `Feature_base`.

  \cgalModels Handle
*/
class Feature_handle { };
#else
//typedef boost::shared_ptr<Feature_base> Feature_handle;

class Feature_set;
  
class Feature_handle
{
  friend Feature_set;
  
  boost::shared_ptr<boost::shared_ptr<Feature_base> > m_base;

  template <typename Feature>
  Feature_handle (Feature* f) : m_base (new boost::shared_ptr<Feature_base>(f)) { }

  template <typename Feature>
  void attach (Feature* f) const
  {
    *m_base = boost::shared_ptr<Feature_base>(f);
  }
public:

  Feature_handle() : m_base (new boost::shared_ptr<Feature_base>()) { }

  Feature_base& operator*() { return **m_base; }

  Feature_base* operator->() { return m_base->get(); }

  const Feature_base& operator*() const { return **m_base; }
  const Feature_base* operator->() const { return m_base->get(); }

  bool operator< (const Feature_handle& other) const { return *m_base < *(other.m_base); }
  bool operator== (const Feature_handle& other) const { return *m_base == *(other.m_base); }
};
  
#endif


} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURE_BASE_H
