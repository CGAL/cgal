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

#ifndef CGAL_CLASSIFICATION_FEATURE_BASE_H
#define CGAL_CLASSIFICATION_FEATURE_BASE_H

#include <CGAL/license/Classification.h>

#include <memory>
#include <string>
#include <vector>
#include <stdexcept>

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
    \brief returns the name of the feature (initialized to
    `abstract_feature` for `Feature_base`).
  */
  const std::string& name() const { return m_name; }

  /*!
    \brief changes the name of the feature.
  */
  void set_name (const std::string& name) { m_name = name; }

  /*!
    \brief returns the value taken by the feature for at the item for
    the item at position `index`. This method must be implemented by
    inherited classes.
  */
  virtual float value (std::size_t index) = 0;

  /*!
    \brief returns the value taken by the multi-dimensional feature
    for at the item for the item at position `index` and in dimension
    `dim`. This method must be implemented by inherited classes for
    multi-dimensional features.

    \note The function `value(std::size_t index)` must still be implemented,
    but is not used. Simply defining it to return zero is sufficient.
  */
  virtual float value (std::size_t index, std::size_t dim) { return 0; }

};


#ifdef DOXYGEN_RUNNING
/*!
  \ingroup PkgClassificationFeature

  \brief %Handle to a `Feature_base`.

  \cgalModels{Handle}
*/
class Feature_handle { };
#else

class Feature_set;

class Feature_handle
{
  friend Feature_set;

  using Feature_base_ptr = std::unique_ptr<Feature_base>;
  std::shared_ptr<Feature_base_ptr> m_base;

  template <typename Feature_ptr>
  Feature_handle (Feature_ptr f)
    : m_base (std::make_shared<Feature_base_ptr>(std::move(f)))
  {
  }

  template <typename Feature_ptr>
  void attach (Feature_ptr f)
  {
    *m_base = std::move(f);
  }
public:

  Feature_handle() : m_base (std::make_shared<Feature_base_ptr>()) { }

  Feature_base& operator*() { return **m_base; }

  Feature_base* operator->() { return m_base->get(); }

  const Feature_base& operator*() const { return **m_base; }
  const Feature_base* operator->() const { return m_base->get(); }

  bool operator< (const Feature_handle& other) const { return *m_base < *(other.m_base); }
  bool operator== (const Feature_handle& other) const { return *m_base == *(other.m_base); }
};

#endif

// Use a throw to prevent the user from accidentally calling an irrelevant member function
class Missing_argument_exception : public std::exception { };

/*!
  \ingroup PkgClassificationFeature

  \brief %An extension of a `Feature_base` to support multi-dimensional features.

  \cgalModels{Feature_base}
*/

class Feature_base_multi_dim : public Feature_base
{
public:

  /// \cond SKIP_IN_MANUAL
  virtual float value(std::size_t index) final {
    throw Missing_argument_exception();
    return 0;
  }
  /// \endcond

  /*!
  \brief returns the value taken by the multi-dimensional feature
  for at the item for the item at position `index` and in dimension
  `dim`. This method must be implemented by inherited classes for
  multi-dimensional features.

  \note The function `value(std::size_t index)` from the base class
  `Feature_base` is not used for multi-dimensional features.
  */
  virtual float value(std::size_t index, std::size_t dim) = 0;
};

namespace Internal {
  class Feature_base_dim : public Feature_base
  {
    std::size_t m_dim;
    Feature_handle m_feature_handle;

  public:
    Feature_base_dim(std::size_t dim, std::string name, Feature_handle feature_handle) : m_dim(dim), m_feature_handle(feature_handle)
    {
      this->set_name(name + '[' + std::to_string(dim) + ']');
    }

    float value(std::size_t pt_index)
    {
      return m_feature_handle->value(pt_index, m_dim);
    }
  };
} // namespace Internal

/*!
  \ingroup PkgClassificationFeature

  \brief casts a feature handle to a specialized feature pointer.
*/
template <typename FeatureType>
FeatureType* feature_cast (Feature_handle fh)
{
  return dynamic_cast<FeatureType*>(&*(fh));
}


} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURE_BASE_H
