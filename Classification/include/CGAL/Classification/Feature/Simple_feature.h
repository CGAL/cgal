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

#ifndef CGAL_CLASSIFICATION_SIMPLE_FEATURE_H
#define CGAL_CLASSIFICATION_SIMPLE_FEATURE_H

#include <CGAL/license/Classification.h>

#include <vector>
#include <CGAL/Classification/Feature_base.h>

namespace CGAL {

namespace Classification {

namespace Feature {

  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on a user-defined scalar field.

    \tparam InputRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam PropertyMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is statically castable to `float`.
  */
template <typename InputRange, typename PropertyMap>
class Simple_feature : public Feature_base
{
  const InputRange& m_input;
  PropertyMap m_pmap;

public:
  /*!
    \brief constructs the feature using an input range and a property map.

    \param input point range.
    \param property_map property map to access scalar field.
    \param name name of the property (no default name is given).
  */
  Simple_feature (const InputRange& input,
                  PropertyMap property_map,
                  const std::string& name)
    : m_input (input), m_pmap (property_map)
  {
    this->set_name (name);
  }
  /// \cond SKIP_IN_MANUAL
  virtual float value (std::size_t pt_index)
  {
    return static_cast<float> (get (m_pmap, *(m_input.begin()+pt_index)));
  }
  /// \endcond
};

} // namespace Feature

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_SIMPLE_FEATURE_H
