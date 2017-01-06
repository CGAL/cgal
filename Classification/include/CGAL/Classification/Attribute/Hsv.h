// Copyright (c) 2016  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_CLASSIFICATION_ATTRIBUTE_HSV_H
#define CGAL_CLASSIFICATION_ATTRIBUTE_HSV_H

#include <vector>

#include <CGAL/Classification/Color.h>

namespace CGAL {

namespace Classification {

namespace Attribute {
  
  /*!
    \ingroup PkgClassificationAttributes

    \brief Attribute based on HSV colorimetric information.

    If the input point cloud has colorimetric information, it can be
    used for classification purposes. This attribute is based on a
    Gaussian probabilistic model on one of the three HSV channels
    (hue, saturation or value).

    It computes the probability of the color of the input point to
    match this specific color channel defined by a mean and a standard
    deviation.

    For example, such an attribute using the channel 0 (hue) with a
    mean of 90 (which corresponds to a green hue) can help identify
    trees.

    \tparam Kernel model of \cgal Kernel.
    \tparam Range range of items, model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam ColorMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `Range` and value type
    is `CGAL::Classification::RGB_Color`.
  */
template <typename Kernel, typename Range, typename ColorMap>
class Hsv : public Attribute_base
{
  typedef typename Classification::RGB_Color RGB_Color;
  typedef typename Classification::HSV_Color HSV_Color;
  
  std::vector<double> color_attribute;
  std::string m_name;
  
public:
  
  /*!

    \brief Constructs an attribute based on the given color channel,
    mean and standard deviation.

    \param input input range.
    \param color_map property map to access the colors of the input points.
    \param channel chosen HSV channel (0 for hue, 1 for saturation, 2 for value).
    \param mean mean value of the specified channel.
    \param sd standard deviation of the specified channel.
  */
  Hsv (const Range& input,
       ColorMap color_map,
       std::size_t channel,
       double mean, double sd)
  {
    this->set_weight(1.);
    for(std::size_t i = 0; i < input.size();i++)
      {
        HSV_Color c = Classification::rgb_to_hsv (get(color_map, *(input.begin()+i)));
        color_attribute.push_back (std::exp (-(c[channel] - mean) * (c[channel] - mean) / (2. * sd * sd)));
      }
    this->compute_mean_max (color_attribute, this->mean, this->max);

    std::ostringstream oss;
    if (channel == 0) oss << "hue";
    else if (channel == 1) oss << "saturation";
    else if (channel == 2) oss << "value";
    oss << "_" << mean;
    m_name = oss.str();
  }

  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return color_attribute[pt_index];
  }

  virtual std::string name() { return m_name; }
  /// \endcond
};

} // namespace Attribute

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_ATTRIBUTE_HSV_H
