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

#ifndef CGAL_DATA_CLASSIFICATION_ATTRIBUTE_COLOR_H
#define CGAL_DATA_CLASSIFICATION_ATTRIBUTE_COLOR_H

#include <vector>

#include <CGAL/Data_classification/Color.h>

namespace CGAL {

namespace Data_classification {

  /*!
    \ingroup PkgDataClassification

    \brief Attribute based on HSV colorimetric information.

    If the input point cloud has colorimetric information, it can be
    used for classification purposes. This attribute is based on
    Gaussian probabilistic model on one of the three HSV channels
    (hue, saturation or value).

    It computes the probability of the color of the input point to
    match this specific color channel defined by a mean and a standard
    deviation.

    For example, such an attribute using the channel 0 (hue) with a
    mean of 90 (which corresponds to a green hue) can help identifying
    trees.

    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam ColorMap is a model of `ReadablePropertyMap` with value type `CGAL::Data_classification::RGB_Color`.A
  */
template <typename Kernel, typename RandomAccessIterator, typename ColorMap>
class Attribute_hsv : public Attribute
{
  typedef typename Data_classification::RGB_Color RGB_Color;
  typedef typename Data_classification::HSV_Color HSV_Color;
  
  std::vector<double> color_attribute;
  std::string m_id;
  
public:
  
  /*!

    \brief Constructs an attribute based on the given color channel,
    mean and standard deviation.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param color_map Property map to access the colors of the input points
    \param channel Chosen HSV channel (0 for hue, 1 for saturation, 2 for value).
    \param mean Mean value of the specified channel
    \param sd Standard deviation of the specified channel
  */
  Attribute_hsv (RandomAccessIterator begin,
                 RandomAccessIterator end,
                 ColorMap color_map,
                 std::size_t channel,
                 double mean, double sd)
  {
    this->weight = 1.;
    for(std::size_t i = 0; i < (std::size_t)(end - begin);i++)
      {
        HSV_Color c = Data_classification::rgb_to_hsv (get(color_map, begin[i]));
        color_attribute.push_back (std::exp (-(c[channel] - mean) * (c[channel] - mean) / (2. * sd * sd)));
      }
    this->compute_mean_max (color_attribute, this->mean, this->max);

    std::ostringstream oss;
    if (channel == 0) oss << "hue";
    else if (channel == 1) oss << "saturation";
    else if (channel == 2) oss << "value";
    oss << "_" << mean;
    m_id = oss.str();
  }

  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return color_attribute[pt_index];
  }

  virtual std::string id() { return m_id; }
  /// \endcond
};

} // namespace Data_classification

} // namespace CGAL

#endif // CGAL_DATA_CLASSIFICATION_ATTRIBUTE_COLOR_H
