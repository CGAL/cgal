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

#ifndef CGAL_CLASSIFICATION_FEATURE_HSV_H
#define CGAL_CLASSIFICATION_FEATURE_HSV_H

#include <vector>

#include <CGAL/Classification/Color.h>

namespace CGAL {

namespace Classification {

namespace Feature {
  
  /*!
    \ingroup PkgClassificationFeatures

    %Feature based on HSV colorimetric information. If the input
    point cloud has colorimetric information, it can be used for
    classification purposes. This feature is based on a Gaussian
    probabilistic model on one of the three HSV channels (hue,
    saturation or value). It computes the probability of the color of
    the input point to match this specific color channel defined by a
    mean and a standard deviation.

    The HSV channels are defined this way:

    - Hue ranges from 0 to 360 and measures the general "tint" of the
      color (green, blue, pink, etc.)

    - Saturation ranges from 0 to 100 and measures the "strength" of the
      color (0 is gray and 100 is the fully saturated color)

    - Value ranges from 0 to 100 and measures the "brightness" of the
      color (0 is black and 100 is the fully bright color)

    For example, such an feature using the channel 0 (hue) with a
    mean of 90 (which corresponds to a green hue) can help to identify
    trees.
    
    \note The user only needs to provide a map to standard (and more common)
    RGB colors, the conversion to HSV is done internally.

    \tparam Geom_traits model of \cgal Kernel.
    \tparam PointRange model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam ColorMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `CGAL::Classification::RGB_Color`.
  */
template <typename Geom_traits, typename PointRange, typename ColorMap>
class Hsv : public Feature_base
{
  typedef typename Classification::RGB_Color RGB_Color;
  typedef typename Classification::HSV_Color HSV_Color;

#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
  std::vector<double> color_feature;
#else
  const PointRange& input;
  ColorMap color_map;
  std::size_t m_channel;
  double m_mean;
  double m_sd;
#endif  
  
public:
  
  /*!

    \brief Constructs an feature based on the given color channel,
    mean and standard deviation.

    \param input input range.
    \param color_map property map to access the colors of the input points.
    \param channel chosen HSV channel (0 for hue, 1 for saturation, 2 for value).
    \param mean mean value of the specified channel.
    \param sd standard deviation of the specified channel.
  */
  Hsv (const PointRange& input,
       ColorMap color_map,
       std::size_t channel,
       double mean, double sd)
#ifndef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
    : input(input), color_map(color_map), m_channel(channel), m_mean(mean), m_sd(sd)
#endif
  {

#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
    for(std::size_t i = 0; i < input.size();i++)
      {
        HSV_Color c = Classification::rgb_to_hsv (get(color_map, *(input.begin()+i)));
        color_feature.push_back (std::exp (-(c[channel] - mean) * (c[channel] - mean) / (2. * sd * sd)));
      }
#endif

    std::ostringstream oss;
    if (channel == 0) oss << "hue";
    else if (channel == 1) oss << "saturation";
    else if (channel == 2) oss << "value";
    oss << "_" << mean;
    this->set_name (oss.str());
  }

  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
#ifdef CGAL_CLASSIFICATION_PRECOMPUTE_FEATURES
    return color_feature[pt_index];
#else
    HSV_Color c = Classification::rgb_to_hsv (get(color_map, *(input.begin()+pt_index)));
    return std::exp (-(c[m_channel] - m_mean) * (c[m_channel] - m_mean) / (2. * m_sd * m_sd));
#endif
  }

  /// \endcond
};

} // namespace Feature

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURE_HSV_H
