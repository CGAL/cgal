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

#ifndef CGAL_CLASSIFICATION_FEATURE_COLOR_CHANNEL_H
#define CGAL_CLASSIFICATION_FEATURE_COLOR_CHANNEL_H

#include <CGAL/license/Classification.h>

#include <vector>

#include <CGAL/Classification/Color.h>
#include <CGAL/Classification/Feature_base.h>

namespace CGAL {

namespace Classification {

namespace Feature {
  
template <typename GeomTraits, typename PointRange, typename ColorMap>
class Color_channel : public Feature_base
{
public:

  /// Selected channel.
  enum Channel
  {
    HUE = 0, ///< 0
    SATURATION = 1, ///< 1
    VALUE = 2 ///< 2
  };

private:
  
  typedef typename Classification::RGB_Color RGB_Color;
  typedef typename Classification::HSV_Color HSV_Color;

  const PointRange& input;
  ColorMap color_map;
  Channel m_channel;
  
public:
  
  Color_channel (const PointRange& input,
                 ColorMap color_map,
                 Channel channel)
    : input(input), color_map(color_map), m_channel (channel)
  {
    if (channel == HUE) this->set_name ("color_hue");
    else if (channel == SATURATION) this->set_name ("color_saturation");
    else if (channel == VALUE) this->set_name ("color_value");
  }

  virtual float value (std::size_t pt_index)
  {
    HSV_Color c = Classification::rgb_to_hsv (get(color_map, *(input.begin()+pt_index)));
    return c[std::size_t(m_channel)];
  }

  /// \endcond
};

} // namespace Feature

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_FEATURE_COLOR_CHANNEL_H
