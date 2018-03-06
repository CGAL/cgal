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

#ifndef CGAL_CLASSIFICATION_COLOR_H
#define CGAL_CLASSIFICATION_COLOR_H

#include <CGAL/license/Classification.h>
#include <CGAL/number_utils.h>
#include <CGAL/int.h>
#include <CGAL/array.h>

namespace CGAL {
namespace Classification {

  /*!
    \ingroup PkgClassificationColor

    %Color described in red/green/blue space. Each component is stored
    as an unsigned char ranging from 0 (no color) to 255 (full color).
   */  
typedef CGAL::cpp11::array<unsigned char, 3> RGB_Color;
  /*!
    \ingroup PkgClassificationColor

    %Color described in hue/saturation/value space. Each component is stored
    as a float:

    - `hue` ranges from 0° to 360° (corresponding to the color tint)
    - `saturation` ranges from 0.0 (gray) to 100.0 (full saturation)
    - `value` ranges from 0.0 (black) to 100.0 (white)
   */  
typedef CGAL::cpp11::array<float, 3> HSV_Color;


  /// \cond SKIP_IN_MANUAL
inline HSV_Color rgb_to_hsv (const RGB_Color& c)
{
  double r = (double)(c[0]) / 255.;
  double g = (double)(c[1]) / 255.;
  double b = (double)(c[2]) / 255.;
  double Cmax = (std::max) (r, (std::max) (g, b));
  double Cmin = (std::min) (r, (std::min) (g, b));
  double delta = Cmax - Cmin;
  double H = 0.;
  
  if (delta != 0.)
  {
    if (Cmax == r)
      H = 60. * ((g - b) / delta);
    else if (Cmax == g)
      H = 60. * (((b - r) / delta) + 2.);
    else
      H = 60. * (((r - g) / delta) + 4.);
  }
  if (H < 0.) H += 360.;
  double S = (Cmax == 0. ? 0. : 100. * (delta / Cmax));
  double V = 100. * Cmax;
  HSV_Color out = {{ float(H), float(S), float(V) }};
  return out;
}

inline RGB_Color hsv_to_rgb (const HSV_Color& c)
{
  double h = c[0];
  double s = c[1];
  double v = c[2];
  
  s /= 100.;
  v /= 100.;
  double C = v*s;
  int hh = (int)(h/60.);
  double X = C * (1-CGAL::abs (hh % 2 - 1));
  double r = 0, g = 0, b = 0;
  
  if( hh>=0 && hh<1 )
  {
    r = C;
    g = X;
  }
  else if( hh>=1 && hh<2 )
  {
    r = X;
    g = C;
  }
  else if( hh>=2 && hh<3 )
  {
    g = C;
    b = X;
  }
  else if( hh>=3 && hh<4 )
  {
    g = X;
    b = C;
  }
  else if( hh>=4 && hh<5 )
  {
    r = X;
    b = C;
  }
  else
  {
    r = C;
    b = X;
  }
  double m = v-C;
  r += m;
  g += m;
  b += m;
  r *= 255.0;
  g *= 255.0;
  b *= 255.0;

  RGB_Color out = {{ (unsigned char)r, (unsigned char)g, (unsigned char)b }};
  return out;
}
  /// \endcond

} // namespace Classification
} // namespace CGAL



#endif // CGAL_CLASSIFICATION_COLOR_H
