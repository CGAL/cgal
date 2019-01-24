// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_COLOR_H
#define CGAL_COLOR_H

#include <CGAL/config.h>
#include <CGAL/array.h>

namespace CGAL {


/*!
  \ingroup PkgStreamSupportRef

  An object of the class `Color` is a color available for drawing
  operations in many \cgal output streams.  Each color is defined by a
  4 unsigned chars `(r,g,b,a)` with 0 \f$\le\f$ r,g,b,a \f$ \le \f$
  255, the so-called <I>rgba-value</I> of the color.

  The alpha parameter (representing transparency) is often ignored and
  left to its default value (255 = no transparency), which is why we
  often refer to the <I>rgb-value</I> of the color.

  \sa `CGAL::Geomview_stream` 

*/

class Color
{
private:

  cpp11::array<unsigned char, 4> m_data;
  
public:

  /// \name Creation 
  /// @{

  /*!
    creates a color with rgba-value `(0,0,0,255)`, i.e.\ black. 
  */ 
  Color()
  {
    set_rgb (0, 0, 0, 255);
  }

  /*!
    creates a color with rgba-value `(red,green,blue,alpha)`. 
  */ 
  Color(unsigned char red,
        unsigned char green,
        unsigned char blue,
        unsigned char alpha = 255)
  {
    set_rgb (red, green, blue, alpha);
  }

  /// @} 

  /// \name Component Access
  /// @{
  
  /*!
    returns the red component. 
  */ 
  unsigned char red() const { return m_data[0]; }

  /*!
    returns a reference on the red component. 
  */ 
  unsigned char& red() { return m_data[0]; }

  /*!
    returns the green component. 
  */ 
  unsigned char green() const { return m_data[1]; }

  /*!
    returns a reference on the green component. 
  */ 
  unsigned char& green() { return m_data[1]; }

  /*!
    returns the blue component. 
  */ 
  unsigned char blue() const { return m_data[2]; }
  
  /*!
    returns a reference on the blue component. 
  */ 
  unsigned char& blue() { return m_data[2]; }
  
  /*!
    returns the alpha component. 
  */ 
  unsigned char alpha() const { return m_data[3]; }

  /*!
    returns a reference on the alpha component. 
  */ 
  unsigned char& alpha() { return m_data[3]; }

  /// @} 

  /// \name Array Access
  /// @{

  /*!
    returns the i^{th} component of the rgb color (the 0^{th} is red,
    the 1^{st} is blue, etc.).
  */ 
  unsigned char operator[] (std::size_t i) const { return m_data[i]; }
  
  /*!
    returns a reference on the i^{th} component of `c` (the 0^{th} is
    red, the 1^{st} is blue, etc.).
  */ 
  unsigned char& operator[] (std::size_t i)  { return m_data[i]; }

  /*!
    returns the array with rgba values.
  */ 
  const cpp11::array<unsigned char, 4>& to_rgba() const { return m_data; }

  /*!
    returns the array with rgb values.
  */ 
  const cpp11::array<unsigned char, 3>& to_rgb() const
  {
    return reinterpret_cast<const cpp11::array<unsigned char, 3>&>(m_data);
  }

  /*!  
    computes the hsv (hue, saturation, value) values and returns an
    array representing them as float values between 0 and 1.
  */ 
  cpp11::array<double, 3> to_hsv() const
  {
    double r = (double)(m_data[0]) / 255.;
    double g = (double)(m_data[1]) / 255.;
    double b = (double)(m_data[2]) / 255.;
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

    return make_array(H,S,V);
  }
  
  /// \name Modification 
  /// @{

  /*!  
    replaces the rgb values of the colors by the one given as parameters.
  */ 
  void set_rgb (unsigned char red,
                unsigned char green,
                unsigned char blue,
                unsigned char alpha = 255)
  {
    m_data[0] = red;
    m_data[1] = green;
    m_data[2] = blue;
    m_data[3] = alpha;
  }

  /*!  
    replaces the rgb values of the colors by the conversion to rgb of
    the hsv values given as parameters.

    Double values given as parameters should take range between 0 and 1.
  */ 
  void set_hsv (double hue,
                double saturation,
                double value,
                unsigned char alpha = 255)
  {
    saturation /= 100.;
    value /= 100.;
    double C = value*saturation;
    int hh = (int)(hue/60.);
    double X = C * (1-std::abs (hh % 2 - 1));
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
    double m = value-C;
    r += m;
    g += m;
    b += m;
    r *= 255.0;
    g *= 255.0;
    b *= 255.0;

    m_data[0] = (unsigned char)r;
    m_data[1] = (unsigned char)g;
    m_data[2] = (unsigned char)b;
    m_data[3] = alpha;
  }

  /// @}

};
  

/*!
  Constructs Color(0,0,0).
  \relates Color
*/
inline Color black() { return Color(0,0,0); }

/*!
  Constructs Color(0,0,255).
  \relates Color
*/
inline Color blue() { return Color(0,0,255); }

/*!
  Constructs Color(10,0,100).
  \relates Color
*/
inline Color deep_blue() { return Color(10,0,100); }

/*!
  Constructs Color(100,100,100).
  \relates Color
*/
inline Color gray() { return Color(100,100,100); }

/*!
  Constructs Color(0,255,0).
  \relates Color
*/
inline Color green() { return Color(0,255,0); }

/*!
  Constructs Color(235,150,0).
  \relates Color
*/
inline Color orange() { return Color(235,150,0); }

/*!
  Constructs Color(100,0,70).
  \relates Color
*/
inline Color purple() { return Color(100,0,70); }

/*!
  Constructs Color(255,0,0).
  \relates Color
*/
inline Color red() { return Color(255,0,0); }

/*!
  Constructs Color(255,0,255).
  \relates Color
*/
inline Color violet() { return Color(255,0,255); }

/*!
  Constructs Color(255,255,255).
  \relates Color
*/
inline Color white() { return Color(255,255,255); }

/*!
  Constructs Color(255,255,0).
  \relates Color
*/
inline Color yellow() { return Color(255,255,0); }


};


} //namespace CGAL

#endif  // CGAL_COLOR_H
