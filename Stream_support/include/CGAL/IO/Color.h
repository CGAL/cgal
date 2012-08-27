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
// 
//
// Author(s)     : Andreas Fabri

#include <CGAL/config.h>

#ifndef CGAL_COLOR_H
#define CGAL_COLOR_H

namespace CGAL {

class Color {
public:
  Color() {}
  Color(unsigned char red, 
	unsigned char green, 
	unsigned char blue, 
	unsigned char alpha = 120)
    : _red(red), _green(green), _blue(blue), _alpha(alpha)
  {}

  unsigned char r() const {return _red;}
  unsigned char g() const {return _green;}
  unsigned char b() const {return _blue;}

  unsigned char red() const {return _red;}
  unsigned char green() const {return _green;}
  unsigned char blue() const {return _blue;}
  unsigned char alpha() const {return _alpha;}
  void set_alpha(unsigned char a) {_alpha=a;}
  bool operator==(const Color &c) const
  {
    return ( (red() == c.red()) &&
             (green() == c.green()) &&
             (blue() == c.blue()) );
  }

  bool operator!=(const Color &c) const
  {
    return !( (*this) == c);
  }

  Color& operator=(const Color &c)
  {
    _red = c.red();
    _green = c.green();
    _blue = c.blue();
    _alpha = c.alpha();
    return *this;
  }

private:
  unsigned char _red;
  unsigned char _green;
  unsigned char _blue;
  unsigned char _alpha;
};


CGAL_EXPORT extern const Color BLACK  ;
CGAL_EXPORT extern const Color WHITE  ;
CGAL_EXPORT extern const Color GRAY  ;

CGAL_EXPORT extern const Color RED    ;
CGAL_EXPORT extern const Color GREEN  ;

CGAL_EXPORT extern const Color DEEPBLUE  ;
CGAL_EXPORT extern const Color BLUE   ;
CGAL_EXPORT extern const Color PURPLE ;
CGAL_EXPORT extern const Color VIOLET ;

CGAL_EXPORT extern const Color ORANGE ;
CGAL_EXPORT extern const Color YELLOW ;


} //namespace CGAL

#endif  // CGAL_COLOR_H
