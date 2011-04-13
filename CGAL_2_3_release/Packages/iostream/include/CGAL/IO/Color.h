// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/IO/Color.h
// package       : iostream
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : Mariette Yvinec (Mariette.Yvinec@sophia.inria.fr)
//
// ============================================================================

#include <CGAL/config.h>

#ifndef CGAL_COLOR_H
#define CGAL_COLOR_H

CGAL_BEGIN_NAMESPACE

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


extern const Color BLACK  ;
extern const Color WHITE  ;
extern const Color GRAY  ;

extern const Color RED    ;
extern const Color GREEN  ;

extern const Color DEEPBLUE  ;
extern const Color BLUE   ;
extern const Color PURPLE ;
extern const Color VIOLET ;

extern const Color ORANGE ;
extern const Color YELLOW ;


CGAL_END_NAMESPACE

#endif  // CGAL_COLOR_H
