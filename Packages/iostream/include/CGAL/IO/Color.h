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
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Hervé Brönnimann
//
// coordinator   : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// ============================================================================


#ifndef CGAL_COLOR_H
#define CGAL_COLOR_H

CGAL_BEGIN_NAMESPACE

class Color {
public:
  Color() {}
  Color(int red, int green, int blue)
    : _red(red), _green(green), _blue(blue)
  {}

  int r() const {return _red;}
  int g() const {return _green;}
  int b() const {return _blue;}

  int red() const {return _red;}
  int green() const {return _green;}
  int blue() const {return _blue;}

  bool operator==(const Color &c) const
  {
    return ( (red() == c.red()) &&
             (green() == c.green()) &&
             (blue() == c.blue()) );
  }

  bool operator!=(Color &c) const
  {
    return !( (*this) == c);
  }

  Color& operator=(const Color &c)
  {
    _red = c.red();
    _green = c.green();
    _blue = c.blue();
    return *this;
  }

private:
  int _red;
  int _green;
  int _blue;
};


extern const Color BLACK;
extern const Color WHITE;
extern const Color RED;
extern const Color GREEN;
extern const Color BLUE;
extern const Color VIOLET;
extern const Color ORANGE;

CGAL_END_NAMESPACE

#endif  // CGAL_COLOR_H
