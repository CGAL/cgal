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
// release       : $CGAL_Revision: CGAL-0.9-I-04 $
// release_date  : $CGAL_Date: 1997/12/15 $
//
// file          : include/CGAL/IO/Color.h
// source        :
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Hervé Brönnimann
//
// coordinator   : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// ============================================================================


#ifndef CGAL_COLOR_H
#define CGAL_COLOR_H

class CGAL_Color {
public:
  CGAL_Color() {}
  CGAL_Color(int red, int green, int blue)
    : _red(red), _green(green), _blue(blue)
  {}

  int r() const {return _red;}
  int g() const {return _green;}
  int b() const {return _blue;}

  int red() const {return _red;}
  int green() const {return _green;}
  int blue() const {return _blue;}

  bool operator==(const CGAL_Color &c) const
  {
    return ( (red() == c.red()) &&
             (green() == c.green()) &&
             (blue() == c.blue()) );
  }

  bool operator!=(CGAL_Color &c) const
  {
    return !( (*this) == c);
  }

  CGAL_Color& operator=(const CGAL_Color &c)
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


extern const CGAL_Color CGAL_BLACK;
extern const CGAL_Color CGAL_WHITE;
extern const CGAL_Color CGAL_RED;
extern const CGAL_Color CGAL_GREEN;
extern const CGAL_Color CGAL_BLUE;
extern const CGAL_Color CGAL_VIOLET;
extern const CGAL_Color CGAL_ORANGE;


#endif  // CGAL_COLOR_H
