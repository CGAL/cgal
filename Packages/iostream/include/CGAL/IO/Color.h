// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium

// This software and related documentation is part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation is provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind. 
//
// Every use of CGAL requires a license. 
//
// Academic research and teaching license
// - For academic research and teaching purposes, permission to use and copy
//   the software and its documentation is hereby granted free of charge,
//   provided that it is not a component of a commercial product, and this
//   notice appears in all copies of the software and related documentation. 
//
// Commercial licenses
// - A commercial license is available through Algorithmic Solutions, who also
//   markets LEDA (http://www.algorithmic-solutions.de). 
// - Commercial users may apply for an evaluation license by writing to
//   Algorithmic Solutions (contact@algorithmic-solutions.com). 
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Free University of Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// release       : CGAL-2.1
// release_date  : 2000, January 11
//
// file          : include/CGAL/IO/Color.h
// package       : iostream (2.5)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Hervé Brönnimann
//
// coordinator   : Herve Bronnimann
//
// email         : cgal@cs.uu.nl
//
// ======================================================================


#include <CGAL/config.h>
#ifndef CGAL_COLOR_H
#define CGAL_COLOR_H

CGAL_BEGIN_NAMESPACE

class Color {
public:
  Color() {}
  Color(unsigned char red, unsigned char green, unsigned char blue, unsigned char alpha = 120)
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
