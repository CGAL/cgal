//  -*- Mode: c++ -*-
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
// release       : $CGAL_Revision: CGAL-1.0 $
// release_date  : $CGAL_Date: 1998/09/12 $
//
// file          : demo/BooleanOperations/include/CGAL/IO-old/Color.h
// source        : demo/BooleanOperations/include/CGAL/IO-old/Color.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                        Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

#ifndef CGAL_COLOR_H
#define CGAL_COLOR_H


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

  bool operator==(const Color &c)
  {
    return ( (_red == c.red()) &&
             (_green == c.green()) &&
             (_blue == c.blue()) );
  }

  bool operator!=(Color &c)
  {
    return !( (*this) == c);
  }

private:
  int _red;
  int _green;
  int _blue;
};

const Color BLACK  = Color(0, 0, 0);
const Color WHITE  = Color(255, 255, 255);
const Color RED    = Color(255, 0, 0);
const Color GREEN  = Color(0, 255, 0);
const Color BLUE   = Color(0, 0, 255);
const Color VIOLET = Color(255, 0, 255);
const Color ORANGE = Color(255, 170, 0);
#endif  // CGAL_COLOR_H

