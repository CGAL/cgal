// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 01
//
// file          : include/CGAL/bops_dcel_element.h
// package       : bops (2.2)
// source        : include/CGAL/bops_dcel_element.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL__DCEL_ELEMENT_H
#define CGAL__DCEL_ELEMENT_H

#include <CGAL/bops_dcel_defs.h>

CGAL_BEGIN_NAMESPACE

/*
  ELEMENT in the DCEL:
  ---------------------
*/


class _Dcel_element_type {
public:
  _Dcel_element_type( int ind, _Dcel_Color col = _NO_COLOR )
        : _index(ind), _color(col) {}
  _Dcel_element_type( _Dcel_Color col = _NO_COLOR )
        : _index(-1), _color(col) {}

  _Dcel_Color    color() const  { return _color; }
  _Dcel_Color    set_color(const _Dcel_Color& c) {
    _color= c;
    return c;
  }
  bool    has_color(const _Dcel_Color& c) const {
    return color()& c ? true : false;
  }

  _Dcel_Color&   color()   { return _color; }
  int              index() const  { return _index; }

private:
  int _index;
  _Dcel_Color _color;
};

CGAL_END_NAMESPACE

#endif /* CGAL__DCEL_ELEMENT_H */
