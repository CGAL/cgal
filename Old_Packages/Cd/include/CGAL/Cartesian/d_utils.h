// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/d_utils.h
// revision      : $Revision$
// revision_date : $Date$
// author        : Herve Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_D_UTILS_H
#define CGAL_CARTESIAN_D_UTILS_H

#include <iterator>
#include <algorithm>
#include <functional>

CGAL_BEGIN_NAMESPACE

// Small object utility for printing objects
// with separators depending on the ostream mode
template < class Object >
struct print_d
{
  char *       _separator;
  std::ostream*_os;
  bool         _print_sep;

  print_d(std::ostream *os) : _os(os), _print_sep(false)
  {
    if (_os->iword(IO::mode)==IO::ASCII) _separator = " ";
    else if (_os->iword(IO::mode)==IO::BINARY) _separator = "";
    else _separator = ", ";
  }
  void reset()
  {
    _print_sep = false;
  }
  void operator()(const Object &x) {
    if (_print_sep && _os->iword(IO::mode) != IO::BINARY)
      (*_os) << _separator;
    (*_os) << x;
    _print_sep = true;
  }
  void operator()(const int &x) {
    if (_print_sep && _os->iword(IO::mode) != IO::BINARY)
      (*_os) << _separator;
    (*_os) << x;
    _print_sep = true;
  }
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_D_UTILS_H
