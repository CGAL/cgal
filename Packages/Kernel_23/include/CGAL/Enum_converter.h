// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// file          : Enum_converter.h
// package       : Kernel_23
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
// ======================================================================

#ifndef CGAL_ENUM_CONVERTER_H
#define CGAL_ENUM_CONVERTER_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

struct Enum_converter
{
  bool              operator()(bool b) const { return b; }

  Sign              operator()(Sign s) const { return s; }

  Oriented_side     operator()(Oriented_side os) const { return os; }

  Bounded_side      operator()(Bounded_side bs) const { return bs; }

  Comparison_result operator()(Comparison_result cr) const {
    return cr;
  }

  Angle operator()(Angle a) const { return a; }
};


CGAL_END_NAMESPACE


#endif // CGAL_ENUM_CONVERTER_H
