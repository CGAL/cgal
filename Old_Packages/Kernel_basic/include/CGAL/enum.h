// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : enum.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_ENUM_H
#define CGAL_ENUM_H

CGAL_BEGIN_NAMESPACE

enum  Sign
      {
        NEGATIVE   = -1,
        ZERO,
        POSITIVE
      };

typedef Sign Orientation;

const Orientation  LEFTTURN   = POSITIVE;
const Orientation  RIGHTTURN  = NEGATIVE;

const Orientation  CLOCKWISE  = NEGATIVE;
const Orientation  COUNTERCLOCKWISE = POSITIVE;

const Orientation  COLLINEAR  = ZERO;
const Orientation  COPLANAR   = ZERO;
const Orientation  DEGENERATE = ZERO;

enum  Oriented_side
      {
        ON_NEGATIVE_SIDE = -1,
        ON_ORIENTED_BOUNDARY,
        ON_POSITIVE_SIDE
      };

enum  Bounded_side
      {
        ON_UNBOUNDED_SIDE = -1,
        ON_BOUNDARY,
        ON_BOUNDED_SIDE
      };

enum  Comparison_result
      {
        SMALLER   = -1,
        EQUAL,
        LARGER
      };

CGAL_END_NAMESPACE

#include <CGAL/functions_on_enums.h>

#endif // CGAL_ENUM_H
