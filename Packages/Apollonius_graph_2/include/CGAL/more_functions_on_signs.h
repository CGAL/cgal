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
// release_date  : 
//
// file          : include/CGAL/more_functions_on_signs.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



#ifndef CGAL_MORE_FUNCTIONS_ON_SIGNS_H
#define CGAL_MORE_FUNCTIONS_ON_SIGNS_H

#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE
inline
Sign
sign_of_first_root(const Sign &sb, const Sign &sc)
{
  // this method computes the sign of the second root of the quadratic
  // equation: a x^2 - b x + c == 0. It is assumed that the sign of a
  // is +1.

  if ( sc == NEGATIVE )
    return NEGATIVE;

  if ( sb == POSITIVE )  return sc;
  if ( sb == NEGATIVE )  return NEGATIVE;
  return opposite(sc);
}


inline
Sign
sign_of_second_root(const Sign &sb, const Sign &sc)
{
  // this method computes the sign of the first root of the quadratic
  // equation: a x^2 - b x + c == 0. It is assumed that the sign of a
  // is +1.

  if ( sc == NEGATIVE )
    return POSITIVE;

  if ( sb == POSITIVE )  return POSITIVE;
  if ( sb == NEGATIVE )  return opposite(sc);
  return sc;
}

CGAL_END_NAMESPACE

#endif // CGAL_MORE_FUNCTIONS_ON_SIGNS_H
