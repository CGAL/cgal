// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// source        : Test.fw
// file          : _test_sign.h
// revision      : 4.2
// revision_date : 13 Dec 2000 
// author(s)     : Stefan Schirra
//
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de> 
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/basic.h>
#include <cassert>

CGAL_BEGIN_NAMESPACE

template <class NT>
bool
_test_sign(const NT& n)
{
  return (    (sign( NT( 0)) == ZERO)
           && (sign( NT( 1)) == POSITIVE)
           && (sign( NT(-1)) == NEGATIVE) );
}
CGAL_END_NAMESPACE

