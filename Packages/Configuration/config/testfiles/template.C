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
// release       : $CGAL_Revision: CGAL-0.9-I-03 $
// release_date  : $CGAL_Date: 1997/11/13 $
//
// file          : config/testfiles/template.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ============================================================================

#include "template.h"

template <class T>
T A<T>::square(const T& t)
{
  return t*t;
}

// EOF //
