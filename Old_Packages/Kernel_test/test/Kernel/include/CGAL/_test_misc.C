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
// file          : _test_misc.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_MISC_C
#define CGAL__TEST_MISC_C

#include <CGAL/_test_cls_object.C>


template <class R> bool _test_misc(const R& r);


template <class R>
bool
_test_misc(const R& r)
{
  return
    _test_cls_object(r);
}


#endif // CGAL__TEST_MISC_C
