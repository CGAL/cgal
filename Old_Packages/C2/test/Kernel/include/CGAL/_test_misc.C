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
// source        : test_kernel_misc.fw
// file          : _test_misc.C
// revision      : 2.1
// revision_date : 05 Aug 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
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
