// ============================================================================
//
// Copyright (c) 1997, 1998, 1999, 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : copy_n.h
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
// package       : $CGAL_Package: STL_Extension $
// source        : stl_extension.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@cs.unc.edu>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// STL like copy that copies n elements
// ============================================================================


// This file is obsolete and exists only for backwards-compatibility.
// Include <CGAL/algorithm.h> instead.

#ifndef CGAL_COPY_N_H
#define CGAL_COPY_N_H 1
#include <cstddef>

// copy_n is usually in the STL as well, but not in the official
// standard. We provide out own copy_n. Only on Gnu g++ 2.8.1
// the hack to work with namespaces gives a name clash, which
// we avoid using the follwing workaround.
#ifndef CGAL_CFG_NO_NAMESPACE
CGAL_BEGIN_NAMESPACE
template <class InputIterator, class Size, class OutputIterator>
OutputIterator copy_n( InputIterator first,
                       Size n,
                       OutputIterator result) {
  // copies the first `n' items from `first' to `result'. Returns
  // the value of `result' after inserting the `n' items.
  while( n--) {
    *result = *first;
    first++;
    result++;
  }
  return result;
}
CGAL_END_NAMESPACE
#endif // CGAL_CFG_NO_NAMESPACE //
#endif // CGAL_COPY_N_H //
// EOF //
