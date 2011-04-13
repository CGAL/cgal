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
// release_date  : 1999, July 28
//
// file          : copy_n.h
// package       : STL_Extension (2.7)
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
// source        : stl_extension.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// STL like copy that copies n elements
// ======================================================================

#ifndef CGAL_COPY_N_H
#define CGAL_COPY_N_H 1
#ifndef CGAL_PROTECT_CSTDDEF
#include <cstddef>
#define CGAL_PROTECT_CSTDDEF
#endif

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
