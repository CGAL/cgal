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
// file          : Random_access_value_adaptor.h
// package       : STL_Extension (2.7)
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
// source        : stl_extension.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// Random Access Value Adaptor provides random access for sequences.
// ======================================================================

#ifndef CGAL_RANDOM_ACCESS_VALUE_ADAPTOR_H
#define CGAL_RANDOM_ACCESS_VALUE_ADAPTOR_H 1
#ifndef CGAL_RANDOM_ACCESS_ADAPTOR_H
#include <CGAL/Random_access_adaptor.h>
#endif // CGAL_RANDOM_ACCESS_ADAPTOR_H

CGAL_BEGIN_NAMESPACE

template < class IC, class T >
class Random_access_value_adaptor : public Random_access_adaptor<IC> {
public:
    Random_access_value_adaptor() {}
        // invalid index.

    Random_access_value_adaptor( const IC& i)
        : Random_access_adaptor<IC>(i) {}
        // empty random access index initialized to start at i.

    Random_access_value_adaptor( const IC& i, const IC& j)
        : Random_access_adaptor<IC>(i,j) {}
        // random access index initialized with range [i,j).

// OPERATIONS

    T& operator[]( size_type n) const {
        // returns inverse index of k.
        return *(Random_access_adaptor<IC>::operator[](n));
    }
};

CGAL_END_NAMESPACE
#endif // CGAL_RANDOM_ACCESS_VALUE_ADAPTOR_H //
// EOF //
