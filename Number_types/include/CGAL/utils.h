// Copyright (c) 2006-2008  Max-Planck-Institute Saarbruecken (Germany)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer <hemmer@mpi-sb.mpg.de>

#ifndef CGAL_UTILS_H
#define CGAL_UTILS_H

#include <CGAL/utils_classes.h>

namespace CGAL {

// use Min
template< class T >
inline T min BOOST_PREVENT_MACRO_SUBSTITUTION
( const T& x , const T& y) {
    return Min<T>()( x , y );
}

template< class T , class Compare >
inline T min BOOST_PREVENT_MACRO_SUBSTITUTION
( const T& x , const T& y, const Compare& c) {
    return Min<T, Compare>(c)( x , y );
}

// use Max
template< class T >
inline T max BOOST_PREVENT_MACRO_SUBSTITUTION
( const T& x , const T& y) {
    return Max<T>()( x , y );
}

template< class T , class Compare >
inline T max BOOST_PREVENT_MACRO_SUBSTITUTION
( const T& x , const T& y, const Compare& c) {
    return Max<T, Compare>(c)( x , y );
}

// use Is_valid
template< class T >
inline bool is_valid( const T& x ) {
  return Is_valid< T >()( x );
}

} //namespace CGAL

#endif // CGAL_UTILS_H
