// ============================================================================
//
// Copyright (c) 1997,1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-wip $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/workaround_stl.h
// chapter       : $CGAL_Chapter: Configuration $
// package       : $CGAL_Package: Workarounds $
//
// source        : web/workarounds.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sven Schönherr <sven@inf.fu-berlin.de>
//
// coordinator   : Utrecht University (Wieger Wesselink <wieger@cs.ruu.nl>)
//
// implementation: Workarounds for STL
// ============================================================================

#ifndef CGAL_WORKAROUND_STL_H
#define CGAL_WORKAROUND_STL_H 1

// fix for GNU's type unification bug in STL's `distance' function
// ---------------------------------------------------------------
template < class BidirectionalIterator, class Distance >
void
distance( BidirectionalIterator first,
          BidirectionalIterator last,
          Distance& n);

template < class BidirectionalIterator, class Distance >
inline
void
CGAL__distance( BidirectionalIterator first,
                BidirectionalIterator last,
                Distance& n)
{
#ifdef CGAL_STL_GCC
    while ( first != last) {
        ++first;
        ++n;
    }
#else
    distance( first, last, n);
#endif // CGAL_STL_GCC
}

#endif // CGAL_WORKAROUND_STL_H

// ===== EOF ==================================================================
