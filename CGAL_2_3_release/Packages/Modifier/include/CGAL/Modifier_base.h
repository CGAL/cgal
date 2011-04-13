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
// release_date  : 1999, April 07
//
// file          : Modifier_base.h
// package       : Modifier (1.2)
// chapter       : $CGAL_Chapter: Protected access to the internal repr.$
// source        : modifier.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// Modifier providing protected access to internal representations.
// ======================================================================

#ifndef CGAL_MODIFIER_BASE_H
#define CGAL_MODIFIER_BASE_H 1

CGAL_BEGIN_NAMESPACE

template <class R>
class Modifier_base {
    // Abstract Base Class for protected internal access.
    // It defines the common interface for all modifiers.
public:
    typedef R Representation;
    virtual void operator()( R& rep) = 0;
        // Postcondition: `rep' is a valid representation.
    virtual ~Modifier_base() {}
};

CGAL_END_NAMESPACE

#endif // CGAL_MODIFIER_BASE_H //
// EOF //
