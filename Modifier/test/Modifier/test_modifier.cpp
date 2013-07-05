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
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : test_modifier.C
// chapter       : $CGAL_Chapter: Protected access to the internal repr.$
// package       : $CGAL_Package: Modifier 1.2 (07 Apr 1999) $
// source        : modifier.fw
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// Protected access to the internal representation.
// ============================================================================


#include <CGAL/basic.h>
#include <CGAL/Modifier_base.h>

#include <cassert>

using CGAL::Modifier_base;

class A {
    int i;  // protected internal representation
public:
    void delegate( Modifier_base<int>& modifier) {
        modifier(i);
        // check validity
        CGAL_postcondition( i > 0);
    }
    int get_i() const { return i;}  // read access
};

struct Modifier : public Modifier_base<int> {
    void operator()( int& rep) { rep = 42;}
};

int main() {
    A a;
    Modifier m;
    a.delegate(m);
    assert( a.get_i() == 42);
    return 0;
}
// EOF //
