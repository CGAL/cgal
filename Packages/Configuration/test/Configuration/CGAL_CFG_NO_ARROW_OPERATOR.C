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
// release       : $CGAL_Revision: V $
// release_date  : $CGAL_Date: 1998, July 17 $
//
// file          : config/testfiles/CGAL_CFG_NO_ARROW_OPERATOR.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_ARROW_OPERATOR.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The arrow operator 'operator->()' could not be overloaded by some
//| compilers within a class template (SunPro CC 4.2 complains if the 
//| return value is not a struct). The following definition is set 
//| if it cannot be overloaded. Note that the arrow operator is mandatory
//| for iterators according to the Dec. 1996 C++ Standard draft.

#include <assert.h>

struct A {
    int i;
    A( int _i) : i(_i) {}
};

template< class T>
struct B {
    T* a;
    B( T* _a) : a(_a) {}
    T* operator->() {
	return a;
    }
    T  operator*() {
	return *a;
    }
};

int main() {
    A a( 5);
    B<A> b( &a);
    assert( b->i == 5);
    int i = 3;
    B<int> c(&i);
    assert( *c == 3);
    return 0;
}

// EOF //
