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
// release       : $CGAL_Revision: CGAL-0.9-I-03 $
// release_date  : $CGAL_Date: 1997/11/13 $
//
// file          : config/testfiles/CGAL_CFG_NO_LAZY_INSTANTIATION.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ============================================================================

// CGAL_CFG_NO_LAZY_INSTANTIATION.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| Implicit instantiation of a class template does only instantiate member
//| functions when needed (Dec. 1996 C++ Standard draft, 14.7.1). 
//| This implies that member functions that are not instanciated in a 
//| certain context are allowed to use functionality from the template 
//| arguments that are not provided by the actual argument. For example
//| the Gnu g++ 2.7.2 does not comply to this. The following definition
//| is set if the implicit instantiation does not work in this lazy fashion.

#include <assert.h>

struct A {
    int i;
    A( int _i) : i(_i) {}
    A& operator++() {
	++i;
	return *this;
    }
};

template< class T>
struct B {
    T a;
    B( T _a) : a(_a) {}
    B& operator++() {
	++a;
	return *this;
    }
    B& operator--() {  // This one is not used nor supported by A.
	--a;
	return *this;
    }
};

int main() {
    A a( 5);
    B<A> b( a);
    ++b;
    assert( b.a.i == 6);
    return 0;
}

// EOF //
