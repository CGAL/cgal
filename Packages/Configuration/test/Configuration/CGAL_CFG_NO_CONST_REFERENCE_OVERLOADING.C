// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_CONST_REFERENCE_OVERLOADING.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_CONST_REFERENCE_OVERLOADING.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The Gnu g++ 2.7.2.1 isn't able to distinguish cleanly between 
//| overloading with a reference and a const reference parameter
//| if he looks fpor a match with a derived class and reports an 
//| ambigious overloading resolution. One workaround could be to
//| remove the const reference declaration. g++ will warn a converion
//| from const X& to X& but compiles. Another workaround would be the
//| explicit cast to the base class (if known).
//| The following definition is set if the compiler fails to
//| distinguish those functions.


#include <assert.h>

struct Y {
    int i;
    Y(int j) : i(j) {}
};

struct X : public Y {
    X(int j) : Y(j) {}
};

template < class T>
struct A {
    T&       foo( T& t)       const { return t; }
    const T& foo( const T& t) const { return t; }
};

int main() {
    A<Y> a;
    X       i(42);
    const X j(43);
    assert( (a.foo(i)).i == 42);
    assert( (a.foo(j)).i == 43);
    return 0;
}

// EOF //


