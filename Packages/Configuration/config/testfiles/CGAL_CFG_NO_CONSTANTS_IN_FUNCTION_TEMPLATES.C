// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_CONSTANTS_IN_FUNCTION_TEMPLATES.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_CONSTANTS_IN_FUNCTION_TEMPLATES.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| Constants are not accepted by certain compilers in a template
//| argument list for template functions. The following definition is set 
//| if they are not accepted.


#include <assert.h>

template < class T, int N>
struct A {
    T i;
    A() : i(N) {}
};

template< class T, int N>
int foo( A<T,N> a) {
    return a.i;
}

int main() {
    A<int,42> a;
    assert( foo(a) == 42);
    return 0;
}

// EOF //
