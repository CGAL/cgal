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
// file          : config/testfiles/CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The parameter types of member functions might contain a scope
//| operator. This works as long as the member function is implemented
//| inline in the class. If the member function is implemented external
//| not all compilers are able to parse the scope operators correctly.
//| The following definition is set if the compiler fails parsing.


#include <assert.h>

template < class T>
struct A {
    typedef T X;
};

template< class T>
struct B {
    typename T::X  foo( typename T::X i);
};


template< class T>
typename T::X  B<T>::foo( typename T::X i) {
    return i + 2;
}

int main() {
    B<A<int> > b;
    assert( b.foo(40) == 42);
    return 0;
}

// EOF //

