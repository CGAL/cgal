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
// file          : config/testfiles/CGAL_CFG_NO_DEFAULT_TEMPLATE_ARGUMENTS.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ============================================================================

// CGAL_CFG_NO_DEFAULT_TEMPLATE_ARGUMENTS.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| Default template arguments like 'template< class A = int>' are not
//| supported by any compiler. The following definition is set if they
//| are not supported.

template< class T = int>
struct A {
    T i;
    A( T _i) : i(_i) {}
};

template< class S, class T = int>
struct B {
    S i;
    T j;
    B( S _i, T _j) : i(_i), j(_j) {}
};

int main() {
    A<> a( 5);
    B<double> b( 0.4, 8);
    return 0;
}

// EOF //
