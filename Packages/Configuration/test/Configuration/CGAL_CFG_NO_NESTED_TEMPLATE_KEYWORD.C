// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_NESTED_TEMPLATE_KEYWORD.C
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_NESTED_TEMPLATE_KEYWORD.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set, if a compiler does not accept the template
//| keyword when referring to nested template classes.
//| E.g. if the templated class C is defined within class A,
//| one refers to it by A::template C< >.


struct A {
    template < class T >
    struct C {};
};

template < class T >
struct B {
    typedef typename T::template C< int > C; 
};

int main() {
    B< A > b;
    (void)(b);
    return 0;
}

// EOF //
