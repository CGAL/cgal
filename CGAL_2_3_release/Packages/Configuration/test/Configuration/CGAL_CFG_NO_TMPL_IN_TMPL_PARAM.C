// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_TMPL_IN_TMPL_PARAM.C
// author(s)     : Lutz Kettner
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_TMPL_IN_TMPL_PARAM.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| Nested templates in template parameter, such as 'template <
//| template <class T> class A>' are not supported by any compiler. 
//| The following definition is set if they are not supported.

template< class X>
struct A {
    X i;
    A( X j) : i(j) {}
};

template< template < class T> class B>
struct C {
    B<int> b;
    C( int i) : b(i) {}
};

template< class N, template < class T> class B>
struct D {
    B<N> b;
    D( N i) : b(i) {}
};

int main() {
    C<A> c(1);
    D< double, A> d( 3.8);
    (void) c;
    (void) d;
    return 0;
}
