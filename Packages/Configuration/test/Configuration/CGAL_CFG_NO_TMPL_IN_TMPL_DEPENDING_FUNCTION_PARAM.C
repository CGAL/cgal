// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_TMPL_IN_TMPL_DEPENDING_FUNCTION_PARAM.C
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_TMPL_IN_TMPL_DEPENDING_FUNCTION_PARAM.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| G++ 2.95.2 has problems with member functions implemented outside of
//| the class body if this member function has a parameter type that is
//| dependant on a template in the template parameter list of the class. A
//| workaround would be to implement the member function inline in the class.
//| The following definition is set if this error error occurs.

template < template < class T> class T_HDS>
struct Container {
    typedef T_HDS<int>           HDS;
    typedef typename HDS::Handle Handle;
    void foo( Handle h);
    //void foo( Handle h) {} // workaround: implement foo inline.
};

template < template < class T> class T_HDS >
void Container<T_HDS>::foo( Handle e) {}

template <class T>
struct A {
    typedef int* Handle;
};

int main() {
    Container<A> c;
    c.foo(0);
    return 0;
}

// EOF //
