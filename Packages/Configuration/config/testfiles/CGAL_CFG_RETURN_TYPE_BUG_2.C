// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_RETURN_TYPE_BUG_2.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_RETURN_TYPE_BUG_2.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The following flag is set if a compiler does not allow the return
//| type A<T>::B of a function or class method.

template <class T>
class A
{
  public:
    typedef int B;
    B f();
    B g() { return B(); } // this doesn't give problems
};

// class method defined outside the class
template <class T>
A<T>::B A<T>::f()
{
  return A<T>::B();
}

// global function with return type A<T>::B
template <class T>
A<T>::B f(const T&)
{
  return A<T>::B();
}

int main()
{
  return 0;
}

// EOF //
