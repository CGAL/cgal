// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_KOENIG_LOOKUP.C
// source        :
// revision      : 1
// revision_date : 17 Oct 1999
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set if the compiler doesn't support the operator Koenig
//| lookup. That is, it does not search in the namespace of the arguments for
//| the function.

namespace A {
  struct S {};

  void foo(const S &) {}
}

namespace B {
//  using A::foo;

  void bar()
  {
    A::S s;
    foo(s);
  }

}

int main()
{
    B::bar();
    return 0;
}


// EOF //
