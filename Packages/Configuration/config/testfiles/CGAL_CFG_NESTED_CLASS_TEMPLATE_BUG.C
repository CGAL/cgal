// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// file          : config/testfiles/CGAL_CFG_NESTED_CLASS_TEMPLATE_BUG.C
// package       : Configuration (2.3)
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NESTED_CLASS_TEMPLATE_BUG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag ist set, if a compiler cannot parse an instantiation of
//| a nested class template B which is defined in another class
//| template A. This bug shows up on VC7. A workaround is to pass the
//| instantiation of A to a helper class which extracts B. (see
//| CGAL/Kernel_d/Iso_box_d.h)

struct II {};

template < bool b > struct A;
template <> struct A<true> { 
  template <typename T> struct B {
    typedef int type;
  };
};
template <> struct A<false> { 
  template <typename T> struct B {
    typedef II type;
  };
};

template <bool b>
struct C {
  typedef typename A<b>::template B<int>::type D;
};

int main()
{
  C<true>::D i = 0;
  ++i;
  return 0;
}
