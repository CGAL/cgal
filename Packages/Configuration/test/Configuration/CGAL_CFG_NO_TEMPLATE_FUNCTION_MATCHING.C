// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_TEMPLATE_FUNCTION_MATCHING.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_TEMPLATE_FUNCTION_MATCHING.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| Some compilers (like g++ 2.7.2) follow the old rules for argument matching
//| stated in the ARM: template functions have to match exactly. In particular,
//| derived classes don't match to base class parameters in a template
//| function. The following flag is set if this is the case.

template <class T>
class B
{
  public: 
    T data;
};

template <class T>
class D : public B<T>
{
};

template <class T>
void 
foo(const B<T>& b1, const B<T>& b2)
{
  b1.data + b2.data; 
}

int
main()
{
  D<int>  d1,d2;
  foo(d1,d2);

  return 0;
}

