// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  :
//
// file          : config/testfiles/CGAL_CFG_NO_TEMPLATE_FRIEND_DISTINCTION.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_TEMPLATE_FRIEND_DISTINCTION.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| Checks whether the compiler wants to have a <> in friend declarations
//| of template functions.

template < class T >
inline int y(T t);

struct A {
  A(int i_) : i(i_) {}
  friend int y<>(A);
private:
  int i;
};

template < class T >
int y(T t) { return t.i - 1; }

int main()
{
  A a(3);
  return (y(a) == 2)?0:1;
}

