// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// A short test program to evaluate a C++ compiler whether it distinguishes
// a template instantiation and a non-templated function with the same
// signature in a friend declaration.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// Perhaps this file should be called CGAL_CFG_NO_EMPTY_TEMPLATE_ARGUMENT_LIST
// ---------------------------------------------------------------------

//| If a compiler doesn't distinguish a template instantiation and 
//| a non-templated function with the same signature, 
//| CGAL_CFG_NO_TEMPLATE_FRIEND_DISTINCTION is set.

int y(int i) { return i + 1; }

template < class T >
int y(T t) { return t - 1; }

int main()
{
  int i = 3;
  if (y(i) == 4 && y<>(i) == 2)
    return 0;
  return 1;
}

// EOF //

