// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_LONGNAME_PROBLEM.C
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_LONGNAME_PROBLEM.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler (or assembler or linker) has problems with long names
//| CGAL_CFG_NO_LONGNAME_PROBLEM is set.

#define LONG_NAME \
Wwwwwwwwwooooooooo_vvvvvvveeeeeerrrryyyy_llllooooonnnnnnggggg_nnnnnaaaammmmmeeeeWwwwwwwwwooooooooo_vvvvvvveeeeeerrrryyyy_llllooooonnnnnnggggg_nnnnnaaaammmmmeeee

template < class A >
struct LONG_NAME
{
  LONG_NAME (int i) : a(i) {}
  A a;
};

int main ()
{
  LONG_NAME< LONG_NAME< LONG_NAME< LONG_NAME< LONG_NAME< LONG_NAME<
  LONG_NAME< LONG_NAME< LONG_NAME< LONG_NAME< int > > > > > > > > > > a (1);
  return 0;
}
