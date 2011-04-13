// ======================================================================
//
// Copyright (c) 1999,2000 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.2-I-2 $
// release_date  : $CGAL_Date: 2000/02/04 $
//
// file          : config/testfiles/CGAL_CFG_NO_LONGNAME_PROBLEM.C
// package       : Configuration (2.1)
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
Wwwwwwwwwooooooooo_vvvvvvveeeeeerrrryyyy_llllooooonnnnnnggggg_nnnnnaaaammmmmeeeeWwwwwwwwwooooooooo_vvvvvvveeeeeerrrryyyy_llllooooonnnnnnggggg_nnnnnaaaammmmmeeeeWwwwwwwwwooooooooo_vvvvvvveeeeeerrrryyyy_llllooooonnnnnnggggg_nnnnnaaaammmmmeeeeWwwwwwwwwooooooooo_vvvvvvveeeeeerrrryyyy_llllooooonnnnnnggggg_nnnnnaaaammmmmeeeeWwwwwwwwwooooooooo_vvvvvvveeeeeerrrryyyy_llllooooonnnnnnggggg_nnnnnaaaammmmmeeeeWwwwwwwwwooooooooo_vvvvvvveeeeeerrrryyyy_llllooooonnnnnnggggg_nnnnnaaaammmmmeeee

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
  (void) a;
  return 0;
}
