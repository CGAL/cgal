// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_MEMBER_TEMPLATES.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_MEMBER_TEMPLATES.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler whether it supports
// member templates or not.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler doesn't support member templates, the flag
//| CGAL_CFG_NO_MEMBER_TEMPLATES is set.

class A
{
  public:
    template <class T> T square(T t) { return t*t; }
};

int main()
{
  A a;
  int i = a.square(3);  

  return 0;
}

// EOF //
