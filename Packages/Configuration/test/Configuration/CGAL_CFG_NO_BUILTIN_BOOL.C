// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: V $
// release_date  : $CGAL_Date: 1998, July 17 $
//
// file          : config/testfiles/CGAL_CFG_NO_BUILTIN_BOOL.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_BUILTIN_BOOL.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler whether it knows
// the keyword explicit or not.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler doesn't know the keyword bool, the flag
//| CGAL_CFG_NO_BUILTIN_BOOL is set.

int main()
{
  bool b = true;

  return 0;
}

// EOF //
