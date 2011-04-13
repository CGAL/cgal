// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_LOCALE.C
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_LOCALE.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler doesn't know the locale classic
//| CGAL_CFG_NO_LOCALE is set. 

#include <locale>


int main()
{
  return std::isspace(' ', std::locale::classic()) ? 0 : 1;
}


