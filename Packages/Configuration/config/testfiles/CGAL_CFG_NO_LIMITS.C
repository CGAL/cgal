// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_LIMITS.C
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_LIMITS.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler doesn't know <limits> (g++-2.95)
//| or has a bug in the implementation (Sun CC 5.4, MipsPro CC)
//| CGAL_CFG_NO_LIMITS is set. 

#include <limits>

int main()
{
  double d = std::numeric_limits<double>::denorm_min();
  double e = std::numeric_limits<double>::min();
  if (d == 0 || d == e)
    return 1;
  return 0;
}
