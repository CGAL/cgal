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

//| If a compiler doesn't know the limits (g++-2.95)
//| or has a bug in the implementation (Sun CC 5.4)
//| CGAL_CFG_NO_LIMITS is set. 

#include <limits>


int main()
{
  if(std::numeric_limits<double>::denorm_min() == 0){
    return 1;
  }
  return 0;
}


