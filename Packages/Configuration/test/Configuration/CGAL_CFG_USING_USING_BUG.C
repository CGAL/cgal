// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.2-I-20 $
// release_date  : $CGAL_Date: 2000/06/02 $
//
// file          : config/testfiles/CGAL_CFG_USING_USING_BUG.C
// package       : Configuration (2.3)
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_USING_USING_BUG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler does not accept a using declaration referring to a 
//| symbol that is again declared by a using declaration, the flag
//| CGAL_CFG_USING_USING_BUG is set.

namespace L { int foo() { return 0; } }
namespace M { using L::foo; }
namespace N { using M::foo; }

int main()
{
  return N::foo();
}

