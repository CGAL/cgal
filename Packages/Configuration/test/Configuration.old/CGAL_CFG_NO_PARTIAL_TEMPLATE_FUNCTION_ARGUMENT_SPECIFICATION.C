// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : 
// release_date  :
//
// file          : config/testfiles/CGAL_CFG_NO_PARTIAL_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ============================================================================

// CGAL_CFG_NO_PARTIAL_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler whether it likes
// partial specification of template arguments in template function calls.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler doesn't like explicit partial specification of 
//| template arguments in template function calls, the flag
//| CGAL_CFG_NO_PARTIAL_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION is set.

template <class To, class From>
To
convert(const From& f)
{ return (To)(f); }

int
main()
{
  double d; 
  int    i = 1;
  d = convert<double>(i);
  return 0;
}

// EOF //
