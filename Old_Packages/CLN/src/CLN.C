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
// release       :
// release_date  :
// 
// file          : src/CLN.C
// revision      : $Revision$
// revision_date : $Date$
// package       : CLN
// author(s)     : Sylvain Pion
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================
 
#ifdef CGAL_USE_CLN

#include <CGAL/basic.h>
#include <cl_output.h>

CGAL_BEGIN_NAMESPACE

// It's a workaround for a bug that happens on Solaris 2.6 with gcc 2.95,
// and libcln.so (not .a).
// It doesn't happen on Linux with gcc 2.95.

// Namely, the default base for printing should be 10, but it's not
// initialized as it should for some reason...

// So we make a static object that we initialize here instead.

struct workaround_4_CLN
{
  workaround_4_CLN()
  {
    cl_default_print_flags.rational_base = 10;
  }
};

static workaround_4_CLN w;

CGAL_END_NAMESPACE

#endif // CGAL_USE_CLN
