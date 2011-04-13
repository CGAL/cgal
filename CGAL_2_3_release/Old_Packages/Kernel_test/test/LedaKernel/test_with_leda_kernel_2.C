// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : test_with_leda_kernel_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifdef CGAL_USE_LEDA
#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/rat_leda.h>
#include <CGAL/rat_leda_in_CGAL_2.h>
#include <CGAL/predicates_on_points_rat_leda_2.h>
#include <CGAL/_test_fct_point_2.C>

int
main()
{
  _test_fct_point_2( CGAL::use_rat_leda_kernel() );
  return 0;
}
#else
int main() { return 0; }
#endif // CGAL_USE_LEDA
