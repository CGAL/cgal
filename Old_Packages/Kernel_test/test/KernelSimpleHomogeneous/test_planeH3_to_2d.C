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
// file          : test_planeH3_to_2d.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/basic.h>
#include <cassert>

#ifdef CGAL_USE_GMP
# include <CGAL/Gmpz.h>
typedef CGAL::Gmpz    Precise_integer;
#else
# ifdef CGAL_USE_LEDA
#  include <CGAL/leda_integer.h>
typedef leda_integer  Precise_integer;
# endif // CGAL_USE_LEDA
#endif // CGAL_USE_GMP

#include <CGAL/Simple_homogeneous.h>
#include <CGAL/_test_mf_plane_3_to_2d.C>

int
main()
{
  typedef   CGAL::Simple_homogeneous< Precise_integer >                 H_Cls;
  std::cout << "Testing with Simple_homogeneous<Precise_integer> :";
  std::cout << std::endl;
  _test_mf_plane_3_to_2d( H_Cls() );
  return 0;
}
