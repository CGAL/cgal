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
// source        : test_kernel_programs.fw
// file          : test_kernelH3.C
// revision      : 3.8
// revision_date : 08 Oct 2000 
// author(s)     : Stefan Schirra
//
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de> 
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

#include <CGAL/Homogeneous.h>
#include <CGAL/_test_3.C>

int
main()
{
  typedef   CGAL::Homogeneous<Precise_integer>     Cls;
  std::cout << "Testing 3d with Homogeneous<Precise_integer> :" << std::endl;
  _test_3( Cls() );
  return 0;
}
