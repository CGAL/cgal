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
// source        : NewKernel/test_programs.fw
// file          : test/KernelHomogeneous/test_new_kernelH3.C
// revision      : 3.8
// revision_date : 08 Oct 2000 
// author(s)     : Stefan Schirra
//
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de> 
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
typedef leda_integer                  Precise_integer;
typedef leda_rational                 Precise_rational;
#else
#  ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#include <CGAL/Quotient.h>
typedef CGAL::Gmpz                    Precise_integer;
typedef CGAL::Quotient<CGAL::Gmpz>    Precise_rational;
#  endif // CGAL_USE_GMP
#endif // CGAL_USE_LEDA



#include <CGAL/_test_new_3.h>


typedef CGAL::Homogeneous<Precise_integer> Kernel;

int main()
{
  test_new_3( Kernel() );
  return 0;
}

