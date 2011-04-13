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
// file          : test/NewKernel/test_new_kernelC3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/config.h>
#undef CGAL_CFG_NO_ADVANCED_KERNEL
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
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


typedef CGAL::Cartesian<Precise_rational>  Kernel;

int main()
{
  test_new_3( Kernel() );
  return 0;
}

