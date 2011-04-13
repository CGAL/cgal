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
// file          : test/KernelCartesian/test_new_kernelC3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Precise_numbers.h>
#include <CGAL/_test_new_3.h>

typedef CGAL::Cartesian<Precise_rational>  Kernel;

int main()
{
  test_new_3( Kernel() );
  return 0;
}

