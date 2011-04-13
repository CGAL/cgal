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
// file          : test/KernelSimpleCartesian/test_new_kernelS2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Precise_numbers.h>
#include <CGAL/_test_new_2.h>

typedef CGAL::Simple_cartesian<Precise_rational>  Kernel;

int main()
{
  test_new_2( Kernel() );
  return 0;
}
