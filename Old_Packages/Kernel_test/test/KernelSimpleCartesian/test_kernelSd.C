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
// file          : test_kernelSd.C
// revision      : 3.8
// revision_date : 08 Oct 2000 
// author(s)     : Stefan Schirra
//
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de> 
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <iostream>
#include <strstream>
#include <cassert>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_d.h>
#include <CGAL/_test_cls_point_d.h>

typedef CGAL::Simple_cartesian<double>   RC;

int
main()
{
  _test_cls_point_d( RC() );
  return 0;
}
