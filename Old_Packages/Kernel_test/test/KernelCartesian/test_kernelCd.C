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
// file          : test_kernelCd.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <iostream>
#include <strstream>
#include <cassert>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_d.h>
#include <CGAL/_test_cls_point_d.h>

typedef CGAL::Cartesian<double>   RC;

int
main()
{
  _test_cls_point_d( RC() );
  return 0;
}
