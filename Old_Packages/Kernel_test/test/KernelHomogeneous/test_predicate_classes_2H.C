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
// source        :
// file          : test_predicate_classes_2H.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Precise_numbers.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/_test_fctobj_point_line_2.C>

int
main()
{
  typedef   CGAL::Homogeneous<Precise_integer>                H_Cls;
  typedef   CGAL::Point_2<H_Cls>  H_Point;
  typedef   CGAL::Line_2<H_Cls>   H_Line;
  std::cout << "Testing 2d with Homogeneous<Precise_integer> :  ";
  std::cout << std::endl;
  _test_fctobj_point_line_2( H_Point(), H_Line());

  return 0;
}
