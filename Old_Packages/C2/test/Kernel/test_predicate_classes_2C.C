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
// file          : test_predicate_classes_2C.C
// revision      : 2.0.5
// revision_date : 24 Mar 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/leda_real.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/_test_fctobj_point_line_2.C>

int
main()
{
  typedef   CGAL::Cartesian<leda_real>     C_Cls;
  typedef   CGAL::Point_2<C_Cls>  C_Point;
  typedef   CGAL::Line_2<C_Cls>   C_Line;
  cout << "Testing 2d with Cartesian<leda_real> > :" << endl;
  _test_fctobj_point_line_2( C_Point(), C_Line());

  return 0;
}
