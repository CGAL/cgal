// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// file          : test_triangulation_3.C
// revision      : 
// revision_date : 
// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/basic.h>
#include <cassert>

#include <list>
#include <vector>

#include <CGAL/_test_types.h>
#include <CGAL/triple.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/_test_cls_delaunay_3.C>

bool del=true;

int main()
{
  typedef Test_rep_cartesian                                        traits;
  typedef CGAL::Delaunay_triangulation_3<traits>                Cls;

  _test_cls_delaunay_3( Cls() );

  return 0;
}
