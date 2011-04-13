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

#include <CGAL/Triangulation_3.h>
#include <CGAL/_test_cls_triangulation_3.C>

bool del = false;

int main()
{
  typedef Test_rep_cartesian                                        traits;
  // Using vertex_base_pointer induces a memory leak (not a bug, but the test
  // program is not adapted), so we use the normal vertex.
  // typedef CGAL::Triangulation_vertex_base_pointer_3<traits>        Vb;
  typedef CGAL::Triangulation_3<traits>                         Cls3;

  _test_cls_triangulation_3( Cls3() );

  return 0;
}
