// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// release       : 
// release_date  : 
//
// file          : 
// source        : chull_traits.lw
// revision      : 2.3  
// revision_date : 01 Feb 2000
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================

#include <CGAL/Homogeneous.h>
#include <CGAL/leda_integer.h>
#include <CGAL/_test_cls_chull_traits_3.C>

typedef CGAL::Homogeneous<leda_integer>    RepCls;

int
main()
{
  __test_cls_chull_traits_3( RepCls() );
  return 0;
}
