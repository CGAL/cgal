// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 2000, August 11
//
// file          : include/CGAL/ExternalMemoryStructures/R_tree_internal_storage.h
// package       : ExternalMemoryStructures (0.631)
// maintainer    : Philipp Kramer <kramer@inf.ethz.ch>
// chapter       : $CGAL_Chapter: Basic / External Data Structures $
// source        : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Gabriele Neyer<neyer@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Peter Widmayer <widmayer@inf.ethz.ch>)
//
// Implementation of the R_tree_internal_storage.
// ======================================================================
#ifndef __R_tree_internal_storage_H__
#define __R_tree_internal_storage_H__
#include<CGAL/R_tree_internal_db.h>
#include<CGAL/R_tree.h>

CGAL_BEGIN_NAMESPACE

class R_tree_internal_storage{
 public:
  typedef CGAL::R_tree_internal_db IO_tree_traits_nodes;
  typedef CGAL::R_tree_internal_db IO_tree_traits_leaves;
  const static  int IO_min_cap_nodes=2;
  const static  int IO_max_cap_nodes=4;
  const static  int IO_min_cap_leaves=2;
  const static  int IO_page_size=70;
  const static bool headerextern=false;
};

CGAL_END_NAMESPACE
#endif
