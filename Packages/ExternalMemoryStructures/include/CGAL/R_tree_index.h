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
// file          : include/CGAL/ExternalMemoryStructures/R_tree_index.h
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
// Implementation of the R_tree_index.
// ======================================================================
#ifndef __R_tree_index_H__
#define __R_tree_index_H__
#include<CGAL/R_tree_index_implementations.h>
#include<CGAL/R_tree.h>

CGAL_BEGIN_NAMESPACE

template<class R_tree_traits>
class R_tree_index{
 public:
  typedef CGAL::quadratic_split_node<CGAL::R_tree_value<
    R_tree_traits>, 
    R_tree_traits> Split_node;
  typedef CGAL::quadratic_split_leaf<CGAL::R_tree_leaf_data
    <R_tree_traits>, R_tree_traits> Split_leaf;
  typedef CGAL::choose_subtree<CGAL::R_tree_value<R_tree_traits>, 
    R_tree_traits> Choose_subtree;
  typedef  CGAL::dummy_reinsertion_node<CGAL::R_tree_value<R_tree_traits>, 
    R_tree_traits> Reinsertion_node;
  typedef CGAL::dummy_reinsertion_leaf<CGAL::R_tree_leaf_data<
    R_tree_traits>, R_tree_traits> 
    Reinsertion_leaf;
  const static bool reinsertions=false;
};

CGAL_END_NAMESPACE
#endif
