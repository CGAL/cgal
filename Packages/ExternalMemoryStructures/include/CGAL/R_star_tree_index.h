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
// file          : include/CGAL/ExternalMemoryStructures/R_star_tree_index.h
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
// Implementation of the R_star_tree_index.
// ======================================================================
#ifndef __R_star_tree_index_H__
#define __R_star_tree_index_H__
#include<CGAL/R_tree_index_implementations.h>
#include<CGAL/R_tree.h>
#include<CGAL/R_tree_key.h>

CGAL_BEGIN_NAMESPACE

template<class R_tree_traits>
class R_star_tree_index{
 public:
  typedef CGAL::star_split_node<CGAL::R_tree_value<R_tree_traits>, 
    R_tree_traits,
    CGAL::sort_axis_key_2_dim<CGAL::R_tree_value<R_tree_traits> > > 
    Split_node;
  typedef CGAL::star_split_leaf<CGAL::R_tree_leaf_data<
    R_tree_traits>, R_tree_traits,
    CGAL::sort_axis_key_2_dim<CGAL::R_tree_leaf_data<R_tree_traits> > > 
    Split_leaf;
  typedef CGAL::star_choose_subtree<CGAL::R_tree_value<
    R_tree_traits>, 
    R_tree_traits> Choose_subtree;
  typedef  CGAL::reinsertion_node<CGAL::R_tree_value<
    R_tree_traits>, 
    R_tree_traits> Reinsertion_node;
  typedef CGAL::reinsertion_leaf<CGAL::R_tree_leaf_data<
    R_tree_traits>, R_tree_traits> 
    Reinsertion_leaf;
  const static bool reinsertions=true;
};

CGAL_END_NAMESPACE
#endif
