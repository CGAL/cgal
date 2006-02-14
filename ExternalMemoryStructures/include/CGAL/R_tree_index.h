// Copyright (c) 1998  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Gabriele Neyer<neyer@inf.ethz.ch>
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
