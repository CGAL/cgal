//
//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
//
// $URL$
// $Id$
//
// NOTE: this file have been taken from boost 1.46.1 for using
//       with Modificable_priority_queue (to enhance the 
//       non-documented mutable_queue).
//       original file is <boost/pending/mutable_heap.hpp>
//
#ifndef CGAL_INTERNAL_BOOST_MUTABLE_HEAP_H
#define CGAL_INTERNAL_BOOST_MUTABLE_HEAP_H

/*
  There are a few things wrong with this set of functions.

  ExternalData should be removed, it is not part of the core
  algorithm. It can be handled inside the tree nodes.

  The swap() should be replaced by assignment since its use is causing
  the number of memory references to double.

  The min_element should be replaced by a fixed length loop
  (fixed at d for d-heaps).

  The member functions of TreeNode should be changed to global
  functions.

  These functions will be replaced by those in heap_tree.h

 */

namespace CGAL{ namespace internal {
  
namespace boost_ {

  template <class TreeNode, class Compare, class ExternalData>
  inline TreeNode up_heap(TreeNode x, const Compare& comp, ExternalData& edata) {
    while (x.has_parent() && comp(x, x.parent()))
      x.swap(x.parent(), edata);
    return x;
  }

  template <class TreeNode, class Compare, class ExternalData>
  inline TreeNode down_heap(TreeNode x, const Compare& comp, ExternalData& edata) {
    while (x.children().size() > 0) {
      typename TreeNode::children_type::iterator 
        child_iter = std::min_element(x.children().begin(),
                                      x.children().end(), 
                                      comp);
      if (comp(*child_iter, x))
        x.swap(*child_iter, edata);
      else
        break;
    }
    return x;
  }

  template <class TreeNode, class Compare, class ExternalData>
  inline void update_heap(TreeNode x, const Compare& comp, ExternalData& edata) {
    x = down_heap(x, comp, edata);
    (void)up_heap(x, comp, edata);
  }

} } } //namespace CGAL::internal::boost_
#endif //CGAL_INTERNAL_BOOST_MUTABLE_HEAP_H
