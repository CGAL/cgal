//
//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// This file is part of the Generic Graph Component Library
//
// You should have received a copy of the License Agreement for the
// Generic Graph Component Library along with the software;  see the
// file LICENSE.  If not, contact Office of Research, University of Notre
// Dame, Notre Dame, IN  46556.
//
// Permission to modify the code and to distribute modified code is
// granted, provided the text of this NOTICE is retained, a notice that
// the code was modified is included with the above COPYRIGHT NOTICE and
// with the COPYRIGHT NOTICE in the LICENSE file, and that the LICENSE
// file is distributed with the modified code.
//
// LICENSOR MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
// By way of example, but not limitation, Licensor MAKES NO
// REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY
// PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS
// OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS
// OR OTHER RIGHTS.
//=======================================================================
//
#ifndef BOOST_GRAPH_DETAIL_MUTABLE_HEAP_H
#define BOOST_GRAPH_DETAIL_MUTABLE_HEAP_H

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

namespace boost {

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

}
#endif
