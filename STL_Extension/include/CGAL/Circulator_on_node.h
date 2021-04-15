// Copyright (c) 2003
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion

#ifndef CGAL_CIRCULATOR_ON_NODE_H
#define CGAL_CIRCULATOR_ON_NODE_H 1

#include <CGAL/circulator.h>

namespace CGAL {

template < class Node,
           class Next,
           class Prev,
           class Ref = Node&,
           class Ptr = Node*,
           class Ctg = Bidirectional_circulator_tag>
class Circulator_on_node {
protected:
  Ptr      nt;    // The internal node ptr.
public:
  typedef  Circulator_on_node<Node,Next,Prev,Ref,Ptr,Ctg> Self;

  typedef  Ctg                iterator_category;
  typedef  Node               value_type;
  typedef  std::ptrdiff_t     difference_type;
  typedef  std::size_t        size_type;
  typedef  Ref                reference;
  typedef  Ptr                pointer;

  // CREATION
  // --------

  Circulator_on_node() : nt(0) {}
  Circulator_on_node( Ptr p) : nt(p) {}

  // OPERATIONS Forward Category
  // ---------------------------

  Ptr  ptr() const { return nt;}

  bool operator==( std::nullptr_t p) const {
    CGAL_assertion( p == 0);
    CGAL_USE(p);
    return ( nt == 0);
  }
  bool  operator!=( std::nullptr_t p) const { return !(*this == p); }
  bool  operator==( const Self& i) const { return ( nt == i.nt); }
  bool  operator!=( const Self& i) const { return !(*this == i); }
  Ref   operator*()  const { return *nt; }
  Ptr   operator->() const { return nt; }
  Self& operator++() {
    Next next;
    nt = next(nt);
    return *this;
  }
  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  // OPERATIONS Bidirectional Category
  // ---------------------------------

  Self& operator--() {
    Prev prev;
    nt = prev(nt);
    return *this;
  }
  Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
  }
};

} //namespace CGAL
#endif // CGAL_CIRCULATOR_ON_NODE_H //
// EOF //
