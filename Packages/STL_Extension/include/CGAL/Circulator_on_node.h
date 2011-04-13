// ============================================================================
//
// Copyright (c) 1997, 1998, 1999, 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : Circulator_on_node.h
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
// package       : $CGAL_Package: STL_Extension $
// source        : stl_extension.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@cs.unc.edu>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// An circulator built over a linked node structure.
// ============================================================================

#ifndef CGAL_CIRCULATOR_ON_NODE_H
#define CGAL_CIRCULATOR_ON_NODE_H 1
#include <CGAL/circulator.h>

CGAL_BEGIN_NAMESPACE

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

  bool operator==( CGAL_NULL_TYPE p) const {
    CGAL_assertion( p == 0);
    return ( nt == 0);
  }
  bool  operator!=( CGAL_NULL_TYPE p) const { return !(*this == p); }
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

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
#ifndef CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
  friend inline  iterator_category
  iterator_category( const Self&) { return iterator_category(); }
  friend inline  value_type*
  value_type( const Self&) { return (value_type*)(0); }
  friend inline  difference_type*
  distance_type( const Self&) { return (difference_type*)(0); }
  friend inline  Circulator_tag
  query_circulator_or_iterator( const Self&) { return Circulator_tag(); }
#endif // CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};

CGAL_END_NAMESPACE
#endif // CGAL_CIRCULATOR_ON_NODE_H //
// EOF //
