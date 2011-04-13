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
// file          : Random_access_adaptor.h
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
// Random Access Adaptor provides random access for sequences.
// ============================================================================

#ifndef CGAL_RANDOM_ACCESS_ADAPTOR_H
#define CGAL_RANDOM_ACCESS_ADAPTOR_H 1
#include <vector>
#include <CGAL/circulator.h>

CGAL_BEGIN_NAMESPACE

template < class IC>
class Random_access_adaptor {

  // DEFINITION
  //
  // The class Random_access_adaptor<IC> provides a random access
  // for data structures. Either the data structure supports random access
  // iterators or circulators where this class maps function calls to the
  // iterator or circulator, or a STL `vector' is used to provide the random
  // access. The iterator or circulator of the data structure are of type
  // `IC'.
  //
  // CREATION

protected:
  typedef std::vector< IC> Index;
  Index   index;
  IC      start;

public:
  typedef typename Index::size_type  size_type;

  void init_index( IC i, const IC& j, std::forward_iterator_tag);
  void init_index( const IC& i, const IC& j,
                   std::bidirectional_iterator_tag){
    init_index( i, j, std::forward_iterator_tag());
  }
  void init_index( const IC& i, const IC&,
                   std::random_access_iterator_tag){
    start = i;
  }
  void init_index( const IC& i, const IC& j) {
#if !defined(CGAL_CFG_NO_ITERATOR_TRAITS) || \
    defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
    typedef typename std::iterator_traits<IC>::iterator_category ICC;
    init_index( i, j, ICC());
#else
    init_index( i, j, std::iterator_category( i));
#endif
  }


  void reserve( size_type r, std::forward_iterator_tag) {
    index.reserve( r);
  }
  void reserve( size_type r, std::bidirectional_iterator_tag){
    reserve( r, std::forward_iterator_tag());
  }
  void reserve( size_type, std::random_access_iterator_tag){}


  void push_back( const IC& k, std::forward_iterator_tag) {
    index.push_back(k);
  }
  void push_back( const IC& k, std::bidirectional_iterator_tag){
    push_back( k, std::forward_iterator_tag());
  }
  void push_back( const IC&, std::random_access_iterator_tag){}


  const IC& find( size_type n, std::forward_iterator_tag) const {
    // returns inverse index of k.
    CGAL_assertion( n < index.size());
    return index[n];
  }
  const IC& find( size_type n, std::bidirectional_iterator_tag) const {
    return find( n, std::forward_iterator_tag());
  }
  IC  find( size_type n, std::random_access_iterator_tag) const {
    return start + n;
  }

  typedef IC   iterator;
  typedef IC   Circulator;

  Random_access_adaptor() : start(IC()) {}
  // invalid index.

  Random_access_adaptor( const IC& i) : start(i) {}
  // empty random access index initialized to start at i.

  Random_access_adaptor( const IC& i, const IC& j) : start(i) {
    // random access index initialized with range [i,j).
    init_index( i, j);
  }

  void reserve( size_type r) {
    // reserve r entries, if a `vector' is used internally.
#if !defined(CGAL_CFG_NO_ITERATOR_TRAITS) || \
    defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
    typedef typename std::iterator_traits<IC>::iterator_category ICC;
    reserve( r, ICC());
#else
    reserve( r, std::iterator_category( IC()));
#endif
  }

  // OPERATIONS

  IC  find( size_type n) const {
    // returns inverse index of k.
#if !defined(CGAL_CFG_NO_ITERATOR_TRAITS) || \
    defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
    typedef typename std::iterator_traits<IC>::iterator_category ICC;
    return find( n, ICC());
#else
    return find( n, std::iterator_category( IC()));
#endif
  }

  IC  operator[]( size_type n) const { return find(n); }

  void push_back( const IC& k) {
    // adds k at the end of the indices.
#if !defined(CGAL_CFG_NO_ITERATOR_TRAITS) || \
    defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
    typedef typename std::iterator_traits<IC>::iterator_category ICC;
    push_back( k, ICC());
#else
    push_back( k, std::iterator_category( k));
#endif
  }
};

template < class IC>
void
Random_access_adaptor< IC>::init_index( IC i, const IC& j,
                                        std::forward_iterator_tag) {
  if ( ! is_empty_range( i, j)) {
    do {
      index.push_back( i);
    } while ((++i) != (j));
  }
}

CGAL_END_NAMESPACE
#endif // CGAL_RANDOM_ACCESS_ADAPTOR_H //
// EOF //
