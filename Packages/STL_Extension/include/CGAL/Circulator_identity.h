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
// file          : Circulator_identity.h
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
// package       : $CGAL_Package: STL_Extension $
// source        : stl_extension.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@cs.unc.edu>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// An circulator adaptor for the identity function.
// ============================================================================

#ifndef CGAL_CIRCULATOR_IDENTITY_H
#define CGAL_CIRCULATOR_IDENTITY_H 1
#include <CGAL/circulator.h>

CGAL_BEGIN_NAMESPACE

template < class C,
           class Ref = typename C::reference,
           class Ptr = typename C::pointer>
class Circulator_identity {
protected:
  C        nt;    // The internal circulator.
public:
  typedef C  Circulator;
  typedef Circulator_identity<C,Ref,Ptr> Self;

  typedef typename  C::iterator_category   iterator_category;
  typedef typename  C::value_type          value_type;
  typedef typename  C::difference_type     difference_type;
  typedef typename  C::size_type           size_type;
  typedef typename  C::reference           reference;
  typedef typename  C::pointer             pointer;

  // CREATION
  // --------

  Circulator_identity() {}
  Circulator_identity( Circulator j) : nt(j) {}

  // OPERATIONS Forward Category
  // ---------------------------

  Circulator  current_circulator() const { return nt;}

  bool operator==( CGAL_NULL_TYPE p) const {
    CGAL_assertion( p == 0);
    return ( nt == 0);                                    //###//
  }
  bool operator!=( CGAL_NULL_TYPE p) const {
    return !(*this == p);
  }
  bool operator==( const Self& i) const {
    return ( nt == i.nt);                                    //###//
  }
  bool operator!=( const Self& i) const {
    return !(*this == i);
  }
  Ref  operator*() const {
    return *nt;                                              //###//
  }
  Ptr  operator->() const {
    return nt.operator->();                                  //###//
  }
  Self& operator++() {
    ++nt;                                                    //###//
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
    --nt;                                                    //###//
    return *this;
  }
  Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
  }

  // OPERATIONS Random Access Category
  // ---------------------------------

  Self  min_circulator() const {
    return Self( nt.min_circulator());                       //###//
  }
  Self& operator+=( difference_type n) {
    nt += n;                                                 //###//
    return *this;
  }
  Self  operator+( difference_type n) const {
    Self tmp = *this;
    return tmp += n;
  }
  Self& operator-=( difference_type n) {
    return operator+=( -n);
  }
  Self  operator-( difference_type n) const {
    Self tmp = *this;
    return tmp += -n;
  }
  difference_type  operator-( const Self& i) const {
    return nt - i.nt;                                        //###//
  }
  Ref  operator[]( difference_type n) const {
    Self tmp = *this;
    tmp += n;
    return tmp.operator*();
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

template < class Dist, class C, class Ref, class Ptr>
inline
Circulator_identity<C,Ref,Ptr>
operator+( Dist n, Circulator_identity<C,Ref,Ptr> i) {
  return i += n;
}

CGAL_END_NAMESPACE
#endif // CGAL_CIRCULATOR_IDENTITY_H //
// EOF //
