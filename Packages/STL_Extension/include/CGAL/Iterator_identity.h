#line 97 "stl_extension.aw"

#line 20 "cgal_header.awi"
// ============================================================================
//
// Copyright (c) 1997, 1998, 1999 The CGAL Consortium
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
// file          : Iterator_identity.h
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
// package       : $CGAL_Package: STL_Extension $
// source        : stl_extension.fw
#line 37 "cgal_header.awi"
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// An iterator adaptor for the identity function.
// ============================================================================

#line 101 "stl_extension.aw"
#ifndef CGAL_ITERATOR_IDENTITY_H
#define CGAL_ITERATOR_IDENTITY_H 1
#line 1990 "stl_extension.aw"
#ifndef CGAL_CIRCULATOR_H
#include <CGAL/circulator.h>
#endif

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
template < class I, class Ref, class Ptr, class Val, class Dist, class Ctg>
#else
template < class I,
           class Ref  = typename std::iterator_traits<I>::reference,
           class Ptr  = typename std::iterator_traits<I>::pointer,
           class Val  = typename std::iterator_traits<I>::value_type,
           class Dist = typename std::iterator_traits<I>::difference_type,
           class Ctg = typename std::iterator_traits<I>::iterator_category>
#endif
class Iterator_identity {
private:
  I        nt;    // The internal iterator.
public:
  typedef I      Iterator;
  typedef Iterator_identity<I,Ref,Ptr,Val,Dist,Ctg>   Self;
  typedef Ctg    iterator_category;
  typedef Val    value_type;
  typedef Dist   difference_type;
  typedef Ref    reference;
  typedef Ptr    pointer;

  // CREATION
  // --------

  Iterator_identity() {}
  Iterator_identity( Iterator j) : nt(j) {}

  // OPERATIONS Forward Category
  // ---------------------------

  Iterator  current_iterator() const { return nt;}

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
  bool operator<( const Self& i) const {
    return ( nt < i.nt);                                     //###//
  }
  bool operator>( const Self& i) const {
    return i < *this;
  }
  bool operator<=( const Self& i) const {
    return !(i < *this);
  }
  bool operator>=( const Self& i) const {
    return !(*this < i);
  }
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  friend inline  value_type*
  value_type( const Self&) {
    return (value_type*)(0);
  }
  friend inline  Ctg
  iterator_category( const Self&){
    return Ctg();
  }
  friend inline  Dist*
  distance_type( const Self&) {
    return (Dist*)(0);
  }
  friend inline  Iterator_tag
  query_circulator_or_iterator( const Self&) {
    return Iterator_tag();
  }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};

template < class D, class I, class Ref, class Ptr, class Val,
           class Dist, class Ctg>
inline
Iterator_identity<I,Ref,Ptr,Val,Dist,Ctg>
operator+( D n, Iterator_identity<I,Ref,Ptr,Val,Dist,Ctg> i) {
  return i += Dist(n);
}

CGAL_END_NAMESPACE
#line 104 "stl_extension.aw"
#endif // CGAL_ITERATOR_IDENTITY_H //
// EOF //
