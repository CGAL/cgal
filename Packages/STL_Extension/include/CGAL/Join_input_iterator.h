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
// file          : Join_input_iterator.h
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
// Join_input_iterator.
// ============================================================================

#ifndef CGAL_JOIN_INPUT_ITERATOR_H
#define CGAL_JOIN_INPUT_ITERATOR_H 1
#include <CGAL/circulator.h>

CGAL_BEGIN_NAMESPACE

template < class I1, class  Creator >
class Join_input_iterator_1 {
  // the join of one iterator `i1'. Applies `Creator' with
  // one argument `*i1'. `value_type' is equal to
  // `Creator::result_type'.
public:
  typedef Join_input_iterator_1<I1,Creator>  Self;
  typedef std::input_iterator_tag            iterator_category;
  typedef typename Creator::result_type      value_type;
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  typedef std::ptrdiff_t                     difference_type;
#else
  typedef std::iterator_traits<I1>           ITraits;
  typedef typename ITraits::difference_type  difference_type;
#endif
  typedef  const value_type&                 reference;
  typedef  const value_type*                 pointer;

protected:
  I1          j1;    // The 1st internal iterator.
  value_type  val;   // The current (internal) value.

public:
  // CREATION
  // --------

  Join_input_iterator_1() {}
  Join_input_iterator_1( I1 i1) : j1(i1), val(Creator()(*j1)) {}

  // OPERATIONS Forward Category
  // ---------------------------

  I1  current_iterator1() const { return j1;}

  bool operator==( const Self& i) const { return ( j1 == i.j1); }
  bool operator!=( const Self& i) const { return !(*this == i); }
  reference  operator*()    const { return val; }
  pointer    operator->()   const { return &val; }
  Self& operator++() {
    ++j1;
    val = Creator()(*j1);
    return *this;
  }
  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
#ifndef CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
  friend inline  value_type*
  value_type( const Self&) {
    return (value_type*)(0);
  }
  friend inline  std::input_iterator_tag
  iterator_category( const Self&){
    return std::input_iterator_tag();
  }
  friend inline  difference_type*
  distance_type( const Self&) {
    return (difference_type*)(0);
  }
  friend inline  Iterator_tag
  query_circulator_or_iterator( const Self&) {
    return Iterator_tag();
  }
#endif // CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};

template < class I1, class I2, class  Creator >
class Join_input_iterator_2 {
  // the join of two iterators `i1' and `i2'. Applies `Creator' with
  // two arguments `*i1' and `*i2'. `value_type' is equal to
  // `Creator::result_type'.
public:
  typedef Join_input_iterator_2<I1,I2,Creator> Self;

  typedef std::input_iterator_tag              iterator_category;
  typedef typename Creator::result_type        value_type;
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  typedef std::ptrdiff_t                       difference_type;
#else
  typedef std::iterator_traits<I1>             ITraits;
  typedef typename ITraits::difference_type    difference_type;
#endif
  typedef const value_type&                    reference;
  typedef const value_type*                    pointer;

protected:
  I1          j1;    // The 1st internal iterator.
  I2          j2;    // The 2nd internal iterator.
  value_type  val;   // The current (internal) value.

public:
  // CREATION
  // --------

  Join_input_iterator_2() {}
  Join_input_iterator_2( I1 i1, I2 i2)
  : j1(i1), j2(i2), val(Creator()(*j1,*j2)) {}

  // OPERATIONS Forward Category
  // ---------------------------

  I1  current_iterator1() const { return j1;}
  I2  current_iterator2() const { return j2;}

  bool operator==( const Self& i) const {
    return ( j1 == i.j1 && j2 == i.j2);
  }
  bool operator!=( const Self& i) const { return !(*this == i); }
  reference operator*() const { return val; }
  pointer   operator->() const { return &val; }
  Self& operator++() {
    ++j1;
    ++j2;
    val = Creator()(*j1,*j2);
    return *this;
  }
  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
#ifndef CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
  friend inline  value_type*
  value_type( const Self&) {
    return (value_type*)(0);
  }
  friend inline  std::input_iterator_tag
  iterator_category( const Self&){
    return std::input_iterator_tag();
  }
  friend inline  difference_type*
  distance_type( const Self&) {
    return (difference_type*)(0);
  }
  friend inline  Iterator_tag
  query_circulator_or_iterator( const Self&) {
    return Iterator_tag();
  }
#endif // CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};

template < class I1, class I2, class I3, class  Creator >
class Join_input_iterator_3 {
  // the join of two iterators `i1' up to `i3'. Applies `Creator' with
  // three arguments `*i1' up to `*i3'. `value_type' is equal to
  // `Creator::result_type'.
public:
  typedef Join_input_iterator_3<I1,I2,I3,Creator> Self;

  typedef std::input_iterator_tag                 iterator_category;
  typedef typename Creator::result_type           value_type;
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  typedef std::ptrdiff_t                          difference_type;
#else
  typedef std::iterator_traits<I1>                ITraits;
  typedef typename ITraits::difference_type       difference_type;
#endif
  typedef const value_type&                       reference;
  typedef const value_type*                       pointer;

protected:
  I1          j1;    // The 1st internal iterator.
  I2          j2;    // The 2nd internal iterator.
  I3          j3;    // The 3rd internal iterator.
  value_type  val;   // The current (internal) value.

public:
  // CREATION
  // --------

  Join_input_iterator_3() {}
  Join_input_iterator_3( I1 i1, I2 i2, I3 i3)
  : j1(i1), j2(i2), j3(i3), val(Creator()(*j1,*j2,*j3)) {}

  // OPERATIONS Forward Category
  // ---------------------------

  I1  current_iterator1() const { return j1;}
  I2  current_iterator2() const { return j2;}
  I3  current_iterator3() const { return j3;}

  bool operator==( const Self& i) const {
    return ( j1 == i.j1 && j2 == i.j2 && j3 == i.j3);
  }
  bool operator!=( const Self& i) const { return !(*this == i); }
  reference operator*() const { return val; }
  pointer   operator->() const { return &val; }
  Self& operator++() {
    ++j1;
    ++j2;
    ++j3;
    val = Creator()(*j1,*j2,*j3);
    return *this;
  }
  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
#ifndef CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
  friend inline  value_type*
  value_type( const Self&) {
    return (value_type*)(0);
  }
  friend inline  std::input_iterator_tag
  iterator_category( const Self&){
    return std::input_iterator_tag();
  }
  friend inline  difference_type*
  distance_type( const Self&) {
    return (difference_type*)(0);
  }
  friend inline  Iterator_tag
  query_circulator_or_iterator( const Self&) {
    return Iterator_tag();
  }
#endif // CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};

template < class I1, class I2, class I3, class I4, class  Creator >
class Join_input_iterator_4 {
  // the join of two iterators `i1' up to `i4'. Applies `Creator' with
  // four arguments `*i1' up to `*i4'. `value_type' is equal to
  // `Creator::result_type'.
public:
  typedef Join_input_iterator_4<I1,I2,I3,I4,Creator> Self;

  typedef std::input_iterator_tag             iterator_category;
  typedef typename Creator::result_type       value_type;
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  typedef std::ptrdiff_t                      difference_type;
#else
  typedef std::iterator_traits<I1>            ITraits;
  typedef typename ITraits::difference_type   difference_type;
#endif
  typedef const value_type&                   reference;
  typedef const value_type*                   pointer;

protected:
  I1          j1;    // The 1st internal iterator.
  I2          j2;    // The 2nd internal iterator.
  I3          j3;    // The 3rd internal iterator.
  I4          j4;    // The 4th internal iterator.
  value_type  val;   // The current (internal) value.

public:
  // CREATION
  // --------

  Join_input_iterator_4() {}
  Join_input_iterator_4( I1 i1, I2 i2, I3 i3, I4 i4)
  : j1(i1), j2(i2), j3(i3), j4(i4), val(Creator()(*j1,*j2,*j3,*j4)){}

  // OPERATIONS Forward Category
  // ---------------------------

  I1  current_iterator1() const { return j1;}
  I2  current_iterator2() const { return j2;}
  I3  current_iterator3() const { return j3;}
  I4  current_iterator4() const { return j4;}

  bool operator==( const Self& i) const {
    return ( j1 == i.j1 &&
             j2 == i.j2 &&
             j3 == i.j3 &&
             j4 == i.j4);
  }
  bool operator!=( const Self& i) const { return !(*this == i); }
  reference operator*() const { return val; }
  pointer   operator->() const { return &val; }
  Self& operator++() {
    ++j1;
    ++j2;
    ++j3;
    ++j4;
    val = Creator()(*j1,*j2,*j3,*j4);
    return *this;
  }
  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
#ifndef CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
  friend inline  value_type*
  value_type( const Self&) {
    return (value_type*)(0);
  }
  friend inline  std::input_iterator_tag
  iterator_category( const Self&){
    return std::input_iterator_tag();
  }
  friend inline  difference_type*
  distance_type( const Self&) {
    return (difference_type*)(0);
  }
  friend inline  Iterator_tag
  query_circulator_or_iterator( const Self&) {
    return Iterator_tag();
  }
#endif // CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};

template < class I1, class I2, class I3, class I4, class I5,
           class  Creator >
class Join_input_iterator_5 {
  // the join of two iterators `i1' up to `i5'. Applies `Creator' with
  // five arguments `*i1' up to `*i5'. `value_type' is equal to
  // `Creator::result_type'.
public:
  typedef Join_input_iterator_5<I1,I2,I3,I4,I5,Creator> Self;

  typedef std::input_iterator_tag             iterator_category;
  typedef typename Creator::result_type       value_type;
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  typedef std::ptrdiff_t                      difference_type;
#else
  typedef std::iterator_traits<I1>            ITraits;
  typedef typename ITraits::difference_type   difference_type;
#endif
  typedef const value_type&                   reference;
  typedef const value_type*                   pointer;

protected:
  I1          j1;    // The 1st internal iterator.
  I2          j2;    // The 2nd internal iterator.
  I3          j3;    // The 3rd internal iterator.
  I4          j4;    // The 4th internal iterator.
  I5          j5;    // The 5th internal iterator.
  value_type  val;   // The current (internal) value.

public:
  // CREATION
  // --------

  Join_input_iterator_5() {}
  Join_input_iterator_5( I1 i1, I2 i2, I3 i3, I4 i4, I5 i5)
  : j1(i1), j2(i2), j3(i3), j4(i4), j5(i5),
    val(Creator()(*j1,*j2,*j3,*j4,*j5)) {}

    // OPERATIONS Forward Category
    // ---------------------------

    I1  current_iterator1() const { return j1;}
    I2  current_iterator2() const { return j2;}
    I3  current_iterator3() const { return j3;}
    I4  current_iterator4() const { return j4;}
    I5  current_iterator5() const { return j5;}

    bool operator==( const Self& i) const {
      return ( j1 == i.j1 &&
               j2 == i.j2 &&
               j3 == i.j3 &&
               j4 == i.j4 &&
               j5 == i.j5);
    }
    bool operator!=( const Self& i) const { return !(*this == i); }
    reference operator*() const { return val; }
    pointer   operator->() const { return &val; }
    Self& operator++() {
      ++j1;
      ++j2;
      ++j3;
      ++j4;
      ++j5;
      val = Creator()(*j1,*j2,*j3,*j4,*j5);
      return *this;
    }
    Self  operator++(int) {
      Self tmp = *this;
      ++*this;
      return tmp;
    }
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
#ifndef CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
    friend inline  value_type*
    value_type( const Self&) {
      return (value_type*)(0);
    }
    friend inline  std::input_iterator_tag
    iterator_category( const Self&){
      return std::input_iterator_tag();
    }
    friend inline  difference_type*
    distance_type( const Self&) {
      return (difference_type*)(0);
    }
    friend inline  Iterator_tag
    query_circulator_or_iterator( const Self&) {
      return Iterator_tag();
    }
#endif // CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};

CGAL_END_NAMESPACE
#endif // CGAL_JOIN_INPUT_ITERATOR_H //
// EOF //
