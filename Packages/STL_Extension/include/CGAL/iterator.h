// ============================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// file          : iterator.h
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
// package       : $CGAL_Package: STL_Extension $
// source        : stl_extension.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion <Sylvain.Pion@mpi-sb.mpg.de>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// Iterators and Iterator Adaptors
// ============================================================================


#ifndef CGAL_ITERATOR_H
#define CGAL_ITERATOR_H 1
#include <CGAL/circulator.h>
#include <vector>
#include <map>

CGAL_BEGIN_NAMESPACE

// +----------------------------------------------------------------+
// | Emptyset_iterator
// +----------------------------------------------------------------+
// |  sends everything to /dev/null
// +----------------------------------------------------------------+

struct Emptyset_iterator
#if defined(__GNUC__) && (__GNUC__ < 3)
  : public std::output_iterator
#else
  : public std::iterator< std::output_iterator_tag, void, void, void*, void >
#endif // defined(__GNUC__) && (__GNUC__ < 3)
{
  Emptyset_iterator() {}
  Emptyset_iterator(const Emptyset_iterator&) {}

  template< class T >
  Emptyset_iterator& operator=(const T&) { return *this; }

  Emptyset_iterator& operator++()        { return *this; }
  Emptyset_iterator& operator++(int)     { return *this; }

  Emptyset_iterator& operator*()         { return *this; }
};

// +---------------------------------------------------------------------+
// | Insert_iterator
// +---------------------------------------------------------------------+
// | Insert output iterator, which calls insert(value) on the container.
// | Similar to std::insert_iterator<> except it doesn't pass an iterator.
// +---------------------------------------------------------------------+

template < class Container >
class Insert_iterator
#if defined(__GNUC__) && (__GNUC__ < 3)
: public std::output_iterator
#else
: public std::iterator< std::output_iterator_tag, void, void, void*, void >
#endif // defined(__GNUC__) && (__GNUC__ < 3)
{
protected:
  Container *container;
public:
  typedef Container container_type;

  explicit Insert_iterator(Container &c)
  : container(&c) {}

  Insert_iterator&
  operator=(typename Container::const_reference value)
  {
    container->insert(value);
    return *this;
  }

  Insert_iterator&
  operator*() { return *this; }

  Insert_iterator&
  operator++() { return *this; }

  Insert_iterator
  operator++(int) { return *this; }
};

template < class Container >
inline Insert_iterator<Container>
inserter(Container &x)
{ return Insert_iterator<Container>(x); }

// +----------------------------------------------------------------+
// | Oneset_iterator
// +----------------------------------------------------------------+
// |  stores a reference to an object of type T
// |  which will be affected by operator*().
// +----------------------------------------------------------------+

template < class T >
class Oneset_iterator
#if defined(__GNUC__) && (__GNUC__ < 3)
  : public std::output_iterator
#else
  : public std::iterator< std::output_iterator_tag, void, void, void*, void >
#endif // defined(__GNUC__) && (__GNUC__ < 3)
{
  T& t;
public:
  Oneset_iterator(T& tt) : t(tt) {}

  Oneset_iterator& operator++()    { return *this; }
  Oneset_iterator& operator++(int) { return *this; }

  T& operator*() { return t; }
};

#if defined(CGAL_CFG_NO_ITERATOR_TRAITS) && \
!defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
template < class I, class Val>
#else
template < class I,
           class Val = CGAL_TYPENAME_MSVC_NULL
                       std::iterator_traits<I>::value_type >
#endif
class Counting_iterator {
protected:
  I            nt;    // The internal iterator.
  std::size_t  d_i;   // The internal counter.
public:
  typedef I  Iterator;
  typedef Counting_iterator<I,Val> Self;

  typedef std::input_iterator_tag  iterator_category;
  typedef Val                      value_type;
  typedef std::ptrdiff_t           difference_type;
  typedef const value_type&        reference;
  typedef const value_type*        pointer;

  // CREATION
  // --------

  Counting_iterator( std::size_t i = 0)             : d_i(i) {}
  Counting_iterator( Iterator j, std::size_t i = 0) : nt(j), d_i(i) {}

  // OPERATIONS Forward Category
  // ---------------------------

  Iterator    current_iterator() const { return nt;}
  std::size_t current_counter()  const { return d_i;}

  bool operator==( const Self& i) const { return ( d_i == i.d_i); }
  bool operator!=( const Self& i) const { return !(*this == i);   }
  reference  operator*()  const { return *nt; }
  pointer    operator->() const { return nt.operator->(); }
  Self& operator++() {
    ++nt;
    ++d_i;
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
  value_type( const Self&) { return (value_type*)(0); }
  friend inline  iterator_category
  iterator_category( const Self&){ return iterator_category(); }
  friend inline  difference_type*
  distance_type( const Self&) { return (difference_type*)(0); }
  friend inline  Iterator_tag
  query_circulator_or_iterator( const Self&) { return Iterator_tag(); }
#endif // CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};
#ifdef __SUNPRO_CC
// sunpro 5.3 complains about multiply defined types value_type,
// reference etc. below, if we use std::iterator_traits directly.
namespace CGALi {
  template < class I >
  struct IT_rename {
    typedef typename std::iterator_traits<I>::reference         REF;
    typedef typename std::iterator_traits<I>::pointer           PTR;
    typedef typename std::iterator_traits<I>::value_type        VAL;
    typedef typename std::iterator_traits<I>::difference_type   DIF;
    typedef typename std::iterator_traits<I>::iterator_category CAT;
  };
}
#endif // __SUNPRO_CC

#if defined(CGAL_CFG_NO_ITERATOR_TRAITS) && \
!defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
template < class I, int N, class Ref, class Ptr,
           class Val, class Dist, class Ctg>
#else
#ifndef __SUNPRO_CC
template < class I,
           int N,
           class Ref  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::reference,
           class Ptr  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::pointer,
           class Val  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::value_type,
           class Dist = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::difference_type,
           class Ctg  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::iterator_category >
#else
template < class I,
           int N,
           class Ref  = typename CGALi::IT_rename<I>::REF,
           class Ptr  = typename CGALi::IT_rename<I>::PTR,
           class Val  = typename CGALi::IT_rename<I>::VAL,
           class Dist = typename CGALi::IT_rename<I>::DIF,
           class Ctg  = typename CGALi::IT_rename<I>::CAT >
#endif // __SUNPRO_CC
#endif
class N_step_adaptor {
protected:
  I        nt;    // The internal iterator.
public:
  typedef I                                        Iterator;
  typedef N_step_adaptor<I,N,Ref,Ptr,Val,Dist,Ctg> Self;
  typedef Ctg                                      iterator_category;
  typedef Val                                      value_type;
  typedef Dist                                     difference_type;
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  typedef Ref                                      reference;
  typedef Ptr                                      pointer;
#else
  typedef typename std::iterator_traits<I>::reference reference;
  typedef typename std::iterator_traits<I>::pointer   pointer;
#endif
  // Special for circulators.
  typedef I_Circulator_size_traits<iterator_category,I> C_S_Traits;
  typedef typename  C_S_Traits::size_type               size_type;

  // CREATION
  // --------

  N_step_adaptor() {}
  N_step_adaptor( Iterator j) : nt(j) {}

  template <class II>
  N_step_adaptor( const N_step_adaptor<II,N>& j)
    : nt( j.current_iterator()) {}

  // OPERATIONS Forward Category
  // ---------------------------

  // Circulator stuff.
  typedef  I  Circulator;
  Circulator  current_circulator() const { return nt;}

  Iterator  current_iterator() const { return nt;}
  bool operator==( CGAL_NULL_TYPE p) const {
    CGAL_assertion( p == 0);
    return ( nt == 0);
  }
  bool  operator!=( CGAL_NULL_TYPE p) const { return !(*this == p); }
  bool  operator==( const Self& i) const { return ( nt == i.nt); }
  bool  operator!=( const Self& i) const { return !(*this == i); }
  Ref   operator*()  const { return *nt; }
  Ptr   operator->() const { return nt.operator->(); }
  Self& operator++() {
    std::advance( nt, N);
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
    std::advance( nt, -N);
    return *this;
  }
  Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
  }

  // OPERATIONS Random Access Category
  // ---------------------------------

  Self  min_circulator() const { return Self( nt.min_circulator()); }
  Self& operator+=( difference_type n) {
    nt += difference_type(N * n);
    return *this;
  }
  Self  operator+( difference_type n) const {
    Self tmp = *this;
    tmp.nt += difference_type(N * n);
    return tmp;
  }
#ifdef CGAL_CFG_NO_CONSTANTS_IN_FUNCTION_TEMPLATES
  friend inline
  Self
  operator+( difference_type n, Self i) {
    i = i + n;
    return i;
  }
#endif // CGAL_CFG_NO_CONSTANTS_IN_FUNCTION_TEMPLATES //
  Self& operator-=( difference_type n) {
    return operator+=( -n);
  }
  Self  operator-( difference_type n) const {
    Self tmp = *this;
    return tmp += -n;
  }
  difference_type  operator-( const Self& i) const { return (nt-i.nt)/N;}
  Ref  operator[]( difference_type n) const {
    Self tmp = *this;
    tmp += n;
    return tmp.operator*();
  }
  bool operator<( const Self& i) const { return ( nt < i.nt); }
  bool operator>( const Self& i) const { return i < *this; }
  bool operator<=( const Self& i) const { return !(i < *this); }
  bool operator>=( const Self& i) const { return !(*this < i); }
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
#ifndef CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
  friend inline  iterator_category
  iterator_category( const Self&) { return iterator_category(); }
  friend inline  value_type*
  value_type( const Self&) { return (value_type*)(0); }
  friend inline  difference_type*
  distance_type( const Self&) { return (difference_type*)(0); }
  typedef _Circulator_traits<iterator_category> C_Traits;
  typedef typename  C_Traits::category  category;
  friend inline  category
  query_circulator_or_iterator( const Self&) { return category(); }
#endif // CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};
#ifndef CGAL_CFG_NO_CONSTANTS_IN_FUNCTION_TEMPLATES
template < class I, int N, class Ref, class Ptr,
           class Val, class Dist, class Ctg>
inline
N_step_adaptor<I,N,Ref,Ptr,Val,Dist,Ctg>
operator+( Dist n, N_step_adaptor<I,N,Ref,Ptr,Val,Dist,Ctg> i)
{ return i += n; }
#endif // CGAL_CFG_NO_CONSTANTS_IN_FUNCTION_TEMPLATES //
template < class I, int N>
class N_step_adaptor_derived : public I {
public:
    typedef I                               Iterator;
    typedef I                               Circulator;
    typedef N_step_adaptor_derived<I,N>     Self;
    typedef typename I::iterator_category   iterator_category;
    typedef typename I::value_type          value_type;
    typedef typename I::difference_type     difference_type;
    typedef typename I::reference           reference;
    typedef typename I::pointer             pointer;
    // Special for circulators.
    typedef I_Circulator_size_traits<iterator_category,I> C_S_Traits;
    typedef typename  C_S_Traits::size_type               size_type;

// CREATION
// --------

    N_step_adaptor_derived() {}
    N_step_adaptor_derived( Iterator j) : I(j) {}

    template <class II>
    N_step_adaptor_derived( const N_step_adaptor_derived<II,N>& j)
        : I( j.current_iterator()) {}

// OPERATIONS Forward Category
// ---------------------------

    Circulator current_circulator() const { return *this;}
    Iterator   current_iterator()   const { return *this;}

    Self& operator++() {
        std::advance( (I&)*this, N);
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
        std::advance( (I&)*this, -N);
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }

// OPERATIONS Random Access Category
// ---------------------------------

    Self  min_circulator() const { return Self( I::min_circulator()); }
    Self& operator+=( difference_type n) {
        I::operator+=( difference_type(N * n));
        return *this;
    }
    Self  operator+( difference_type n) const {
        Self tmp = *this;
        tmp += n;
        return tmp;
    }
#ifdef CGAL_CFG_NO_CONSTANTS_IN_FUNCTION_TEMPLATES
    friend inline
    Self
    operator+( difference_type n, Self i) {
        i = i + n;
        return i;
    }
#endif // CGAL_CFG_NO_CONSTANTS_IN_FUNCTION_TEMPLATES //
    Self& operator-=( difference_type n) {
        return operator+=( -n);
    }
    Self  operator-( difference_type n) const {
        Self tmp = *this;
        return tmp += -n;
    }
    difference_type  operator-( const Self& i) const {
        return (I::operator-(i)) / N;
    }
    reference  operator[]( difference_type n) const {
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
    typedef _Circulator_traits<iterator_category> C_Traits;
    typedef typename  C_Traits::category  category;
    friend inline  category
    query_circulator_or_iterator( const Self&) { return category(); }
#endif // CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};
#ifndef CGAL_CFG_NO_CONSTANTS_IN_FUNCTION_TEMPLATES
template < class I, int N>
inline
N_step_adaptor_derived<I,N>
operator+( typename N_step_adaptor_derived<I,N>::difference_type n,
           N_step_adaptor_derived<I,N> i)
{ return i += n; }
#endif // CGAL_CFG_NO_CONSTANTS_IN_FUNCTION_TEMPLATES //
#if defined(CGAL_CFG_NO_ITERATOR_TRAITS) && \
!defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
template < class I, class P, class Ref, class Ptr,
           class Val, class Dist, class Ctg >
#else
#ifndef __SUNPRO_CC
template < class I,
           class P,
           class Ref  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::reference,
           class Ptr  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::pointer,
           class Val  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::value_type,
           class Dist = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::difference_type,
           class Ctg  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::iterator_category>
#else
template < class I,
           class P,
           class Ref  = typename CGALi::IT_rename<I>::REF,
           class Ptr  = typename CGALi::IT_rename<I>::PTR,
           class Val  = typename CGALi::IT_rename<I>::VAL,
           class Dist = typename CGALi::IT_rename<I>::DIF,
           class Ctg  = typename CGALi::IT_rename<I>::CAT >
#endif // __SUNPRO_CC
#endif
struct Filter_iterator {
  typedef I                                            Iterator;
  typedef P                                            Predicate;
  typedef Filter_iterator<I,P,Ref,Ptr,Val,Dist,Ctg>    Self;
  typedef Ctg                                          iterator_category;
  typedef Val                                          value_type;
  typedef Dist                                         difference_type;
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  typedef Ref                                          reference;
  typedef Ptr                                          pointer;
#else
  typedef typename std::iterator_traits<I>::reference  reference;
  typedef typename std::iterator_traits<I>::pointer    pointer;
#endif
  // Special for circulators.
  typedef I_Circulator_size_traits<iterator_category,I> C_S_Traits;
  typedef typename  C_S_Traits::size_type               size_type;

protected:
  Iterator b_, e_;   // The range.
  Iterator c_;       // current position.
  Predicate p_;      // Leave out x <==> p_(x).
public:

  Filter_iterator() {}

  Filter_iterator(Iterator b, Iterator e, const Predicate& p)
  : b_(b), e_(e), c_(b), p_(p)
  {
    while (c_ != e_ && p_(c_))
      ++c_;
  }

  Filter_iterator(Iterator b, Iterator e, const Predicate& p, Iterator c)
  : b_(b), e_(e), c_(c), p_(p)
  {
    while (c_ != e_ && p_(c_))
      ++c_;
  }

  bool operator==(const Self& it) const {
    CGAL_precondition(b_ == it.b_ && e_ == it.e_);
    return c_ == it.c_;
  }
  bool operator!=(const Self& it) const { return !(*this == it); }

  Self& operator++() {
    do { ++c_; } while (c_ != e_ && p_(c_));
    return *this;
  }

  Self& operator--() {
    do {
      --c_;
    } while (p_(c_) && c_ != b_);
    return *this;
  }

  Self operator++(int) {
    Self tmp(*this);
    ++(*this);
    return tmp;
  }

  Self operator--(int) {
    Self tmp(*this);
    --(*this);
    return tmp;
  }

  reference operator*() const { return *c_;  }
  pointer operator->() const  { return &*c_; }

};

template < class I, class P >
inline Filter_iterator< I, P >
filter_iterator(I b, I e, const P& p)
{ return Filter_iterator< I, P >(b, e, p); }

template < class I, class P >
inline Filter_iterator< I, P >
filter_iterator(I b, I e, const P& p, I c)
{ return Filter_iterator< I, P >(b, e, p, c); }


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
template < class IC>
class Inverse_index {

  // DEFINITION
  //
  // The class Inverse_index<IC,T> constructs an inverse index for a
  // given range [i,j) of two iterators or circulators of type `IC' with the
  // value type `T'. The first element I in the
  // range [i,j) has the index 0. Consecutive elements are numbered
  // incrementally. The inverse index provides a query for a given iterator
  // or circulator k to retrieve its index number. For random access
  // iterators or circulators, it is done in constant time by subtracting i.
  // For other iterator categories, an STL `map' is used, which results in a
  // log j-i query time. A comparison operator `operator<' is needed for
  // `T*'.
  //
  // CREATION

protected:
  typedef std::map< const void*, std::size_t, std::less<const void*> >
    Index;
  Index   idx;
  IC      start;
  typedef typename Index::iterator        Index_iterator;
  typedef typename Index::const_iterator  Index_const_iterator;
  typedef typename Index::value_type      Item;

protected:
  void ini_idx( IC i, const IC& j, std::input_iterator_tag);
  void ini_idx( const IC& i, const IC& j, std::forward_iterator_tag){
    ini_idx( i, j, std::input_iterator_tag());
  }
  void ini_idx(const IC& i,const IC& j, std::bidirectional_iterator_tag){
    ini_idx( i, j, std::input_iterator_tag());
  }
  void ini_idx( const IC& i, const IC& j, Forward_circulator_tag) {
    ini_idx( i, j, std::input_iterator_tag());
  }
  void ini_idx( const IC& i, const IC& j, Bidirectional_circulator_tag){
    ini_idx( i, j, std::input_iterator_tag());
  }
  void ini_idx( const IC&, const IC&, std::random_access_iterator_tag){}
  void ini_idx( const IC&, const IC&, Random_access_circulator_tag){}

public:
  void init_index( const IC& i, const IC& j) {
#if !defined(CGAL_CFG_NO_ITERATOR_TRAITS) || \
    defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
    typedef typename std::iterator_traits<IC>::iterator_category ICC;
    ini_idx( i, j, ICC());
#else
    ini_idx( i, j, std::iterator_category( i));
#endif
  }

protected:
  void push_back( const IC& k, std::input_iterator_tag) {
    std::size_t d = idx.size();
    idx[ &*k] = d;
  }
  void push_back( const IC& k, std::forward_iterator_tag){
    push_back( k, std::input_iterator_tag());
  }
  void push_back( const IC& k, std::bidirectional_iterator_tag){
    push_back( k, std::input_iterator_tag());
  }
  void push_back( const IC& k, Forward_circulator_tag){
    push_back( k, std::input_iterator_tag());
  }
  void push_back( const IC& k, Bidirectional_circulator_tag){
    push_back( k, std::input_iterator_tag());
  }
  void push_back( const IC&, std::random_access_iterator_tag){}
  void push_back( const IC&, Random_access_circulator_tag){}

public:
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

  std::size_t find( const IC& k, std::random_access_iterator_tag) const {
    return std::size_t(k - start);
  }
  std::size_t find( const IC& k, Random_access_circulator_tag) const {
    return std::size_t(k - start);
  }
  std::size_t find( const IC& k, std::input_iterator_tag) const {
    // returns inverse index of k.
    Index_const_iterator i = idx.find( &*k);
    CGAL_assertion( i != idx.end());
    return (*i).second;
  }
  std::size_t find( const IC& k, std::forward_iterator_tag) const {
    return find( k, std::input_iterator_tag());
  }
  std::size_t find( const IC& k, std::bidirectional_iterator_tag
#if defined( _MSC_VER ) && (_MSC_VER < 1300 )
  , int dummy=0
#endif //_MSC_VER
  ) const {
    return find( k, std::input_iterator_tag());
  }
  std::size_t find( const IC& k, Forward_circulator_tag) const {
    return find( k, std::input_iterator_tag());
  }
  std::size_t find( const IC& k, Bidirectional_circulator_tag) const {
    return find( k, std::input_iterator_tag());
  }

  typedef IC           iterator;
  typedef IC           Circulator;
  typedef std::size_t  size_type;

  Inverse_index() : start(IC()) {}
  // invalid index.

  Inverse_index( const IC& i) : start(i) {};
  // empty inverse index initialized to start at i.

  Inverse_index( const IC& i, const IC& j) : start(i) {
    // inverse index initialized with range [i,j).
    init_index( i, j);
  }

  // OPERATIONS

  std::size_t operator[]( const IC& k) const {
    // returns inverse index of k.
#if !defined(CGAL_CFG_NO_ITERATOR_TRAITS) || \
    defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
    typedef typename std::iterator_traits<IC>::iterator_category
      category;
    return find( k, category());
#else
    return find( k, std::iterator_category( k));
#endif
  }
};

#if (defined(__GNUC__) && (__GNUC__ >= 3))
template < class IC>
void
Inverse_index< IC>::ini_idx( IC i, const IC& j, std::input_iterator_tag) {
  std::size_t n = 0;
  if ( ! is_empty_range( i, j)) {
    do {
      idx.insert(Item( &*i, n));
      n++;
    } while ((++i) != (j));
  }
}
#else
template < class IC>
void
Inverse_index< IC>::ini_idx( IC i, const IC& j, std::input_iterator_tag) {
  std::size_t n = 0;
  Index_iterator hint = idx.begin();
  if ( ! is_empty_range( i, j)) {
    do {
      hint = idx.insert( hint, Item( &*i, n));
      n++;
    } while ((++i) != (j));
  }
}
#endif // (__GNUC__ >= 3)
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
template < class IC, class T >
class Random_access_value_adaptor : public Random_access_adaptor<IC> {
public:
  typedef typename Random_access_adaptor<IC>::size_type size_type;

  Random_access_value_adaptor() {}
  // invalid index.

  Random_access_value_adaptor( const IC& i)
  : Random_access_adaptor<IC>(i) {}
  // empty random access index initialized to start at i.

  Random_access_value_adaptor( const IC& i, const IC& j)
  : Random_access_adaptor<IC>(i,j) {}
  // random access index initialized with range [i,j).

  // OPERATIONS

  T& operator[]( size_type n) const {
    // returns inverse index of k.
    return *(Random_access_adaptor<IC>::operator[](n));
  }
};

CGAL_END_NAMESPACE
#endif // CGAL_ITERATOR_H //
// EOF //
