// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, September 02
//
// file          : circulator_compat.h
// package       : Circulator (3.4)
// chapter       : $CGAL_Chapter: Circulators $
// source        : circulator_compat.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// Old Circulator File for Compatibility with pre iterator_traits.
// ======================================================================

#ifndef CGAL_CIRCULATOR_COMPAT_H
#define CGAL_CIRCULATOR_COMPAT_H 1
// This file replaces circulator.h
#define CGAL_CIRCULATOR_H 1

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifndef CGAL_PROTECT_CSTDDEF
#include <cstddef>
#define CGAL_PROTECT_CSTDDEF
#endif
#ifndef CGAL_PROTECT_FUNCTIONAL
#include <functional>
#define CGAL_PROTECT_FUNCTIONAL
#endif
#ifndef CGAL_PROTECT_ITERATOR
#include <iterator>
#define CGAL_PROTECT_ITERATOR
#endif
#ifndef CGAL_CIRCULATOR_BASES_H
#include <CGAL/circulator_bases.h>
#endif

// CGAL defined Tags.
#ifndef CGAL_CFG_NO_ITERATOR_TRAITS
#define CGAL__CIRC_STL_ITERATOR_TRAITS
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //

#ifndef CGAL__CIRC_STL_ITERATOR_TRAITS

   // Try to figure out what STL is used.
#  ifdef __GNUG__
#      if __GNUC__ == 2 && __GNUC_MINOR__ == 7
#          define CGAL__CIRC_STL_GCC
#      else
#          define CGAL__CIRC_STL_SGI_3_0
#      endif

#  else // __GNUG__ //
        // Try to distinguish between HP and newest SGI STL (June 1997).
        // CGAL standard multiple inclusion protection does not harm here.
#ifndef CGAL_PROTECT_LIST
#include <list>
#define CGAL_PROTECT_LIST
#endif
#      ifdef LIST_H
#          ifdef  __sgi
               // Assume that SGI don't use HP STL any more. Then LIST_H
               // indicates an earlier SGI STL delivered with the C++
               // compiler release 7.1.
#              define CGAL__CIRC_STL_SGI_CC_1996
#          else  // __sgi //
               // On non SGI systems assume HP STL as default.
#              define CGAL__CIRC_STL_HP
#          endif // __sgi //
#      else  // LIST_H //
#          ifdef __SGI_STL_LIST_H
               // No else part. If its not the SGI STL, we don't know what
               // it is. Check if partial specialization and iterator
               // traits are supported. Check also for new 3.0 STL.
#              ifdef __STL_CLASS_PARTIAL_SPECIALIZATION
#                  define CGAL__CIRC_STL_ITERATOR_TRAITS
#              else  // __STL_CLASS_PARTIAL_SPECIALIZATION //
#                  if defined( __SGI_STL_INTERNAL_ITERATOR_H) || \
                     (defined( __SGI_STL_PORT) && __SGI_STL_PORT >= 0x2031)
#                      define CGAL__CIRC_STL_SGI_3_0
#                  else
#                      define CGAL__CIRC_STL_SGI_JUNE_1997
#                  endif // Release 3.0 //
#              endif // __STL_CLASS_PARTIAL_SPECIALIZATION //
#          endif // __SGI_STL_LIST_H //
#      endif // LIST_H //
#  endif // __GNUG__ //

#endif

#ifndef CGAL_NULL_TYPE
#if defined( __GNUG__ )
    // (__GNUC__ < 2 || (__GNUC__ == 2 && __GNUC_MINOR__ < 91))
#define CGAL_NULL_TYPE const void*
#define CGAL_CIRC_NULL 0
#else // __GNUG__ //
#define CGAL_NULL_TYPE int
#define CGAL_CIRC_NULL NULL
#endif // __GNUG__ //
#endif // CGAL_NULL_TYPE //

CGAL_BEGIN_NAMESPACE

template <class C>
struct _Circulator_traits {
    typedef  Iterator_tag  category;
};
CGAL_TEMPLATE_NULL
struct _Circulator_traits<Forward_circulator_tag> {
    typedef  Circulator_tag  category;
};
CGAL_TEMPLATE_NULL
struct _Circulator_traits<Bidirectional_circulator_tag> {
    typedef  Circulator_tag  category;
};
CGAL_TEMPLATE_NULL
struct _Circulator_traits<Random_access_circulator_tag> {
    typedef  Circulator_tag  category;
};

template <class Tag, class IC>
struct _Circulator_size_traits {
    typedef  std::size_t  size_type;
};
#ifndef CGAL_CFG_NO_ITERATOR_TRAITS
template <class C>
struct _Circulator_size_traits<Forward_circulator_tag,C> {
    typedef  typename  C::size_type  size_type;
};
template <class C>
struct _Circulator_size_traits<Bidirectional_circulator_tag,C> {
    typedef  typename  C::size_type  size_type;
};
template <class C>
struct _Circulator_size_traits<Random_access_circulator_tag,C> {
    typedef  typename  C::size_type  size_type;
};
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
#ifdef CGAL__CIRC_STL_ITERATOR_TRAITS
template <class C>
struct Circulator_traits {
    typedef std::iterator_traits<C>                 traits;
    typedef typename traits::iterator_category      iterator_category;
    typedef _Circulator_traits< iterator_category>  C_traits;
    typedef typename C_traits::category             category;
};

template <class C>
typename Circulator_traits<C>::category
query_circulator_or_iterator( const C&) {
    typedef typename Circulator_traits<C>::category category;
    return category();
}

template <class C>
struct _Iterator_from_circulator_traits {};

CGAL_TEMPLATE_NULL
struct _Iterator_from_circulator_traits< Forward_circulator_tag> {
    typedef  std::forward_iterator_tag        iterator_category;
};
CGAL_TEMPLATE_NULL
struct _Iterator_from_circulator_traits<
    Bidirectional_circulator_tag> {
    typedef  std::bidirectional_iterator_tag  iterator_category;
};
CGAL_TEMPLATE_NULL
struct _Iterator_from_circulator_traits<
    Random_access_circulator_tag> {
    typedef  std::random_access_iterator_tag  iterator_category;
};

template <class C>
struct _Circulator_from_iterator_traits {
};
CGAL_TEMPLATE_NULL
struct _Circulator_from_iterator_traits< std::forward_iterator_tag> {
    typedef  Forward_circulator_tag        iterator_category;
};
CGAL_TEMPLATE_NULL
struct _Circulator_from_iterator_traits< std::bidirectional_iterator_tag> {
    typedef  Bidirectional_circulator_tag  iterator_category;
};
CGAL_TEMPLATE_NULL
struct _Circulator_from_iterator_traits< std::random_access_iterator_tag> {
    typedef  Random_access_circulator_tag  iterator_category;
};    
#else  // CGAL__CIRC_STL_ITERATOR_TRAITS //
// No iterator traits.
// ===================
// Iterators
// ---------
template< class T, class D> inline
Iterator_tag
query_circulator_or_iterator( const std::input_iterator<T,D>&){
    return Iterator_tag();
}
inline
Iterator_tag
query_circulator_or_iterator(  const std::output_iterator&){
    return Iterator_tag();
}
template< class T, class D>   inline
Iterator_tag
query_circulator_or_iterator(  const std::forward_iterator<T,D>&){
    return Iterator_tag();
}
template< class T, class D>   inline
Iterator_tag
query_circulator_or_iterator(  const std::bidirectional_iterator<T,D>&){
    return Iterator_tag();
}
template< class T, class D>   inline
Iterator_tag
query_circulator_or_iterator(  const std::random_access_iterator<T,D>&){
    return Iterator_tag();
}
template< class T>   inline
Iterator_tag
query_circulator_or_iterator( T*){
    return Iterator_tag();
}
template< class T>   inline
Iterator_tag
query_circulator_or_iterator( const T*){
    return Iterator_tag();
}

// Circulators
// -----------
template< class T, class D, class S> inline
Circulator_tag
query_circulator_or_iterator(
        const Forward_circulator_base<T,D,S>&){
    return Circulator_tag();
}
template< class T, class D, class S> inline
Circulator_tag
query_circulator_or_iterator(
        const Bidirectional_circulator_base<T,D,S>&){
    return Circulator_tag();
}
template< class T, class D, class S> inline
Circulator_tag
query_circulator_or_iterator(
        const Random_access_circulator_base<T,D,S>&){
    return Circulator_tag();
}

// variant base classes
// --------------------
template< class T, class D, class S> inline
Circulator_tag
query_circulator_or_iterator(
        const Forward_circulator_ptrbase<T,D,S>&){
    return Circulator_tag();
}
template< class T, class D, class S> inline
Circulator_tag
query_circulator_or_iterator(
       const Bidirectional_circulator_ptrbase<T,D,S>&){
    return Circulator_tag();
}
template< class T, class D, class S> inline
Circulator_tag
query_circulator_or_iterator(
        const Random_access_circulator_ptrbase<T,D,S>&){
    return Circulator_tag();
}
// No iterator traits.
// ===================
// variant base classes
// --------------------
// The normal base classes could be handled automatically through inheritance.

template< class T, class D, class S> inline
Forward_circulator_tag
std::iterator_category(  const Forward_circulator_ptrbase<T,D,S>&){
    return Forward_circulator_tag();
}
template< class T, class D, class S> inline
Bidirectional_circulator_tag
std::iterator_category(  const Bidirectional_circulator_ptrbase<T,D,S>&){
    return Bidirectional_circulator_tag();
}
template< class T, class D, class S> inline
Random_access_circulator_tag
std::iterator_category(  const Random_access_circulator_ptrbase<T,D,S>&){
    return Random_access_circulator_tag();
}

template< class T, class D, class S> inline
Forward_circulator_tag
std::iterator_category( const Forward_circulator_base<T,D,S>&){
    return Forward_circulator_tag();
}
template< class T, class D, class S> inline
Bidirectional_circulator_tag
std::iterator_category( const Bidirectional_circulator_base<T,D,S>&){
    return Bidirectional_circulator_tag();
}
template< class T, class D, class S> inline
Random_access_circulator_tag
std::iterator_category( const Random_access_circulator_base<T,D,S>&){
    return Random_access_circulator_tag();
}
template <class T, class Dist, class Size> inline
T* std::value_type( const Forward_circulator_ptrbase<T,Dist,Size>&) {
    return (T*)(0);
}
template <class T, class Dist, class Size> inline
T* std::value_type( const Bidirectional_circulator_ptrbase<T,Dist,Size>&) {
    return (T*)(0);
}
template <class T, class Dist, class Size> inline
T* std::value_type( const Random_access_circulator_ptrbase<T,Dist,Size>&) {
    return (T*)(0);
}
template <class T, class Dist, class Size> inline
Dist* std::distance_type( const Forward_circulator_ptrbase<T,Dist,Size>&) {
    return (Dist*)(0);
}
template <class T, class Dist, class Size> inline
Dist* std::distance_type(
    const Bidirectional_circulator_ptrbase<T,Dist,Size>&) {
    return (Dist*)(0);
}
template <class T, class Dist, class Size> inline
Dist* std::distance_type(
    const Random_access_circulator_ptrbase<T,Dist,Size>&) {
    return (Dist*)(0);
}

template <class T, class Dist, class Size> inline
T* std::value_type( const Forward_circulator_base<T,Dist,Size>&) {
    return (T*)(0);
}
template <class T, class Dist, class Size> inline
T* std::value_type( const Bidirectional_circulator_base<T,Dist,Size>&) {
    return (T*)(0);
}
template <class T, class Dist, class Size> inline
T* std::value_type( const Random_access_circulator_base<T,Dist,Size>&) {
    return (T*)(0);
}
template <class T, class Dist, class Size> inline
Dist* std::distance_type( const Forward_circulator_base<T,Dist,Size>&) {
    return (Dist*)(0);
}
template <class T, class Dist, class Size> inline
Dist* std::distance_type( const Bidirectional_circulator_base<T,Dist,Size>&) {
    return (Dist*)(0);
}
template <class T, class Dist, class Size> inline
Dist* std::distance_type( const Random_access_circulator_base<T,Dist,Size>&) {
    return (Dist*)(0);
}
#endif // CGAL__CIRC_STL_ITERATOR_TRAITS //
#  ifndef CGAL__CIRC_STL_ITERATOR_TRAITS
#      ifdef CGAL__CIRC_STL_SGI_3_0
#ifndef __STL_NON_TYPE_TMPL_PARAM_BUG
template <class T, class Ref, class Ptr, std::size_t BufSiz>
struct __deque_iterator;
template <class T, class Ref, class Ptr, std::size_t BufSiz> inline
Iterator_tag
query_circulator_or_iterator( const __deque_iterator<T, Ref, Ptr,
                                                         BufSiz>&){
    return Iterator_tag();
}
#else /* __STL_NON_TYPE_TMPL_PARAM_BUG */
template <class T, class Ref, class Ptr>
struct __deque_iterator;
template <class T, class Ref, class Ptr> inline
Iterator_tag
query_circulator_or_iterator( const __deque_iterator<T, Ref, Ptr>&){
    return Iterator_tag();
}
#endif /* __STL_NON_TYPE_TMPL_PARAM_BUG */

template <class Value, class Key, class HashFcn,
          class ExtractKey, class EqualKey, class Alloc>
struct __hashtable_iterator;
template< class V, class K, class HF, class ExK, class EqK, class All>
inline
Iterator_tag
query_circulator_or_iterator(
        const __hashtable_iterator<V, K, HF, ExK, EqK, All>&){
    return Iterator_tag();
}
template <class Value, class Key, class HashFcn,
          class ExtractKey, class EqualKey, class Alloc>
struct __hashtable_const_iterator;
template< class V, class K, class HF, class ExK, class EqK, class All>
inline
Iterator_tag
query_circulator_or_iterator(
        const __hashtable_const_iterator<V, K, HF, ExK, EqK, All>&){
    return Iterator_tag();
}

template<class T, class Ref, class Ptr>
struct __list_iterator;
template< class T, class Ref, class Ptr> inline
Iterator_tag
query_circulator_or_iterator( const __list_iterator<T,Ref,Ptr>&){
    return Iterator_tag();
}

template<class CharT, class Alloc> class __rope_iterator;
template< class T, class Alloc> inline
Iterator_tag
query_circulator_or_iterator( const __rope_iterator<T,Alloc>&){
    return Iterator_tag();
}
template<class CharT, class Alloc> class __rope_const_iterator;
template< class T, class Alloc> inline
Iterator_tag
query_circulator_or_iterator( const __rope_const_iterator<T,Alloc>&){
    return Iterator_tag();
}

template <class T, class Ref, class Ptr>
struct __slist_iterator;
template< class T, class Ref, class Ptr> inline
Iterator_tag
query_circulator_or_iterator( const __slist_iterator<T,Ref,Ptr>&){
    return Iterator_tag();
}

template <class Value, class Ref, class Ptr>
struct __rb_tree_iterator;
template< class V, class Ref, class Ptr> inline
Iterator_tag
query_circulator_or_iterator( const __rb_tree_iterator<V,Ref,Ptr>&) {
    return Iterator_tag();
}
#      endif // CGAL__CIRC_STL_SGI_3_0 //
#      ifdef CGAL__CIRC_STL_SGI_JUNE_1997
#ifndef __STL_NON_TYPE_TMPL_PARAM_BUG
template <class T, class Ref, std::size_t BufSiz>
struct __deque_iterator;
template <class T, class Ref, std::size_t BufSiz> inline
Iterator_tag
query_circulator_or_iterator( const __deque_iterator<T, Ref, BufSiz>&){
    return Iterator_tag();
}
#else /* __STL_NON_TYPE_TMPL_PARAM_BUG */
template <class T, class Ref>
struct __deque_iterator;
template <class T, class Ref> inline
Iterator_tag
query_circulator_or_iterator( const __deque_iterator<T, Ref>&){
    return Iterator_tag();
}
#endif /* __STL_NON_TYPE_TMPL_PARAM_BUG */

template <class Value, class Key, class HashFcn,
          class ExtractKey, class EqualKey, class Alloc>
struct __hashtable_iterator;
template< class V, class K, class HF, class ExK, class EqK, class All>
inline
Iterator_tag
query_circulator_or_iterator(
        const __hashtable_iterator<V, K, HF, ExK, EqK, All>&){
    return Iterator_tag();
}
template <class Value, class Key, class HashFcn,
          class ExtractKey, class EqualKey, class Alloc>
struct __hashtable_const_iterator;
template< class V, class K, class HF, class ExK, class EqK, class All>
inline
Iterator_tag
query_circulator_or_iterator(
        const __hashtable_const_iterator<V, K, HF, ExK, EqK, All>&){
    return Iterator_tag();
}

template<class T, class Ref>
struct __list_iterator;
template< class T, class Ref> inline
Iterator_tag
query_circulator_or_iterator( const __list_iterator<T,Ref>&){
    return Iterator_tag();
}

template<class CharT, class Alloc> class __rope_iterator;
template< class T, class Alloc> inline
Iterator_tag
query_circulator_or_iterator( const __rope_iterator<T,Alloc>&){
    return Iterator_tag();
}
template<class CharT, class Alloc> class __rope_const_iterator;
template< class T, class Alloc> inline
Iterator_tag
query_circulator_or_iterator( const __rope_const_iterator<T,Alloc>&){
    return Iterator_tag();
}

template <class T, class Ref>
struct __slist_iterator;
template< class T, class Ref> inline
Iterator_tag
query_circulator_or_iterator( const __slist_iterator<T,Ref>&){
    return Iterator_tag();
}

template <class Value, class Ref>
struct __rb_tree_iterator;
template< class V, class Ref> inline
Iterator_tag
query_circulator_or_iterator( const __rb_tree_iterator<V,Ref>&) {
    return Iterator_tag();
}
#      endif // CGAL__CIRC_STL_SGI_JUNE_1997 //
#      ifdef CGAL__CIRC_STL_SGI_WWW_1996
template <class T>
struct __deque_iterator;
template <class T> inline
Iterator_tag
query_circulator_or_iterator( const __deque_iterator<T>&){
    return Iterator_tag();
}
template <class T>
struct __deque_const_iterator;
template <class T> inline
Iterator_tag
query_circulator_or_iterator( const __deque_const_iterator<T>&){
    return Iterator_tag();
}

template <class Value, class Key, class HashFcn,
          class ExtractKey, class EqualKey, class Alloc>
struct __hashtable_iterator;
template< class V, class K, class HF, class ExK, class EqK, class All>
inline
Iterator_tag
query_circulator_or_iterator(
        const __hashtable_iterator<V, K, HF, ExK, EqK, All>&){
    return Iterator_tag();
}
template <class Value, class Key, class HashFcn,
          class ExtractKey, class EqualKey, class Alloc>
struct __hashtable_const_iterator;
template< class V, class K, class HF, class ExK, class EqK, class All>
inline
Iterator_tag
query_circulator_or_iterator(
        const __hashtable_const_iterator<V, K, HF, ExK, EqK, All>&){
    return Iterator_tag();
}

template<class T>
struct __list_iterator;
template< class T> inline
Iterator_tag
query_circulator_or_iterator( const __list_iterator<T>&){
    return Iterator_tag();
}
template<class T>
struct __list_const_iterator;
template< class T> inline
Iterator_tag
query_circulator_or_iterator( const __list_const_iterator<T>&){
    return Iterator_tag();
}

template <class Value>
struct __rb_tree_iterator;
template< class V> inline
Iterator_tag
query_circulator_or_iterator( const __rb_tree_iterator<V>&) {
    return Iterator_tag();
}
template <class Value>
struct __rb_tree_const_iterator;
template< class V> inline
Iterator_tag
query_circulator_or_iterator( const __rb_tree_const_iterator<V>&) {
    return Iterator_tag();
}
#      endif // CGAL__CIRC_STL_SGI_WWW_1996 //
#      ifdef CGAL__CIRC_STL_SGI_CC_1996
#ifdef DEQUE_H
template< class T, class Alloc> inline
Iterator_tag
query_circulator_or_iterator(const typename deque<T,Alloc>::iterator&){
    return Iterator_tag();
}
template< class T, class Alloc> inline
Iterator_tag
query_circulator_or_iterator(
    const typename deque<T,Alloc>::const_iterator&){
    return Iterator_tag();
}
#endif // DEQUE_H //

#if defined(SGI_STL_HASHTABLE_H)
template< class V, class K, class HF, class ExK, class EqK, class All>
inline
Iterator_tag
query_circulator_or_iterator(
        const typename hashtable<V, K, HF, ExK, EqK, All>::iterator&){
    return Iterator_tag();
}
template< class V, class K, class HF, class ExK, class EqK, class All>
inline
Iterator_tag
query_circulator_or_iterator(
       const typename hashtable<V, K, HF, ExK, EqK, All>::const_iterator&){
    return Iterator_tag();
}
#endif // SGI_STL_HASHTABLE_H //

#ifdef LIST_H
template< class T, class Alloc> inline
Iterator_tag
query_circulator_or_iterator(const typename std::list<T,Alloc>::iterator&){
    return Iterator_tag();
}
template< class T, class Alloc> inline
Iterator_tag
query_circulator_or_iterator(
    const typename std::list<T,Alloc>::const_iterator&){
    return Iterator_tag();
}
#endif // LIST_H //

#if defined(TREE_H)
template< class K, class V, class KoV, class Cmp, class Alloc> inline
Iterator_tag
query_circulator_or_iterator(
        const typename rb_tree<K,V,KoV,Cmp,Alloc>::iterator&){
    return Iterator_tag();
}
template< class K, class V, class KoV, class Cmp, class Alloc> inline
Iterator_tag
query_circulator_or_iterator(
        const typename rb_tree<K,V,KoV,Cmp,Alloc>::const_iterator&){
    return Iterator_tag();
}
#endif // TREE_H //
#      endif // CGAL__CIRC_STL_SGI_CC_1996 //
#  endif // CGAL__CIRC_STL_ITERATOR_TRAITS //
/* A function that asserts a specific compile time tag */
/* forcing its two arguments to have equal type.       */
/* It is encapsulated with #ifdef since it will be defined also elsewhere. */
#ifndef CGAL_ASSERT_COMPILE_TIME_TAG
#define CGAL_ASSERT_COMPILE_TIME_TAG 1
template <class Base>
struct _Assert_tag_class {
    void match_compile_time_tag( const Base&) const {}
};
template< class Tag, class Derived>
inline void Assert_compile_time_tag( const Tag&, const Derived& b) {
    _Assert_tag_class<Tag> x;
    x.match_compile_time_tag(b);
}
#endif

template <class C> inline
void Assert_circulator( const C &c) {
    Assert_compile_time_tag(Circulator_tag(), query_circulator_or_iterator(c));
}
template <class I> inline
void Assert_iterator( const I &i) {
    Assert_compile_time_tag( Iterator_tag(), query_circulator_or_iterator(i));
}
template <class I> inline
void Assert_input_category( const I &i) {
    Assert_compile_time_tag( std::input_iterator_tag(),
                             std::iterator_category(i));
}
template <class I> inline
void Assert_output_category( const I &i) {
    Assert_compile_time_tag( std::output_iterator_tag(),
                             std::iterator_category(i));
}
template <class IC> inline
void Assert_forward_category( const IC &ic) {
    Assert_compile_time_tag( std::forward_iterator_tag(),
                             std::iterator_category(ic));
}
template <class IC> inline
void Assert_bidirectional_category( const IC &ic) {
    Assert_compile_time_tag( std::bidirectional_iterator_tag(),
                             std::iterator_category(ic));
}
template <class IC> inline
void Assert_random_access_category( const IC &ic) {
    Assert_compile_time_tag( std::random_access_iterator_tag(),
                             std::iterator_category(ic));
}

// The assert at-least-category functions use the following
// functions to resolve properly. Note the proper order of the
// arguments: 1st is the to be type, 2nd is the actual type.
inline void _Has_to_be_at_least( std::input_iterator_tag,
                                 std::input_iterator_tag){}
inline void _Has_to_be_at_least( std::input_iterator_tag,
                                 std::forward_iterator_tag){}
inline void _Has_to_be_at_least( std::input_iterator_tag,
                                 std::bidirectional_iterator_tag){}
inline void _Has_to_be_at_least( std::input_iterator_tag,
                                 std::random_access_iterator_tag){}

inline void _Has_to_be_at_least( std::output_iterator_tag,
                                 std::output_iterator_tag){}
inline void _Has_to_be_at_least( std::output_iterator_tag,
                                 std::forward_iterator_tag){}
inline void _Has_to_be_at_least( std::output_iterator_tag,
                                 std::bidirectional_iterator_tag){}
inline void _Has_to_be_at_least( std::output_iterator_tag,
                                 std::random_access_iterator_tag){}

inline void _Has_to_be_at_least( std::forward_iterator_tag,
                                 std::forward_iterator_tag){}
inline void _Has_to_be_at_least( std::forward_iterator_tag,
                                 std::bidirectional_iterator_tag){}
inline void _Has_to_be_at_least( std::forward_iterator_tag,
                                 std::random_access_iterator_tag){}

inline void _Has_to_be_at_least( std::bidirectional_iterator_tag,
                                 std::bidirectional_iterator_tag){}
inline void _Has_to_be_at_least( std::bidirectional_iterator_tag,
                                 std::random_access_iterator_tag){}

inline void _Has_to_be_at_least( std::random_access_iterator_tag,
                                 std::random_access_iterator_tag){}

// The is-at-least assertions.
template <class I> inline
void Assert_is_at_least_input_category( const I& i) {
    _Has_to_be_at_least(std::input_iterator_tag(), std::iterator_category(i));
}
template <class I> inline
void Assert_is_at_least_output_category( const I& i) {
    _Has_to_be_at_least(std::output_iterator_tag(),
                        std::iterator_category(i));
}
template <class IC> inline
void Assert_is_at_least_forward_category( const IC& ic) {
    _Has_to_be_at_least(std::forward_iterator_tag(),
                        std::iterator_category(ic));
}
template <class IC> inline
void Assert_is_at_least_bidirectional_category( const IC& ic) {
    _Has_to_be_at_least(std::bidirectional_iterator_tag(),
                        std::iterator_category(ic)) ;
}
template <class IC> inline
void Assert_is_at_least_random_access_category( const IC& ic) {
    _Has_to_be_at_least(std::random_access_iterator_tag(),
                        std::iterator_category(ic));
}

template< class C> inline
bool _is_empty_range( const C& c1, const C&, Circulator_tag){
    return c1 == CGAL_CIRC_NULL;
}

template< class I> inline
bool _is_empty_range( const I& i1, const I& i2, Iterator_tag){
    return i1 == i2;
}

template< class IC> inline
bool is_empty_range( const IC& ic1, const IC& ic2){
    // is `true' if the range [`ic1, ic2') is empty, `false' otherwise.
    // Precondition: `T' is either a circulator or an iterator type. The
    // range [`ic1, ic2') is valid.
    return _is_empty_range( ic1, ic2, query_circulator_or_iterator( ic1));
}

struct Circulator_or_iterator_tag {};  // any circulator or iterator.

inline
Circulator_or_iterator_tag
check_circulator_or_iterator( Circulator_tag ){
    return Circulator_or_iterator_tag();
}
inline
Circulator_or_iterator_tag
check_circulator_or_iterator( Iterator_tag ){
    return Circulator_or_iterator_tag();
}

template< class IC> inline
void Assert_circulator_or_iterator( const IC &ic){
    Assert_compile_time_tag(
        Circulator_or_iterator_tag(),
        check_circulator_or_iterator(
            query_circulator_or_iterator( ic)
        )
    );
}

#define CGAL_For_all( ic1, ic2) \
    for ( bool _circ_loop_flag = ! ::CGAL::is_empty_range( ic1, ic2); \
          _circ_loop_flag; \
          _circ_loop_flag = ((++ic1) != (ic2)) )

#define CGAL_For_all_backwards( ic1, ic2) \
    for ( bool _circ_loop_flag = ! ::CGAL::is_empty_range( ic1, ic2); \
          _circ_loop_flag; \
          _circ_loop_flag = ((ic1) != (--ic2)) )

// Note: these macros are superfluous now. See the more natural macros above.
// (Note that the macros below avoids problems with dangling else's.)
#define CGAL__For_all_old( ic1, ic2, body) \
    if ( ::CGAL::is_empty_range( ic1, ic2)); \
    else { \
        do \
            body \
        while ((++ic1) != (ic2)); \
    }
#define CGAL__For_all_backwards_old( ic1, ic2, body) \
    if ( ::CGAL::is_empty_range( ic1, ic2)); \
    else { \
        do \
            body \
        while ((ic1) != (--ic2)); \
    }

template <class T>
class Size_type_return_value_proxy {
public:
    typedef typename T::size_type  size_type;
    size_type x;
    Size_type_return_value_proxy( size_type _x) : x(_x) {}
    operator size_type() const { return x; }
};

template <class C> inline
Size_type_return_value_proxy<C>
_min_circulator_size( const C& c) {
    Assert_circulator(c);
    Assert_random_access_category(c);
    typedef typename C::size_type  size_type;
    size_type n = 0;
    if ( c != CGAL_CIRC_NULL) {
        n = (c-1) - c + 1;
        CGAL_assertion(n > 0);
    }
    return n;
}

template <class C>
Size_type_return_value_proxy<C>
_circulator_size( const C& c, Forward_circulator_tag) {
    // Simply count.
    if ( c == CGAL_CIRC_NULL)
        return 0;
    typedef typename C::size_type  size_type;
    size_type n = 0;
    C      d = c;
    do {
        ++n;
        ++d;
    } while( c != d);
    return n;
}
template <class C> inline
Size_type_return_value_proxy<C>
_circulator_size(const C& c, Bidirectional_circulator_tag) {
    return _circulator_size( c, Forward_circulator_tag());
}
template <class C> inline
Size_type_return_value_proxy<C>
_circulator_size(const C& c, Random_access_circulator_tag) {
    return _min_circulator_size( c.min_circulator());
}

template <class C> inline
Size_type_return_value_proxy<C>
circulator_size(const C& c) {
    return _circulator_size( c, std::iterator_category(c));
}
template <class T>
class Difference_type_return_value_proxy {
public:
    typedef typename T::difference_type  difference_type;
    difference_type x;
    Difference_type_return_value_proxy( difference_type _x) : x(_x) {}
    operator difference_type() const { return x; }
};

template <class C>
Difference_type_return_value_proxy<C>
_circulator_distance( C c, const C& d, Forward_circulator_tag) {
    // Simply count.
    if ( c == CGAL_CIRC_NULL)
        return 0;
    typedef typename C::difference_type  difference_type;
    difference_type n = 0;
    do {
        ++n;
    } while( ++c != d);
    return n;
}
template <class C> inline
Difference_type_return_value_proxy<C>
_circulator_distance(const C& c, const C& d, Bidirectional_circulator_tag){
    return _circulator_distance( c, d, Forward_circulator_tag());
}
template <class C> inline
Difference_type_return_value_proxy<C>
_circulator_distance(const C& c, const C& d, Random_access_circulator_tag){
    typedef typename C::difference_type  difference_type;
    typedef typename C::size_type        size_type;
    if ( d - c > 0)
        return (d - c);
    return difference_type(size_type(_min_circulator_size(
               c.min_circulator()))) - (c-d);
}

template <class C> inline
Difference_type_return_value_proxy<C>
circulator_distance(const C& c, const C& d) {
    return _circulator_distance( c, d, std::iterator_category(c));
}
template <class C> inline
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
std::ptrdiff_t
#else  // CGAL_CFG_NO_ITERATOR_TRAITS //
typename std::iterator_traits<C>::difference_type
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
_iterator_distance(const C& c1, const C& c2, Circulator_tag) {
    return circulator_distance( c1, c2);
}

template <class I> inline
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
std::ptrdiff_t
#else  // CGAL_CFG_NO_ITERATOR_TRAITS //
typename std::iterator_traits<I>::difference_type
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
_iterator_distance(const I& i1, const I& i2, Iterator_tag) {
    #ifdef CGAL_CFG_NO_ITERATOR_TRAITS
    std::ptrdiff_t dist = 0;
    std::distance( i1, i2, dist);
    return dist;
    #else
    return std::distance( i1, i2);
    #endif // CGAL_CFG_NO_ITERATOR_TRAITS //
}

template <class IC> inline
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
std::ptrdiff_t
#else  // CGAL_CFG_NO_ITERATOR_TRAITS //
typename std::iterator_traits<IC>::difference_type
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
iterator_distance(const IC& ic1, const IC& ic2) {
    return _iterator_distance( ic1, ic2,
                               query_circulator_or_iterator(ic1));
}
template <class C> inline
C _get_min_circulator( C c, Forward_circulator_tag) { return c; }
template <class C> inline
C _get_min_circulator( C c, Bidirectional_circulator_tag) { return c; }
template <class C> inline
C _get_min_circulator( C c, Random_access_circulator_tag) {
    return c.min_circulator();
}
template <class C> inline
C get_min_circulator( C c) {
    return _get_min_circulator( c, std::iterator_category(c));
}
template<class I, class U> inline
I non_negative_mod(I n, U m) {
    CGAL_assertion( m > 0);
    #if (-1 % 3) > 0
        n = n % m;
    #else
    if (n < 0)
        n = - (( - n - 1) % m) + m - 1;
    else
        n = n % m;
    #endif
    CGAL_assertion( n >= 0);
    return n;
}

template < class C, class Ref, class Ptr>
class Forward_iterator_from_circulator {
private:
    const C*  anchor;
    C         current;
    int       winding;
public:
//
// TYPES

    typedef  C                                           Circulator;
    typedef  Forward_iterator_from_circulator<C,Ref,Ptr> Self;
    typedef  std::forward_iterator_tag                   iterator_category;
    typedef  typename C::value_type                      value_type;
    typedef  typename C::difference_type                 difference_type;
    typedef  Ref                                         reference;
    typedef  Ptr                                         pointer;

//
// CREATION

    Forward_iterator_from_circulator() : anchor(0), winding(0) {}

    Forward_iterator_from_circulator( const C* circ, int n)
        : anchor( circ), current( *circ), winding(n) {}

//
// OPERATIONS

    bool operator==( const Self& i) const {
        CGAL_assertion( anchor == i.anchor);  // same anchor?
        return ( current == i.current) && ( winding == i.winding);
    }
    bool operator!=( const Self& i) const {
        return !(*this == i);
    }
    Ref  operator*() const {
        CGAL_assertion( anchor  != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        return Ref(*current);
    }
    Ptr  operator->() const {
        CGAL_assertion( anchor  != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        return Ptr(current.operator->());
    }
    Self& operator++() {
        CGAL_assertion( anchor  != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        if ( current == *anchor)
            ++winding;
        ++current;
        return *this;
    }
    Self  operator++(int) {
        Self tmp = *this;
        ++*this;
        return tmp;
    }
    Circulator  current_circulator() const { return current;}
};

template < class C, class Ref, class Ptr>
class Bidirectional_iterator_from_circulator {
private:
    const C*  anchor;
    C         current;
    int       winding;
public:
//
// TYPES

    typedef  C                                      Circulator;
    typedef  Bidirectional_iterator_from_circulator<C,Ref,Ptr> Self;
    typedef  std::bidirectional_iterator_tag        iterator_category;
    typedef  typename C::value_type                 value_type;
    typedef  typename C::difference_type            difference_type;
    typedef  Ref                                    reference;
    typedef  Ptr                                    pointer;

//
// CREATION

    Bidirectional_iterator_from_circulator() : anchor(0), winding(0) {}

    Bidirectional_iterator_from_circulator( const C* circ, int n)
        : anchor( circ), current( *circ), winding(n) {}

//
// OPERATIONS

    bool operator==( const Self& i) const {
        CGAL_assertion( anchor == i.anchor);  // same anchor?
        return ( current == i.current) && ( winding == i.winding);
    }
    bool operator!=( const Self& i) const {
        return !(*this == i);
    }
    Ref operator*() const {
        CGAL_assertion( anchor  != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        return Ref(*current);
    }
    Ptr  operator->() const {
        CGAL_assertion( anchor  != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        return Ptr(current.operator->());
    }
    Self& operator++() {
        CGAL_assertion( anchor  != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        if ( current == *anchor)
            ++winding;
        ++current;
        return *this;
    }
    Self  operator++(int) {
        Self tmp = *this;
        ++*this;
        return tmp;
    }

    Self& operator--() {
        CGAL_assertion( anchor != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        --current;
        if ( current == *anchor)
            --winding;
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }
    Circulator  current_circulator() const { return current;}
};

template < class C, class Ref, class Ptr>
class Random_access_iterator_from_circulator {
private:
    // The anchor is normalized to be a minimal circulator.
    const C*  anchor;
    C         current;
    int       winding;
public:
//
// TYPES

    typedef  C                                      Circulator;
    typedef  Random_access_iterator_from_circulator<C,Ref,Ptr> Self;
    typedef  std::random_access_iterator_tag        iterator_category;
    typedef  typename C::value_type                 value_type;
    typedef  typename C::difference_type            difference_type;
    typedef  Ref                                    reference;
    typedef  Ptr                                    pointer;

//
// CREATION

    Random_access_iterator_from_circulator() : anchor(0), winding(0) {}

    Random_access_iterator_from_circulator( const C* circ, int n)
        : anchor( circ), current( *circ), winding(n) {}

//
// OPERATIONS

    bool operator==( const Self& i) const {
        CGAL_assertion( anchor == i.anchor);  // same anchor?
        return ( current == i.current) && ( winding == i.winding);
    }
    bool operator!=( const Self& i) const {
        return !(*this == i);
    }
    Ref operator*() const {
        CGAL_assertion( anchor  != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        return Ref(*current);
    }
    Ptr  operator->() const {
        CGAL_assertion( anchor  != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        return Ptr(current.operator->());
    }
    Self& operator++() {
        CGAL_assertion( anchor  != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        ++current;
        if ( current == *anchor)
            ++winding;
        return *this;
    }
    Self  operator++(int) {
        Self tmp = *this;
        ++*this;
        return tmp;
    }

    Self& operator--() {
        CGAL_assertion( anchor != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        if ( current == *anchor)
            --winding;
        --current;
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }

    Self& operator+=( typename C::difference_type n) {
        CGAL_assertion( anchor != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        if ( n < 0 && current == *anchor)  // We are leaving the anchor.
            --winding;
        current += n;
        if ( n > 0 && current == *anchor)  // Back again at the anchor.
            ++winding;
        return *this;
    }
    Self  operator+( typename C::difference_type n) const {
        Self tmp = *this;
        return tmp += n;
    }
    Self& operator-=( typename C::difference_type n) {
        return operator+=( -n);
    }
    Self  operator-( typename C::difference_type n) const {
        Self tmp = *this;
        return tmp += -n;
    }
    typename C::difference_type  operator-( const Self& i) const;

    Ref  operator[](typename C::difference_type n) const {
        Self tmp = *this;
        tmp += n;
        return tmp.operator*();
    }
    bool operator<( const Self& i) const {
        CGAL_assertion( anchor  != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        CGAL_assertion( anchor  == i.anchor);
        return (     (winding < i.winding)
                 || (    (winding == i.winding)
                      && (current - *anchor) < (i.current - *anchor)
                    )
               );
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
    Circulator  current_circulator() const { return current;}
};

template < class Dist, class  C, class Ref, class Ptr>
Random_access_iterator_from_circulator<C,Ref,Ptr>
operator+( Dist n,
           const Random_access_iterator_from_circulator<C,Ref,Ptr>&
               circ) {
    Random_access_iterator_from_circulator<C,Ref,Ptr> tmp = circ;
    return tmp += n;
}

template < class  C, class Ref, class Ptr>
typename C::difference_type
Random_access_iterator_from_circulator<C,Ref,Ptr>::
operator-( const Random_access_iterator_from_circulator<C,Ref,Ptr>&
           i) const {
    CGAL_assertion( anchor  != CGAL_CIRC_NULL);
    CGAL_assertion( current != CGAL_CIRC_NULL);
    CGAL_assertion( anchor  == i.anchor);
    if ( winding != i.winding) {
        typename C::difference_type s = _min_circulator_size( *anchor);
        return   (current - *anchor) - (i.current - *anchor)
               + s * (winding - i.winding);
    }
    return (current - *anchor) - (i.current - *anchor);
}


template < class  C >
class Forward_container_from_circulator {
private:
    C anchor;
public:

// DEFINITION
//
// The adaptor Forward_container_from_circulator<C> is a class that
// converts any circulator type `C' to a kind of containerclass, i.e. a
// class that provides an iterator and a const iterator type and two
// member functions -- begin() and end() -- that return the appropriate
// forward iterators. In analogy to STL container classes these member
// functions return a const iterator in the case that the container itself
// is constant and a mutable iterator otherwise.
//
// For Forward_container_from_circulator<C> the circulator has to
// fulfill at least the requirements for a forward circulator. The similar
// adaptor `Bidirectional_container_from_circulator<C>' requires a
// bidirectional circulator to provide bidirectional iterators and the
// adaptor `Random_access_container_from_circulator<C>' requires a
// random access circulator to provide random access iterators. In this
// case the adaptor implements a total ordering relation that is currently
// not required for random access circulators. The total order is based on
// the difference value from all circulators to the circulator given at
// construction time. The difference value is a consistent ordering as
// stated in the requirements for random access circulators.
//
// PARAMETERS
//
// `C' is the appropriate circulator type.
//
// CREATION
//
// New creation variable is: `adaptor'

    Forward_container_from_circulator() {}
        // the resulting iterators will have a singular value.

    Forward_container_from_circulator(const C& c) : anchor(c) {}
        // the resulting iterators will have a singular value if the
        // circulator `c' is singular.

//
// TYPES

typedef  C                   Circulator;

typedef  typename C::value_type       value_type;
typedef  typename C::difference_type  difference_type;
typedef  typename C::size_type        size_type;
typedef  value_type&                  reference;
typedef  const value_type&            const_reference;
typedef  value_type*                  pointer;
typedef  const value_type*            const_pointer;

typedef  Forward_iterator_from_circulator<C,reference,pointer>
                             iterator;
typedef  Forward_iterator_from_circulator<C,const_reference,
         const_pointer>      const_iterator;

//
// OPERATIONS

    iterator begin() {
        // the start iterator.
        return iterator( &anchor, 0);
    }

    const_iterator begin() const {
        // the start const iterator.
        return const_iterator( &anchor, 0);
    }

    iterator end() {
        // the past-the-end iterator.
        return anchor == CGAL_CIRC_NULL ?  iterator( &anchor, 0)
                                       :  iterator( &anchor, 1);
    }
    const_iterator end() const {
        // the past-the-end const iterator.
        return anchor == CGAL_CIRC_NULL ?  const_iterator( &anchor, 0)
                                       :  const_iterator( &anchor, 1);
    }
    #if defined( __GNUG__ ) && defined( CGAL__CIRC_STL_GCC )
    friend inline difference_type* distance_type(const iterator&) {
        return (difference_type*)(0);
    }
    friend inline value_type* value_type(const iterator&) {
        return (value_type*)(0);
    }
    friend inline std::forward_iterator_tag
    iterator_category(iterator) {
        return std::forward_iterator_tag();
    }
    friend inline CGAL_Iterator_tag
    CGAL_query_circulator_or_iterator(iterator) {
        return CGAL_Iterator_tag();
    }
    
    friend inline difference_type* distance_type(const const_iterator&) {
        return (difference_type*)(0);
    }
    friend inline value_type* value_type(const const_iterator&) {
        return (value_type*)(0);
    }
    friend inline std::forward_iterator_tag
    iterator_category(const_iterator) {
        return std::forward_iterator_tag();
    }
    friend inline CGAL_Iterator_tag
    CGAL_query_circulator_or_iterator( const_iterator) {
        return CGAL_Iterator_tag();
    }
    #endif
};


template < class  C >
class Bidirectional_container_from_circulator {
private:
    C anchor;
public:
//
// CREATION

    Bidirectional_container_from_circulator() {}
        // the resulting iterators will have a singular value.

    Bidirectional_container_from_circulator(const C& c) : anchor(c) {}
        // the resulting iterators will have a singular value if the
        // circulator `c' is singular.

//
// TYPES

typedef C  Circulator;

typedef  typename C::value_type       value_type;
typedef  typename C::difference_type  difference_type;
typedef  typename C::size_type        size_type;
typedef  value_type&                  reference;
typedef  const value_type&            const_reference;
typedef  value_type*                  pointer;
typedef  const value_type*            const_pointer;

typedef  Bidirectional_iterator_from_circulator<C,reference,pointer>
                             iterator;
typedef  Bidirectional_iterator_from_circulator<C,const_reference,
         const_pointer>      const_iterator;

//
// OPERATIONS

    iterator begin() {
        // the start iterator.
        return iterator( &anchor, 0);
    }

    const_iterator begin() const {
        // the start const iterator.
        return const_iterator( &anchor, 0);
    }

    iterator end() {
        // the past-the-end iterator.
        return anchor == CGAL_CIRC_NULL ?  iterator( &anchor, 0)
                                       :  iterator( &anchor, 1);
    }
    const_iterator end() const {
        // the past-the-end const iterator.
        return anchor == CGAL_CIRC_NULL ?  const_iterator( &anchor, 0)
                                       :  const_iterator( &anchor, 1);
    }
    #if defined( __GNUG__ ) && defined( CGAL__CIRC_STL_GCC )
    friend inline difference_type* distance_type(const iterator&) {
        return (difference_type*)(0);
    }
    friend inline value_type* value_type(const iterator&) {
        return (value_type*)(0);
    }
    friend inline std::bidirectional_iterator_tag
    iterator_category(iterator) {
        return std::bidirectional_iterator_tag();
    }
    friend inline CGAL_Iterator_tag
    CGAL_query_circulator_or_iterator(iterator) {
        return CGAL_Iterator_tag();
    }
    
    friend inline difference_type* distance_type(const const_iterator&) {
        return (difference_type*)(0);
    }
    friend inline value_type* value_type(const const_iterator&) {
        return (value_type*)(0);
    }
    friend inline std::bidirectional_iterator_tag
    iterator_category(const_iterator) {
        return std::bidirectional_iterator_tag();
    }
    friend inline CGAL_Iterator_tag
    CGAL_query_circulator_or_iterator( const_iterator) {
        return CGAL_Iterator_tag();
    }
    #endif
};


template < class  C >
class Random_access_container_from_circulator {
private:
    C anchor;
public:
//
// CREATION

    Random_access_container_from_circulator() {}
        // the resulting iterators will have a singular value.

    Random_access_container_from_circulator(const C& c)
        // The anchor is normalized to be a minimal circulator.
        : anchor(c.min_circulator()) {}
        // the resulting iterators will have a singular value if the
        // circulator `c' is singular.

//
// TYPES

typedef C  Circulator;

typedef  typename C::value_type       value_type;
typedef  typename C::difference_type  difference_type;
typedef  typename C::size_type        size_type;
typedef  value_type&                  reference;
typedef  const value_type&            const_reference;
typedef  value_type*                  pointer;
typedef  const value_type*            const_pointer;

typedef  Random_access_iterator_from_circulator<C,reference,pointer>
                             iterator;
typedef  Random_access_iterator_from_circulator<C,const_reference,
         const_pointer>      const_iterator;

//
// OPERATIONS

    iterator begin() {
        // the start iterator.
        return iterator( &anchor, 0);
    }

    const_iterator begin() const {
        // the start const iterator.
        return const_iterator( &anchor, 0);
    }

    iterator end() {
        // the past-the-end iterator.
        return anchor == CGAL_CIRC_NULL ?  iterator( &anchor, 0)
                                       :  iterator( &anchor, 1);
    }
    const_iterator end() const {
        // the past-the-end const iterator.
        return anchor == CGAL_CIRC_NULL ?  const_iterator( &anchor, 0)
                                       :  const_iterator( &anchor, 1);
    }
    #if defined( __GNUG__ ) && defined( CGAL__CIRC_STL_GCC )
    friend inline difference_type* distance_type(const iterator&) {
        return (difference_type*)(0);
    }
    friend inline value_type* value_type(const iterator&) {
        return (value_type*)(0);
    }
    friend inline std::random_access_iterator_tag
    iterator_category(iterator) {
        return std::random_access_iterator_tag();
    }
    friend inline CGAL_Iterator_tag
    CGAL_query_circulator_or_iterator(iterator) {
        return CGAL_Iterator_tag();
    }
    
    friend inline difference_type* distance_type(const const_iterator&) {
        return (difference_type*)(0);
    }
    friend inline value_type* value_type(const const_iterator&) {
        return (value_type*)(0);
    }
    friend inline std::random_access_iterator_tag
    iterator_category(const_iterator) {
        return std::random_access_iterator_tag();
    }
    friend inline CGAL_Iterator_tag
    CGAL_query_circulator_or_iterator( const_iterator) {
        return CGAL_Iterator_tag();
    }
    #endif
};

#ifdef CGAL__CIRC_STL_ITERATOR_TRAITS

template < class  C, class Ref, class Ptr>
class Iterator_from_circulator {
private:
    // The anchor is normalized to be a minimal circulator.
    const C*  anchor;
    C         current;
    int       winding;

    typedef  std::iterator_traits<C>                      Itraits;
    typedef  typename  Itraits::iterator_category         Iter_cat;
    typedef  _Iterator_from_circulator_traits<Iter_cat>   ICtraits;

public:
//
// TYPES

    typedef C  Circulator;
    typedef Iterator_from_circulator<C,Ref,Ptr> Self;

    typedef typename ICtraits::iterator_category iterator_category;

    typedef typename C::value_type       value_type;
    typedef typename C::difference_type  difference_type;
    typedef typename C::size_type        size_type;
    typedef typename C::reference        reference;
    typedef typename C::pointer          pointer;

//
// CREATION

    Iterator_from_circulator() : anchor(0), winding(0) {}

    Iterator_from_circulator( const C* circ, int n)
        : anchor( circ), current( *circ), winding(n) {}

//
// OPERATIONS

    bool operator==( const Self& i) const {
        CGAL_assertion( anchor == i.anchor);  // same anchor?
        return ( current == i.current) && ( winding == i.winding);
    }
    bool operator!=( const Self& i) const {
        return !(*this == i);
    }
    Ref  operator*() const {
        CGAL_assertion( anchor  != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        return Ref(*current);
    }
    Ptr  operator->() const {
        CGAL_assertion( anchor  != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        return Ptr(current.operator->());
    }
    Self& operator++() {
        CGAL_assertion( anchor  != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        ++current;
        if ( current == *anchor)
            ++winding;
        return *this;
    }
    Self  operator++(int) {
        Self tmp = *this;
        ++*this;
        return tmp;
    }
    Self& operator--() {
        CGAL_assertion( anchor != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        if ( current == *anchor)
            --winding;
        --current;
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }
    Self& operator+=( difference_type n) {
        CGAL_assertion( anchor != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        if ( n < 0 && current == *anchor)  // We are leaving the anchor.
            --winding;
        current += n;
        if ( n > 0 && current == *anchor)  // Back again at the anchor.
            ++winding;
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
    difference_type  operator-( const Self& i) const;

    Ref  operator[](difference_type n) const {
        Self tmp = *this;
        tmp += n;
        return tmp.operator*();
    }
    bool operator<( const Self& i) const {
        CGAL_assertion( anchor  != CGAL_CIRC_NULL);
        CGAL_assertion( current != CGAL_CIRC_NULL);
        CGAL_assertion( anchor  == i.anchor);
        return (     (winding < i.winding)
                 || (    (winding == i.winding)
                      && (current - *anchor) < (i.current - *anchor)
                    )
               );
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
    Circulator  current_circulator() const { return current;}
};

template < class Dist, class  C, class Ref, class Ptr>
Iterator_from_circulator<C,Ref,Ptr>
operator+( Dist n,
           const Iterator_from_circulator<C,Ref,Ptr>& circ) {
    Iterator_from_circulator<C,Ref,Ptr> tmp = circ;
    return tmp += n;
}

template < class  C, class Ref, class Ptr>
typename C::difference_type
Iterator_from_circulator<C,Ref,Ptr>::
operator-( const Iterator_from_circulator<C,Ref,Ptr>& i) const {
    CGAL_assertion( anchor  != CGAL_CIRC_NULL);
    CGAL_assertion( current != CGAL_CIRC_NULL);
    CGAL_assertion( anchor  == i.anchor);
    if ( winding != i.winding) {
        difference_type s = _min_circulator_size( *anchor);
        return   (current - *anchor) - (i.current - *anchor)
               + s * (winding - i.winding);
    }
    return (current - *anchor) - (i.current - *anchor);
}


template < class  C >
class Container_from_circulator {
private:
    C anchor;
public:
//
// CREATION

    Container_from_circulator() {}
        // the resulting iterators will have a singular value.

    Container_from_circulator(const C& c)
        // The anchor is normalized to be a minimal circulator.
        : anchor(get_min_circulator(c)) {}
        // the resulting iterators will have a singular value if the
        // circulator `c' is singular.

//
// TYPES

typedef C  Circulator;

typedef typename C::value_type       value_type;
typedef value_type&                  reference;
typedef const value_type&            const_reference;
typedef value_type*                  pointer;
typedef const value_type*            const_pointer;
typedef typename C::size_type        size_type;
typedef typename C::difference_type  difference_type;

typedef Iterator_from_circulator< C, reference, pointer>
    iterator;
typedef Iterator_from_circulator< C, const_reference, const_pointer>
    const_iterator;
//
// OPERATIONS

    iterator begin() {
        // the start iterator.
        return iterator( &anchor, 0);
    }

    const_iterator begin() const {
        // the start const iterator.
        return const_iterator( &anchor, 0);
    }

    iterator end() {
        // the past-the-end iterator.
        return anchor == CGAL_CIRC_NULL ?  iterator( &anchor, 0)
                                       :  iterator( &anchor, 1);
    }
    const_iterator end() const {
        // the past-the-end const iterator.
        return anchor == CGAL_CIRC_NULL ?  const_iterator( &anchor, 0)
                                       :  const_iterator( &anchor, 1);
    }
};
#endif // CGAL__CIRC_STL_ITERATOR_TRAITS //
template < class  C>
class Forward_circulator_from_container
    : public  Forward_circulator_ptrbase<
            typename C::value_type, typename C::difference_type,
            typename C::size_type> {
private:
    typename C::iterator  i;
public:

// DEFINITION
//
// The adaptor Forward_circulator_from_container<C> is a class that
// provides a forward circulator for a container C as specified by the
// STL. The iterators belonging to the container C are supposed to be at
// least forward iterators. The adaptor for bidirectional circulators is
// `Bidirectional_circulator_from_container<C>' and
// `Random_access_circulator_from_container<C>' for random access
// circulators. Appropriate const circulators are also available.
//
// PARAMETERS
//
// `C' is the container type. The container is supposed to conform to the
// STL requirements for container (i.e. to have a `begin()' and an `end()'
// iterator as well as the local types `value_type', `size_type()', and
// `difference_type').
//
// CREATION

    Forward_circulator_from_container() {}

    Forward_circulator_from_container( C* c)
        :  Forward_circulator_ptrbase<
                typename C::value_type, typename C::difference_type,
                typename C::size_type>(c),
           i(c->begin()) {}

    Forward_circulator_from_container( C* c, typename C::iterator _i)
        :  Forward_circulator_ptrbase<
                typename C::value_type, typename C::difference_type,
                typename C::size_type>(c),
           i(_i) {}

//
// TYPES

    typedef C                              Container;
    typedef typename C::iterator           iterator;
    typedef typename C::const_iterator     const_iterator;
    typedef typename C::value_type         value_type;
    typedef typename C::reference          reference;
    typedef value_type*                    pointer;
    typedef typename C::size_type          size_type;
    typedef typename C::difference_type    difference_type;

//
// OPERATIONS

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        return    (_ptr == NULL)
               || (((C*)_ptr)->begin()==((C*)_ptr)->end());
    }
    bool operator!=( CGAL_NULL_TYPE p) const {
        return !(*this == p);
    }
    bool operator==( const Forward_circulator_from_container
                           <C>& c) const {
        return i == c.i;
    }
    bool operator!=( const Forward_circulator_from_container
                           <C>& c) const {
        return !(*this == c);
    }
    reference  operator*() const {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((C*)_ptr)->end());
        return *i;
    }
    pointer  operator->() const {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((C*)_ptr)->end());
        return &(*i);
    }
    Forward_circulator_from_container<C>&
    operator++() {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((C*)_ptr)->end());
        ++i;
        if ( i == ((C*)_ptr)->end())
            i = ((C*)_ptr)->begin();
        return *this;
    }
    Forward_circulator_from_container<C>
    operator++(int) {
        Forward_circulator_from_container<C> tmp= *this;
        ++*this;
        return tmp;
    }
    iterator  current_iterator() const { return i;}
    #if defined( __GNUG__ ) && defined( CGAL__CIRC_STL_GCC )
    friend inline typename C::difference_type* distance_type(
        const
    Forward_circulator_from_container<C>&) {
        return (typename C::difference_type*)(0);
    }
    friend inline typename C::value_type* value_type(
        const
    Forward_circulator_from_container<C>&) {
        return (typename C::value_type*)(0);
    }
    friend inline Forward_circulator_tag iterator_category(
    Forward_circulator_from_container<C>) {
        return Forward_circulator_tag();
    }
    friend inline CGAL_Circulator_tag CGAL_query_circulator_or_iterator(
    Forward_circulator_from_container<C>) {
        return CGAL_Circulator_tag();
    }
    #endif
};

template < class  C>
class Forward_const_circulator_from_container
    : public  Forward_circulator_ptrbase<
            typename C::value_type, typename C::difference_type,
            typename C::size_type> {
private:
    typename C::const_iterator  i;
public:
//
// CREATION

    Forward_const_circulator_from_container() {}

    Forward_const_circulator_from_container( const C* c)
        :  Forward_circulator_ptrbase<
            typename C::value_type, typename C::difference_type,
            typename C::size_type>((void*)c), i(c->begin()) {}

    Forward_const_circulator_from_container( const C* c,
                                             typename C::const_iterator _i)
        :  Forward_circulator_ptrbase<
                typename C::value_type, typename C::difference_type,
                typename C::size_type>((void*)c),
           i(_i) {}

//
// TYPES

    typedef C                              Container;
    typedef typename C::iterator           iterator;
    typedef typename C::const_iterator     const_iterator;
    typedef typename C::value_type         value_type;
    typedef typename C::const_reference    reference;
    typedef const value_type*              pointer;
    typedef typename C::size_type          size_type;
    typedef typename C::difference_type    difference_type;

//
// OPERATIONS

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        return    (_ptr == NULL)
               || (((const C*)_ptr)->begin()
                 ==((const C*)_ptr)->end());
    }
    bool operator!=( CGAL_NULL_TYPE p) const {
        return !(*this == p);
    }
    bool operator==(
        const Forward_const_circulator_from_container
                           <C>& c) const {
        return i == c.i;
    }
    bool operator!=( const Forward_const_circulator_from_container
                           <C>& c) const {
        return !(*this == c);
    }
    reference  operator*() const {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((const C*)_ptr)->end());
        return *i;
    }
    pointer  operator->() const {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((const C*)_ptr)->end());
        return &(*i);
    }
    Forward_const_circulator_from_container<C>&
    operator++() {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((const C*)_ptr)->end());
        ++i;
        if ( i == ((const C*)_ptr)->end())
            i = ((const C*)_ptr)->begin();
        return *this;
    }
    Forward_const_circulator_from_container<C>
    operator++(int) {
        Forward_const_circulator_from_container<C>
            tmp= *this;
        ++*this;
        return tmp;
    }
    const_iterator  current_iterator() const { return i;}
    #if defined( __GNUG__ ) && defined( CGAL__CIRC_STL_GCC )
    friend inline typename C::difference_type* distance_type(
        const
    Forward_const_circulator_from_container<C>&) {
        return (typename C::difference_type*)(0);
    }
    friend inline typename C::value_type* value_type(
        const
    Forward_const_circulator_from_container<C>&) {
        return (typename C::value_type*)(0);
    }
    friend inline Forward_circulator_tag iterator_category(
    Forward_const_circulator_from_container<C>) {
        return Forward_circulator_tag();
    }
    friend inline CGAL_Circulator_tag CGAL_query_circulator_or_iterator(
    Forward_const_circulator_from_container<C>) {
        return CGAL_Circulator_tag();
    }
    #endif
};


template < class  C>
class Bidirectional_circulator_from_container
    : public  Bidirectional_circulator_ptrbase<
            typename C::value_type, typename C::difference_type,
            typename C::size_type> {
private:
    typename C::iterator  i;
public:
//
// CREATION

    Bidirectional_circulator_from_container() {}

    Bidirectional_circulator_from_container( C* c)
        :  Bidirectional_circulator_ptrbase<
                typename C::value_type, typename C::difference_type,
                typename C::size_type>(c),
           i(c->begin()) {}

    Bidirectional_circulator_from_container( C* c,
                                                 typename C::iterator _i)
        :  Bidirectional_circulator_ptrbase<
                typename C::value_type, typename C::difference_type,
                typename C::size_type>(c),
           i(_i) {}

//
// TYPES

    typedef C                              Container;
    typedef typename C::iterator           iterator;
    typedef typename C::const_iterator     const_iterator;
    typedef typename C::value_type         value_type;
    typedef typename C::reference          reference;
    typedef value_type*                    pointer;
    typedef typename C::size_type          size_type;
    typedef typename C::difference_type    difference_type;

//
// OPERATIONS

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        return    (_ptr == NULL)
               || (((C*)_ptr)->begin()==((C*)_ptr)->end());
    }
    bool operator!=( CGAL_NULL_TYPE p) const {
        return !(*this == p);
    }
    bool operator==( const Bidirectional_circulator_from_container
                           <C>& c) const {
        return i == c.i;
    }
    bool operator!=( const Bidirectional_circulator_from_container
                           <C>& c) const {
        return !(*this == c);
    }
    reference  operator*() const {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((C*)_ptr)->end());
        return *i;
    }
    pointer  operator->() const {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((C*)_ptr)->end());
        return &(*i);
    }
    Bidirectional_circulator_from_container<C>&
    operator++() {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((C*)_ptr)->end());
        ++i;
        if ( i == ((C*)_ptr)->end())
            i = ((C*)_ptr)->begin();
        return *this;
    }
    Bidirectional_circulator_from_container<C>
    operator++(int) {
        Bidirectional_circulator_from_container<C> tmp= *this;
        ++*this;
        return tmp;
    }
    Bidirectional_circulator_from_container<C>&
    operator--() {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((C*)_ptr)->end());
        if ( i == ((C*)_ptr)->begin())
            i = ((C*)_ptr)->end();
        --i;
        return *this;
    }
    Bidirectional_circulator_from_container<C>
    operator--(int) {
        Bidirectional_circulator_from_container<C> tmp = *this;
        --*this;
        return tmp;
    }
    iterator  current_iterator() const { return i;}
    #if defined( __GNUG__ ) && defined( CGAL__CIRC_STL_GCC )
    friend inline typename C::difference_type* distance_type(
        const
    Bidirectional_circulator_from_container<C>&) {
        return (typename C::difference_type*)(0);
    }
    friend inline typename C::value_type* value_type(
        const
    Bidirectional_circulator_from_container<C>&) {
        return (typename C::value_type*)(0);
    }
    friend inline Bidirectional_circulator_tag iterator_category(
    Bidirectional_circulator_from_container<C>) {
        return Bidirectional_circulator_tag();
    }
    friend inline CGAL_Circulator_tag CGAL_query_circulator_or_iterator(
    Bidirectional_circulator_from_container<C>) {
        return CGAL_Circulator_tag();
    }
    #endif
};

template < class  C>
class Bidirectional_const_circulator_from_container
    : public  Bidirectional_circulator_ptrbase<
            typename C::value_type, typename C::difference_type,
            typename C::size_type> {
private:
    typename C::const_iterator  i;
public:
//
// CREATION

    Bidirectional_const_circulator_from_container() {}

    Bidirectional_const_circulator_from_container( const C* c)
        :  Bidirectional_circulator_ptrbase<
                typename C::value_type, typename C::difference_type,
                typename C::size_type>((void*)c), i(c->begin()) {}

    Bidirectional_const_circulator_from_container( const C* c,
                                             typename C::const_iterator _i)
        :  Bidirectional_circulator_ptrbase<
                typename C::value_type, typename C::difference_type,
                typename C::size_type>((void*)c),
           i(_i) {}

//
// TYPES

    typedef C                              Container;
    typedef typename C::iterator           iterator;
    typedef typename C::const_iterator     const_iterator;
    typedef typename C::value_type         value_type;
    typedef typename C::const_reference    reference;
    typedef const value_type*              pointer;
    typedef typename C::size_type          size_type;
    typedef typename C::difference_type    difference_type;

//
// OPERATIONS

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        return    (_ptr == NULL)
               || (((const C*)_ptr)->begin()
                 ==((const C*)_ptr)->end());
    }
    bool operator!=( CGAL_NULL_TYPE p) const {
        return !(*this == p);
    }
    bool operator==(
        const Bidirectional_const_circulator_from_container
                           <C>& c) const {
        return i == c.i;
    }
    bool operator!=(const Bidirectional_const_circulator_from_container
                           <C>& c) const {
        return !(*this == c);
    }
    reference  operator*() const {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((const C*)_ptr)->end());
        return *i;
    }
    pointer  operator->() const {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((const C*)_ptr)->end());
        return &(*i);
    }
    Bidirectional_const_circulator_from_container<C>&
    operator++() {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((const C*)_ptr)->end());
        ++i;
        if ( i == ((const C*)_ptr)->end())
            i = ((const C*)_ptr)->begin();
        return *this;
    }
    Bidirectional_const_circulator_from_container<C>
    operator++(int) {
        Bidirectional_const_circulator_from_container<C>
            tmp= *this;
        ++*this;
        return tmp;
    }
    Bidirectional_const_circulator_from_container<C>&
    operator--() {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((const C*)_ptr)->end());
        if ( i == ((const C*)_ptr)->begin())
            i = ((const C*)_ptr)->end();
        --i;
        return *this;
    }
    Bidirectional_const_circulator_from_container<C>
    operator--(int) {
        Bidirectional_const_circulator_from_container<C>
            tmp = *this;
        --*this;
        return tmp;
    }
    const_iterator  current_iterator() const { return i;}
    #if defined( __GNUG__ ) && defined( CGAL__CIRC_STL_GCC )
    friend inline typename C::difference_type* distance_type(
        const
    Bidirectional_const_circulator_from_container<C>&) {
        return (typename C::difference_type*)(0);
    }
    friend inline typename C::value_type* value_type(
        const
    Bidirectional_const_circulator_from_container<C>&) {
        return (typename C::value_type*)(0);
    }
    friend inline Bidirectional_circulator_tag iterator_category(
    Bidirectional_const_circulator_from_container<C>) {
        return Bidirectional_circulator_tag();
    }
    friend inline CGAL_Circulator_tag CGAL_query_circulator_or_iterator(
    Bidirectional_const_circulator_from_container<C>) {
        return CGAL_Circulator_tag();
    }
    #endif
};


template < class  C>
class Random_access_circulator_from_container
    : public  Random_access_circulator_ptrbase<
            typename C::value_type, typename C::difference_type,
            typename C::size_type> {
private:
    typename C::iterator  i;
public:
//
// CREATION

    Random_access_circulator_from_container() {}

    Random_access_circulator_from_container( C* c)
        :  Random_access_circulator_ptrbase<
                typename C::value_type, typename C::difference_type,
                typename C::size_type>(c),
           i(c->begin()) {}

    Random_access_circulator_from_container( C* c,
                                                 typename C::iterator _i)
        :  Random_access_circulator_ptrbase<
                typename C::value_type, typename C::difference_type,
                typename C::size_type>(c),
           i(_i) {}

//
// TYPES

    typedef C                              Container;
    typedef typename C::iterator           iterator;
    typedef typename C::const_iterator     const_iterator;
    typedef typename C::value_type         value_type;
    typedef typename C::reference          reference;
    typedef value_type*                    pointer;
    typedef typename C::size_type          size_type;
    typedef typename C::difference_type    difference_type;

//
// OPERATIONS

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        return    (_ptr == NULL)
               || (((C*)_ptr)->begin()==((C*)_ptr)->end());
    }
    bool operator!=( CGAL_NULL_TYPE p) const {
        return !(*this == p);
    }
    bool operator==( const Random_access_circulator_from_container
                           <C>& c) const {
        return i == c.i;
    }
    bool operator!=( const Random_access_circulator_from_container
                           <C>& c) const {
        return !(*this == c);
    }
    reference  operator*() const {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((C*)_ptr)->end());
        return *i;
    }
    pointer  operator->() const {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((C*)_ptr)->end());
        return &(*i);
    }
    Random_access_circulator_from_container<C>&
    operator++() {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((C*)_ptr)->end());
        ++i;
        if ( i == ((C*)_ptr)->end())
            i = ((C*)_ptr)->begin();
        return *this;
    }
    Random_access_circulator_from_container<C>
    operator++(int) {
        Random_access_circulator_from_container<C> tmp= *this;
        ++*this;
        return tmp;
    }

    Random_access_circulator_from_container<C>&
    operator--() {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((C*)_ptr)->end());
        if ( i == ((C*)_ptr)->begin())
            i = ((C*)_ptr)->end();
        --i;
        return *this;
    }
    Random_access_circulator_from_container<C>
    operator--(int) {
        Random_access_circulator_from_container<C> tmp = *this;
        --*this;
        return tmp;
    }

    Random_access_circulator_from_container<C>&
    operator+=( difference_type n);

    Random_access_circulator_from_container<C>
    operator+( difference_type n) const {
        Random_access_circulator_from_container<C> tmp = *this;
        return tmp += n;
    }
    Random_access_circulator_from_container<C>&
    operator-=( difference_type n) {
        return operator+=( -n);
    }
    Random_access_circulator_from_container<C>
    operator-( difference_type n) const {
        Random_access_circulator_from_container<C> tmp = *this;
        return tmp += -n;
    }
    difference_type
    operator-( const Random_access_circulator_from_container<
                   C>& c) const {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( c._ptr != NULL);
        return i - c.i;
    }
    reference  operator[](difference_type n) const {
        Random_access_circulator_from_container<C> tmp = *this;
        tmp += n;
        return tmp.operator*();
    }

    friend inline
    Random_access_circulator_from_container<C>
    operator+( difference_type n, const
               Random_access_circulator_from_container<C>& c) {
        Random_access_circulator_from_container<C> tmp = c;
        return tmp += n;
    }
    iterator  current_iterator() const { return i;}

    Random_access_circulator_from_container<C>
    min_circulator() const {
        return Random_access_circulator_from_container<C>((C*)_ptr);
    }
    #if defined( __GNUG__ ) && defined( CGAL__CIRC_STL_GCC )
    friend inline typename C::difference_type* distance_type(
        const
    Random_access_circulator_from_container<C>&) {
        return (typename C::difference_type*)(0);
    }
    friend inline typename C::value_type* value_type(
        const
    Random_access_circulator_from_container<C>&) {
        return (typename C::value_type*)(0);
    }
    friend inline Random_access_circulator_tag iterator_category(
    Random_access_circulator_from_container<C>) {
        return Random_access_circulator_tag();
    }
    friend inline CGAL_Circulator_tag CGAL_query_circulator_or_iterator(
    Random_access_circulator_from_container<C>) {
        return CGAL_Circulator_tag();
    }
    #endif
};

template < class C>
Random_access_circulator_from_container<C>&
Random_access_circulator_from_container<C>::
operator+=( typename C::difference_type n) {
    CGAL_assertion( _ptr != NULL);
    CGAL_assertion( i != ((C*)_ptr)->end());

    typename C::difference_type j    = i - ((C*)_ptr)->begin();
    typename C::difference_type size = ((C*)_ptr)->size();
    CGAL_assertion( j    >= 0);
    CGAL_assertion( size >= 0);
    j = non_negative_mod( j + n, size);
    CGAL_assertion( j >= 0);
    CGAL_assertion( j < size);
    i = ((C*)_ptr)->begin() + j;
    return *this;
}


template < class  C>
class Random_access_const_circulator_from_container
    : public  Random_access_circulator_ptrbase<
            typename C::value_type, typename C::difference_type,
            typename C::size_type> {
private:
    typename C::const_iterator  i;
public:
//
// CREATION

    Random_access_const_circulator_from_container() {}

    Random_access_const_circulator_from_container( const C* c)
        :  Random_access_circulator_ptrbase<
                typename C::value_type, typename C::difference_type,
                typename C::size_type>((void*)c), i(c->begin()) {}

    Random_access_const_circulator_from_container( const C* c,
                                             typename C::const_iterator _i)
        :  Random_access_circulator_ptrbase<
                typename C::value_type, typename C::difference_type,
                typename C::size_type>((void*)c), i(_i) {}

//
// TYPES

    typedef C                              Container;
    typedef typename C::iterator           iterator;
    typedef typename C::const_iterator     const_iterator;
    typedef typename C::value_type         value_type;
    typedef typename C::const_reference    reference;
    typedef const value_type*              pointer;
    typedef typename C::size_type          size_type;
    typedef typename C::difference_type    difference_type;

//
// OPERATIONS

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        return    (_ptr == NULL)
               || (((const C*)_ptr)->begin()
                 ==((const C*)_ptr)->end());
    }
    bool operator!=( CGAL_NULL_TYPE p) const {
        return !(*this == p);
    }
    bool operator==(
        const Random_access_const_circulator_from_container
                           <C>& c) const {
        return i == c.i;
    }
    bool operator!=(const Random_access_const_circulator_from_container
                           <C>& c) const {
        return !(*this == c);
    }
    reference  operator*() const {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((const C*)_ptr)->end());
        return *i;
    }
    pointer  operator->() const {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((const C*)_ptr)->end());
        return &(*i);
    }
    Random_access_const_circulator_from_container<C>&
    operator++() {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((const C*)_ptr)->end());
        ++i;
        if ( i == ((const C*)_ptr)->end())
            i = ((const C*)_ptr)->begin();
        return *this;
    }
    Random_access_const_circulator_from_container<C>
    operator++(int) {
        Random_access_const_circulator_from_container<C>
            tmp= *this;
        ++*this;
        return tmp;
    }

    Random_access_const_circulator_from_container<C>&
    operator--() {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( i != ((const C*)_ptr)->end());
        if ( i == ((const C*)_ptr)->begin())
            i = ((const C*)_ptr)->end();
        --i;
        return *this;
    }
    Random_access_const_circulator_from_container<C>
    operator--(int) {
        Random_access_const_circulator_from_container<C>
            tmp = *this;
        --*this;
        return tmp;
    }

    Random_access_const_circulator_from_container<C>&
    operator+=( difference_type n);

    Random_access_const_circulator_from_container<C>
    operator+( difference_type n) const {
        Random_access_const_circulator_from_container<C>
            tmp = *this;
        return tmp += n;
    }
    Random_access_const_circulator_from_container<C>&
    operator-=( difference_type n) {
        return operator+=( -n);
    }
    Random_access_const_circulator_from_container<C>
    operator-( difference_type n) const {
        Random_access_const_circulator_from_container<C>
            tmp = *this;
        return tmp += -n;
    }
    difference_type
    operator-( const Random_access_const_circulator_from_container<
                   C>& c) const {
        CGAL_assertion( _ptr != NULL);
        CGAL_assertion( c._ptr != NULL);
        return i - c.i;
    }
    reference  operator[](difference_type n) const {
        Random_access_const_circulator_from_container<C>
            tmp = *this;
        tmp += n;
        return tmp.operator*();
    }

    friend inline
    Random_access_const_circulator_from_container<C>
    operator+( difference_type n, const
               Random_access_const_circulator_from_container<
                   C>& c) {
        Random_access_const_circulator_from_container<C>
            tmp = c;
        return tmp += n;
    }
    const_iterator  current_iterator() const { return i;}

    Random_access_const_circulator_from_container<C>
    min_circulator() const {
        return Random_access_const_circulator_from_container<C>
                   ((const C*)_ptr);
    }
    #if defined( __GNUG__ ) && defined( CGAL__CIRC_STL_GCC )
    friend inline typename C::difference_type* distance_type(
        const
    Random_access_const_circulator_from_container<C>&) {
        return (typename C::difference_type*)(0);
    }
    friend inline typename C::value_type* value_type(
        const
    Random_access_const_circulator_from_container<C>&) {
        return (typename C::value_type*)(0);
    }
    friend inline Random_access_circulator_tag iterator_category(
    Random_access_const_circulator_from_container<C>) {
        return Random_access_circulator_tag();
    }
    friend inline CGAL_Circulator_tag CGAL_query_circulator_or_iterator(
    Random_access_const_circulator_from_container<C>) {
        return CGAL_Circulator_tag();
    }
    #endif
};

template < class C>
Random_access_const_circulator_from_container<C>&
Random_access_const_circulator_from_container<C>::
operator+=( typename C::difference_type n) {
    CGAL_assertion( _ptr != NULL);
    CGAL_assertion( i != ((const C*)_ptr)->end());

    typename C::difference_type j    = i - ((const C*)_ptr)
                                            ->begin();
    typename C::difference_type size = ((const C*)_ptr)->size();
    CGAL_assertion( j    >= 0);
    CGAL_assertion( size >= 0);
    j = non_negative_mod( j + n, size);
    CGAL_assertion( j >= 0);
    CGAL_assertion( j < size);
    i = ((const C*)_ptr)->begin() + j;
    return *this;
}
template < class I, class T, class Size, class Dist>
class Forward_circulator_from_iterator
    : public Forward_circulator_base<T,Dist,Size> {
private:
    I _begin;
    I _end;
    I current;
public:

// DEFINITION
//
// The adaptor Forward_circulator_from_iterator< I, T,  Size, Dist>
// is a class that converts two iterators, a begin and a past-the-end
// value, to a forward circulator. The iterators are supposed to be at
// least forward iterators. The adaptor for bidirectional circulators is
// `Bidirectional_circulator_from_iterator< I, T, Size, Dist>' and
// for random access circulators it is
// `Random_access_circulator_from_iterator< I, T, Size, Dist>'.
// Appropriate const circulators are also available.
//
// PARAMETERS
//
// `I' is the appropriate iterator type, `T' its value type, `Size' the
// unsigned integral value to hold the possible number of items in a
// sequence, and `Dist' is a signed integral value, the distance type
// between two iterators of the same sequence.
//
//
// TYPES

    typedef I         iterator;
    typedef T         value_type;
    typedef T&        reference;
    typedef T*        pointer;
    typedef Size      size_type;
    typedef Dist      difference_type;

    typedef Forward_circulator_from_iterator<I,T,Size,Dist> Self;

// CREATION
//
// New creation variable is: `adaptor'

    Forward_circulator_from_iterator()
        : _begin(I()), _end(I()), current(I()) {}
        // a circulator `adaptor' with a singular value.

    Forward_circulator_from_iterator( const I& begin, const I& end)
        : _begin( begin), _end( end), current( begin) {}
        // a circulator `adaptor' initialized to refer to the element
        // `*begin' in a range [`begin' ,`end' ). The circulator `adaptor'
        // contains a singular value if `begin==end'.

    Forward_circulator_from_iterator( const I& begin, const I& end,
                                          const I& cur)
        : _begin( begin), _end( end), current( cur) {}

    Forward_circulator_from_iterator( const Self& c, const I& cur)
        : _begin( c._begin), _end( c._end), current( cur) {}
//
// OPERATIONS

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        CGAL_assertion((_end == _begin) || (current != _end));
        return _end == _begin;
    }
    bool operator!=( CGAL_NULL_TYPE p) const {
        return !(*this == p);
    }
    bool operator==( const Self& c) const {
        // CGAL_assertion((_begin == c._begin) && (_end == c._end));
        return current == c.current;
    }
    bool operator!=( const Self& c) const {
        return !(*this == c);
    }
    reference  operator*() const {
        CGAL_assertion( current != _end);
        return *current;
    }
    pointer  operator->() const {
        CGAL_assertion( current != _end);
        return &(*current);
    }
    Self&  operator++() {
        CGAL_assertion( current != _end);
        ++current;
        if ( current == _end)
            current = _begin;
        return *this;
    }
    Self  operator++(int) {
        Self tmp= *this;
        ++*this;
        return tmp;
    }
    iterator  current_iterator() const { return current;}
    #if defined( __GNUG__ ) && defined( CGAL__CIRC_STL_GCC )
    friend inline  Dist* distance_type(
        const
    Forward_circulator_from_iterator<I, T, Size, Dist>&) {
        return ( Dist*)(0);
    }
    friend inline T* value_type(
        const
    Forward_circulator_from_iterator<I, T, Size, Dist>&) {
        return (T*)(0);
    }
    friend inline Forward_circulator_tag iterator_category(
    Forward_circulator_from_iterator<I, T, Size, Dist>) {
        return Forward_circulator_tag();
    }
    friend inline CGAL_Circulator_tag CGAL_query_circulator_or_iterator(
    Forward_circulator_from_iterator<I, T, Size, Dist>) {
        return CGAL_Circulator_tag();
    }
    #endif
};


template < class  I, class  T, class   Size, class  Dist >
class Forward_const_circulator_from_iterator
    : public Forward_circulator_base<T,Dist,Size> {
private:
    I _begin;
    I _end;
    I current;
public:
//
// TYPES

    typedef I         iterator;
    typedef T         value_type;
    typedef const T&  reference;
    typedef const T*  pointer;
    typedef Size      size_type;
    typedef Dist      difference_type;

    typedef Forward_const_circulator_from_iterator<I,T,Size,Dist> Self;

//
// CREATION

    Forward_const_circulator_from_iterator()
        : _begin(I()), _end(I()), current(I()) {}

    Forward_const_circulator_from_iterator( const I& begin,
                                                const I& end)
        : _begin( begin), _end( end), current( begin) {}

    Forward_const_circulator_from_iterator( const I& begin,
                                                const I& end,
                                                const I& cur)
        : _begin( begin), _end( end), current( cur) {}

    Forward_const_circulator_from_iterator( const Self& c,
                                                const I& cur)
        : _begin( c._begin), _end( c._end), current( cur) {}

//
// OPERATIONS

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        CGAL_assertion((_end == _begin) || (current != _end));
        return _end == _begin;
    }
    bool operator!=( CGAL_NULL_TYPE p) const {
        return !(*this == p);
    }
    bool operator==( const Self& c) const {
        // CGAL_assertion((_begin == c._begin) && (_end == c._end));
        return current == c.current;
    }
    bool operator!=( const Self& c) const {
        return !(*this == c);
    }
    reference  operator*() const {
        CGAL_assertion( current != _end);
        return *current;
    }
    pointer  operator->() const {
        CGAL_assertion( current != _end);
        return &(*current);
    }
    Self& operator++() {
        CGAL_assertion( current != _end);
        ++current;
        if ( current == _end)
            current = _begin;
        return *this;
    }
    Self  operator++(int) {
        Self tmp= *this;
        ++*this;
        return tmp;
    }
    iterator  current_iterator() const { return current;}
    #if defined( __GNUG__ ) && defined( CGAL__CIRC_STL_GCC )
    friend inline  Dist* distance_type(
        const
    Forward_const_circulator_from_iterator<I, T, Size, Dist>&) {
        return ( Dist*)(0);
    }
    friend inline T* value_type(
        const
    Forward_const_circulator_from_iterator<I, T, Size, Dist>&) {
        return (T*)(0);
    }
    friend inline Forward_circulator_tag iterator_category(
    Forward_const_circulator_from_iterator<I, T, Size, Dist>) {
        return Forward_circulator_tag();
    }
    friend inline CGAL_Circulator_tag CGAL_query_circulator_or_iterator(
    Forward_const_circulator_from_iterator<I, T, Size, Dist>) {
        return CGAL_Circulator_tag();
    }
    #endif
};


template < class  I, class  T, class   Size, class  Dist >
class Bidirectional_circulator_from_iterator
    : public Bidirectional_circulator_base<T,Dist,Size> {
private:
    I _begin;
    I _end;
    I current;
public:
//
// TYPES

    typedef I         iterator;
    typedef T         value_type;
    typedef T&        reference;
    typedef T*        pointer;
    typedef Size      size_type;
    typedef Dist      difference_type;

    typedef Bidirectional_circulator_from_iterator<I,T,Size,Dist> Self;
//
// CREATION

    Bidirectional_circulator_from_iterator()
        : _begin(I()), _end(I()), current(I()) {}

    Bidirectional_circulator_from_iterator( const I& begin,
                                                const I& end)
        : _begin( begin), _end( end), current( begin) {}

    Bidirectional_circulator_from_iterator( const I& begin,
                                                const I& end,
                                                const I& cur)
        : _begin( begin), _end( end), current( cur) {}

    Bidirectional_circulator_from_iterator( const Self& c,
                                                const I& cur)
        : _begin( c._begin), _end( c._end), current( cur) {}
//
// OPERATIONS

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        CGAL_assertion((_end == _begin) || (current != _end));
        return _end == _begin;
    }
    bool operator!=( CGAL_NULL_TYPE p) const {
        return !(*this == p);
    }
    bool operator==( const Self& c) const {
        // CGAL_assertion((_begin == c._begin) && (_end == c._end));
        return current == c.current;
    }
    bool operator!=( const Self& c) const {
        return !(*this == c);
    }
    reference  operator*() const {
        CGAL_assertion( current != _end);
        return *current;
    }
    pointer  operator->() const {
        CGAL_assertion( current != _end);
        return &(*current);
    }
    Self& operator++() {
        CGAL_assertion( current != _end);
        ++current;
        if ( current == _end)
            current = _begin;
        return *this;
    }
    Self  operator++(int) {
        Self tmp= *this;
        ++*this;
        return tmp;
    }
    Self& operator--() {
        CGAL_assertion( current != _end);
        if ( current == _begin)
            current = _end;
        --current;
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }
    iterator  current_iterator() const { return current;}
    #if defined( __GNUG__ ) && defined( CGAL__CIRC_STL_GCC )
    friend inline  Dist* distance_type(
        const
    Bidirectional_circulator_from_iterator<I, T, Size, Dist>&) {
        return ( Dist*)(0);
    }
    friend inline T* value_type(
        const
    Bidirectional_circulator_from_iterator<I, T, Size, Dist>&) {
        return (T*)(0);
    }
    friend inline Bidirectional_circulator_tag iterator_category(
    Bidirectional_circulator_from_iterator<I, T, Size, Dist>) {
        return Bidirectional_circulator_tag();
    }
    friend inline CGAL_Circulator_tag CGAL_query_circulator_or_iterator(
    Bidirectional_circulator_from_iterator<I, T, Size, Dist>) {
        return CGAL_Circulator_tag();
    }
    #endif
};

template < class  I, class  T, class   Size, class  Dist >
class Bidirectional_const_circulator_from_iterator
    : public Bidirectional_circulator_base<T,Dist,Size> {
private:
    I _begin;
    I _end;
    I current;
public:
//
// TYPES

    typedef I         iterator;
    typedef T         value_type;
    typedef const T&  reference;
    typedef const T*  pointer;
    typedef Size      size_type;
    typedef Dist      difference_type;

    typedef Bidirectional_const_circulator_from_iterator<I,T,Size,Dist>
            Self;
//
// CREATION

    Bidirectional_const_circulator_from_iterator()
        : _begin(I()), _end(I()), current(I()) {}

    Bidirectional_const_circulator_from_iterator( const I& begin,
                                                      const I& end)
        : _begin( begin), _end( end), current( begin) {}

    Bidirectional_const_circulator_from_iterator( const I& begin,
                                                      const I& end,
                                                      const I& cur)
        : _begin( begin), _end( end), current( cur) {}

    Bidirectional_const_circulator_from_iterator( const Self& c,
                                                      const I& cur)
        : _begin( c._begin), _end( c._end), current( cur) {}
//
// OPERATIONS

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        CGAL_assertion((_end == _begin) || (current != _end));
        return _end == _begin;
    }
    bool operator!=( CGAL_NULL_TYPE p) const {
        return !(*this == p);
    }
    bool operator==( const Self& c) const {
        // CGAL_assertion((_begin == c._begin) && (_end == c._end));
        return current == c.current;
    }
    bool operator!=( const Self& c) const {
        return !(*this == c);
    }
    reference  operator*() const {
        CGAL_assertion( current != _end);
        return *current;
    }
    pointer  operator->() const {
        CGAL_assertion( current != _end);
        return &(*current);
    }
    Self& operator++() {
        CGAL_assertion( current != _end);
        ++current;
        if ( current == _end)
            current = _begin;
        return *this;
    }
    Self  operator++(int) {
        Self tmp= *this;
        ++*this;
        return tmp;
    }
    Self& operator--() {
        CGAL_assertion( current != _end);
        if ( current == _begin)
            current = _end;
        --current;
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }
    iterator  current_iterator() const { return current;}
    #if defined( __GNUG__ ) && defined( CGAL__CIRC_STL_GCC )
    friend inline  Dist* distance_type(
        const
    Bidirectional_const_circulator_from_iterator<I, T, Size, Dist>&) {
        return ( Dist*)(0);
    }
    friend inline T* value_type(
        const
    Bidirectional_const_circulator_from_iterator<I, T, Size, Dist>&) {
        return (T*)(0);
    }
    friend inline Bidirectional_circulator_tag iterator_category(
    Bidirectional_const_circulator_from_iterator<I, T, Size, Dist>) {
        return Bidirectional_circulator_tag();
    }
    friend inline CGAL_Circulator_tag CGAL_query_circulator_or_iterator(
    Bidirectional_const_circulator_from_iterator<I, T, Size, Dist>) {
        return CGAL_Circulator_tag();
    }
    #endif
};

template < class  I, class  T, class   Size, class  Dist >
class Random_access_circulator_from_iterator
    : public Random_access_circulator_base<T,Dist,Size> {
private:
    I _begin;
    I _end;
    I current;
public:
//
// TYPES

    typedef I         iterator;
    typedef T         value_type;
    typedef T&        reference;
    typedef T*        pointer;
    typedef Size      size_type;
    typedef Dist      difference_type;

    typedef Random_access_circulator_from_iterator<I,T,Size,Dist> Self;
//
// CREATION

    Random_access_circulator_from_iterator()
        : _begin(I()), _end(I()), current(I()) {}

    Random_access_circulator_from_iterator( const I& begin,
                                                const I& end)
        : _begin( begin), _end( end), current( begin) {}

    Random_access_circulator_from_iterator( const I& begin,
                                                const I& end,
                                                const I& cur)
        : _begin( begin), _end( end), current( cur) {}

    Random_access_circulator_from_iterator( const Self& c,
                                                const I& cur)
        : _begin( c._begin), _end( c._end), current( cur) {}

//
// OPERATIONS

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        CGAL_assertion((_end == _begin) || (current != _end));
        return _end == _begin;
    }
    bool operator!=( CGAL_NULL_TYPE p) const {
        return !(*this == p);
    }
    bool operator==( const Self& c) const {
        // CGAL_assertion((_begin == c._begin) && (_end == c._end));
        return current == c.current;
    }
    bool operator!=( const Self& c) const {
        return !(*this == c);
    }
    reference  operator*() const {
        CGAL_assertion( current != _end);
        return *current;
    }
    pointer  operator->() const {
        CGAL_assertion( current != _end);
        return &(*current);
    }
    Self& operator++() {
        CGAL_assertion( current != _end);
        ++current;
        if ( current == _end)
            current = _begin;
        return *this;
    }
    Self  operator++(int) {
        Self tmp= *this;
        ++*this;
        return tmp;
    }
    Self& operator--() {
        CGAL_assertion( current != _end);
        if ( current == _begin)
            current = _end;
        --current;
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }
    Self& operator+=( Dist n);

    Self  operator+( Dist n) const {
        Self tmp = *this;
        return tmp += n;
    }
    Self& operator-=( Dist n) {
        return operator+=( -n);
    }
    Self  operator-( Dist n) const {
        Self tmp = *this;
        return tmp += -n;
    }
    Dist  operator-( const Self& i) const {
        CGAL_assertion((_begin == i._begin) && (_end == i._end));
        return current - i.current;
    }
    reference  operator[](Dist n) const {
        Self tmp = *this;
        tmp += n;
        return tmp.operator*();
    }
    iterator  current_iterator() const { return current;}

    Self  min_circulator() const {
        return Self( _begin, _end);
    }
    #if defined( __GNUG__ ) && defined( CGAL__CIRC_STL_GCC )
    friend inline  Dist* distance_type(
        const
    Random_access_circulator_from_iterator<I, T, Size, Dist>&) {
        return ( Dist*)(0);
    }
    friend inline T* value_type(
        const
    Random_access_circulator_from_iterator<I, T, Size, Dist>&) {
        return (T*)(0);
    }
    friend inline Random_access_circulator_tag iterator_category(
    Random_access_circulator_from_iterator<I, T, Size, Dist>) {
        return Random_access_circulator_tag();
    }
    friend inline CGAL_Circulator_tag CGAL_query_circulator_or_iterator(
    Random_access_circulator_from_iterator<I, T, Size, Dist>) {
        return CGAL_Circulator_tag();
    }
    #endif
};

template < class D, class I, class  T, class Size, class Dist> inline
Random_access_circulator_from_iterator< I, T, Size, Dist>
operator+( D n, const
    Random_access_circulator_from_iterator< I, T, Size, Dist>& circ) {
    Random_access_circulator_from_iterator< I, T, Size, Dist>
        tmp = circ;
    return tmp += Dist(n);
}

template < class I, class  T, class Size, class Dist>
Random_access_circulator_from_iterator< I, T, Size, Dist>&
Random_access_circulator_from_iterator< I, T, Size, Dist>::
operator+=( Dist n) {
    CGAL_assertion( current != _end);

    Dist i    = current - _begin;
    Dist size = _end    - _begin;
    CGAL_assertion( i    >= 0);
    CGAL_assertion( size >= 0);
    i = non_negative_mod( i + n, size);
    CGAL_assertion( i >= 0);
    CGAL_assertion( i < size);
    current = _begin + i;
    return *this;
}


template < class  I, class  T, class   Size, class  Dist >
class Random_access_const_circulator_from_iterator
    : public Random_access_circulator_base<T,Dist,Size> {
private:
    I _begin;
    I _end;
    I current;
public:
//
// TYPES

    typedef I         iterator;
    typedef T         value_type;
    typedef const T&  reference;
    typedef const T*  pointer;
    typedef Size      size_type;
    typedef Dist      difference_type;

    typedef Random_access_const_circulator_from_iterator<I,T,Size,Dist>
            Self;
//
// CREATION

    Random_access_const_circulator_from_iterator()
        : _begin(I()), _end(I()), current(I()) {}

    Random_access_const_circulator_from_iterator( const I& begin,
                                                      const I& end)
        : _begin( begin), _end( end), current( begin) {}

    Random_access_const_circulator_from_iterator( const I& begin,
                                                      const I& end,
                                                      const I& cur)
        : _begin( begin), _end( end), current( cur) {}

    Random_access_const_circulator_from_iterator( const Self& c,
                                                      const I& cur)
        : _begin( c._begin), _end( c._end), current( cur) {}

//
// OPERATIONS

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        CGAL_assertion((_end == _begin) || (current != _end));
        return _end == _begin;
    }
    bool operator!=( CGAL_NULL_TYPE p) const {
        return !(*this == p);
    }
    bool operator==( const Self& c) const {
        // CGAL_assertion((_begin == c._begin) && (_end == c._end));
        return current == c.current;
    }
    bool operator!=( const Self& c) const {
        return !(*this == c);
    }
    reference  operator*() const {
        CGAL_assertion( current != _end);
        return *current;
    }
    pointer  operator->() const {
        CGAL_assertion( current != _end);
        return &(*current);
    }
    Self& operator++() {
        CGAL_assertion( current != _end);
        ++current;
        if ( current == _end)
            current = _begin;
        return *this;
    }
    Self  operator++(int) {
        Self tmp= *this;
        ++*this;
        return tmp;
    }
    Self& operator--() {
        CGAL_assertion( current != _end);
        if ( current == _begin)
            current = _end;
        --current;
        return *this;
    }
    Self
    operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }
    Self& operator+=( Dist n);

    Self
    operator+( Dist n) const {
        Self tmp = *this;
        return tmp += n;
    }
    Self& operator-=( Dist n) {
        return operator+=( -n);
    }
    Self  operator-( Dist n) const {
        Self tmp = *this;
        return tmp += -n;
    }
    Dist  operator-( const Self& i) const {
        CGAL_assertion((_begin == i._begin) && (_end == i._end));
        return current - i.current;
    }
    reference  operator[](Dist n) const {
        Self tmp = *this;
        tmp += n;
        return tmp.operator*();
    }
    iterator  current_iterator() const { return current;}

    Self  min_circulator() const {
        return Self( _begin, _end);
    }
    #if defined( __GNUG__ ) && defined( CGAL__CIRC_STL_GCC )
    friend inline  Dist* distance_type(
        const
    Random_access_const_circulator_from_iterator<I, T, Size, Dist>&) {
        return ( Dist*)(0);
    }
    friend inline T* value_type(
        const
    Random_access_const_circulator_from_iterator<I, T, Size, Dist>&) {
        return (T*)(0);
    }
    friend inline Random_access_circulator_tag iterator_category(
    Random_access_const_circulator_from_iterator<I, T, Size, Dist>) {
        return Random_access_circulator_tag();
    }
    friend inline CGAL_Circulator_tag CGAL_query_circulator_or_iterator(
    Random_access_const_circulator_from_iterator<I, T, Size, Dist>) {
        return CGAL_Circulator_tag();
    }
    #endif
};

template < class D, class I, class  T, class Size, class Dist> inline
Random_access_const_circulator_from_iterator< I, T, Size, Dist>
operator+( D n, const
    Random_access_const_circulator_from_iterator< I, T, Size, Dist>&
    circ) {
    Random_access_const_circulator_from_iterator< I, T, Size, Dist>
        tmp = circ;
    return tmp += Dist(n);
}

template < class I, class  T, class Size, class Dist>
Random_access_const_circulator_from_iterator< I, T, Size, Dist>&
Random_access_const_circulator_from_iterator< I, T, Size, Dist>::
operator+=( Dist n) {
    CGAL_assertion( current != _end);

    Dist i    = current - _begin;
    Dist size = _end    - _begin;
    CGAL_assertion( i    >= 0);
    CGAL_assertion( size >= 0);
    i = non_negative_mod( i + n, size);
    CGAL_assertion( i >= 0);
    CGAL_assertion( i < size);
    current = _begin + i;
    return *this;
}

CGAL_END_NAMESPACE

#endif // CGAL_CIRCULATOR_COMPAT_H //
// EOF //
