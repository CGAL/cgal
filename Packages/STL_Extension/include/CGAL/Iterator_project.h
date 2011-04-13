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
// release_date  : 1999, July 28
//
// file          : Iterator_project.h
// package       : STL_Extension (2.7)
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
// source        : stl_extension.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// An iterator adaptor performing a projection on the value type.
// ======================================================================

#ifndef CGAL_ITERATOR_PROJECT_H
#define CGAL_ITERATOR_PROJECT_H 1
#ifndef CGAL_CIRCULATOR_H
#include <CGAL/circulator.h>
#endif

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
template < class I, class Fct, class Ref, class Ptr, class Dist, class Ctg>
#else
template < class I,
           class Fct,
           class Ref  = typename std::iterator_traits<I>::reference,
           class Ptr  = typename std::iterator_traits<I>::pointer,
           class Dist = typename std::iterator_traits<I>::difference_type,
           class Ctg = typename std::iterator_traits<I>::iterator_category>
#endif
class Iterator_project {
protected:
    I        nt;    // The internal iterator.
public:
    typedef  I  Iterator;
    typedef  Iterator_project<I,Fct,Ref,Ptr,Dist,Ctg> Self;

    typedef  Ctg                          iterator_category;
    typedef  typename Fct::argument_type  argument_type;
    typedef  typename Fct::result_type    value_type;
    typedef  Ref                          reference;
    typedef  Ptr                          pointer;
    typedef  Dist                         difference_type;


// CREATION
// --------

    Iterator_project() {}
    Iterator_project( I j) : nt(j) {}

// OPERATIONS Forward Category
// ---------------------------

    Iterator  current_iterator() const { return nt;}
    Ptr       ptr() const {
        Fct fct;
        return &(fct(*nt));
    }
    bool  operator==( const Self& i) const { return ( nt == i.nt); }
    bool  operator!=( const Self& i) const { return !(*this == i); }
    Ref   operator*()  const { return *ptr(); }
    Ptr   operator->() const { return ptr(); }
    Self& operator++() {
        ++nt;
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
        --nt;
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
        nt += n;
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
        return nt - i.nt;
    }
    Ref  operator[]( difference_type n) const {
        Self tmp = *this;
        tmp += n;
        return tmp.operator*();
    }
    bool operator<( const Self& i) const {
        return ( nt < i.nt);
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
    value_type( const Self&) { return (value_type*)(0); }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};

template < class D, class Fct, class I, class Ref, class Ptr,
           class Dist, class Ctg>
inline
Iterator_project<I,Fct,Ref,Ptr,Dist,Ctg>
operator+( D n, Iterator_project<I,Fct,Ref,Ptr,Dist,Ctg> i) {
    return i += Dist(n);
}

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
template < class I, class Fct, class Ref, class Ptr, class Dist, class Ctg>
inline  Ctg
iterator_category( const Iterator_project<I,Fct,Ref,Ptr,Dist,Ctg>&) {
    return Ctg();
}
template < class I, class Fct, class Ref, class Ptr, class Dist, class Ctg>
inline  Dist*
distance_type( const Iterator_project<I,Fct,Ref,Ptr,Dist,Ctg>&) {
    return (Dist*)(0);
}
template < class I, class Fct, class Ref, class Ptr, class Dist, class Ctg>
inline  Iterator_tag
query_circulator_or_iterator(
    const Iterator_project<I,Fct,Ref,Ptr,Dist,Ctg>&) {
    return Iterator_tag();
}
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //


#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
template < class I, class II, class Fct, class Dist, class Ctg>
#else
template < class I,
           class II,
           class Fct,
           class Dist = typename std::iterator_traits<I>::difference_type,
           class Ctg = typename std::iterator_traits<I>::iterator_category>
#endif
class Iterator_const_project {
protected:
    I        nt;    // The internal iterator.
public:
    typedef  I  Iterator;
    typedef  Iterator_const_project<I,II,Fct,Dist,Ctg> Self;

    typedef  Ctg                          iterator_category;
    typedef  typename Fct::argument_type  argument_type;
    typedef  typename Fct::result_type    value_type;
    typedef  const value_type&            reference;
    typedef  const value_type*            pointer;
    typedef  Dist                         difference_type;

    typedef  reference                    Ref;
    typedef  pointer                      Ptr;

    typedef  Iterator_project<II,Fct,value_type&,value_type*,Dist,Ctg>
                                          mutable_iterator;

// CREATION
// --------

    Iterator_const_project() {}
    Iterator_const_project( Iterator j) : nt(j) {}
    Iterator_const_project( mutable_iterator j) : nt( &*j) {}

// OPERATIONS Forward Category
// ---------------------------

    Iterator  current_iterator() const { return nt;}
    Ptr       ptr() const {
        Fct fct;
        return &(fct(*nt));
    }

    bool  operator==( const Self& i) const { return ( nt == i.nt); }
    bool  operator!=( const Self& i) const { return !(*this == i); }
    Ref   operator*()  const { return *ptr(); }
    Ptr   operator->() const { return ptr(); }
    Self& operator++() {
        ++nt;
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
        --nt;
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
        nt += n;
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
        return nt - i.nt;
    }
    Ref  operator[]( difference_type n) const {
        Self tmp = *this;
        tmp += n;
        return tmp.operator*();
    }
    bool operator<( const Self& i) const {
        return ( nt < i.nt);
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
    value_type( const Self&) { return (value_type*)(0); }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};

template < class D, class Fct, class I, class II, class Dist, class Ctg>
inline
Iterator_const_project<I,II,Fct,Dist,Ctg>
operator+( D n, Iterator_const_project<I,II,Fct,Dist,Ctg> i) {
    return i += Dist(n);
}

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
template < class I, class II, class Fct, class Dist, class Ctg>
inline  Ctg
iterator_category( const Iterator_const_project<I,II,Fct,Dist,Ctg>&) {
    return Ctg();
}
template < class I, class II, class Fct, class Dist, class Ctg>
inline  Dist*
distance_type( const Iterator_const_project<I,II,Fct,Dist,Ctg>&) {
    return (Dist*)(0);
}
template < class I, class II, class Fct, class Dist, class Ctg>
inline  Iterator_tag
query_circulator_or_iterator(
    const Iterator_const_project<I,II,Fct,Dist,Ctg>&) {
    return Iterator_tag();
}
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //

CGAL_END_NAMESPACE
#endif // CGAL_ITERATOR_PROJECT_H //
// EOF //
