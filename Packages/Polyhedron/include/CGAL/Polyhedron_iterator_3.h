// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : Polyhedron_iterator_3.h
// chapter       : $CGAL_Chapter: 3D-Polyhedral Surfaces $
// package       : $CGAL_Package: Polyhedron 2.9 (13 Sep 2000) $
// source        : polyhedron.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Iterator and Circulator for Polyhedral Surfaces.
// ============================================================================

#ifndef CGAL_POLYHEDRON_ITERATOR_3_H
#define CGAL_POLYHEDRON_ITERATOR_3_H 1
#ifndef CGAL_CIRCULATOR_H
#include <CGAL/circulator.h>
#endif

// Define shorter names to please linker (g++/egcs)
#define _Polyhedron_iterator                _PhI
#define _Polyhedron_const_iterator          _PhCI
#define _Polyhedron_edge_iterator           _PhEI
#define _Polyhedron_edge_const_iterator     _PhECI
#define _Polyhedron_facet_circ              _PhFC
#define _Polyhedron_facet_const_circ        _PhFCC
#define _Polyhedron_vertex_circ             _PhVC
#define _Polyhedron_vertex_const_circ       _PhVCC

CGAL_BEGIN_NAMESPACE

// The following two iterators are similar to Iterator_project
// and Iterator_const_project, but they implement the arrow
// operator -> in any case. This is here possible since the elements
// a polyhedron consists of are always classes with members.

template < class I, class Val, class Dist, class Ctg>
class _Polyhedron_iterator {
protected:
    I        nt;    // The internal iterator.
public:
    typedef  I  Iterator;
    typedef  _Polyhedron_iterator<I,Val,Dist,Ctg> Self;

    typedef  Ctg                          iterator_category;
    typedef  Val                          value_type;
    typedef  value_type&                  reference;
    typedef  value_type*                  pointer;
    typedef  Dist                         difference_type;

// CREATION
// --------

    _Polyhedron_iterator() {}
    _Polyhedron_iterator( I j) : nt(j) {}
    // _Polyhedron_iterator( Ptr p) : nt(p) {}

// OPERATIONS Forward Category
// ---------------------------

    Iterator  current_iterator() const { return nt;}
    pointer   ptr() const { return (pointer)(&(*nt)); }

    bool  operator==( const Self& i) const { return ( nt == i.nt); }
    bool  operator!=( const Self& i) const { return !(*this == i); }
    reference   operator*() const  { return *ptr(); }
    pointer     operator->() const { return  ptr(); }
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
    Self& operator-=( difference_type n) { return operator+=( -n); }
    Self  operator-( difference_type n) const {
        Self tmp = *this;
        return tmp += -n;
    }
    difference_type  operator-( const Self& i) const { return nt - i.nt; }
    reference  operator[]( difference_type n) const {
        Self tmp = *this;
        tmp += n;
        return tmp.operator*();
    }
    bool operator< ( const Self& i) const { return ( nt < i.nt); }
    bool operator> ( const Self& i) const { return i < *this;    }
    bool operator<=( const Self& i) const { return !(i < *this); }
    bool operator>=( const Self& i) const { return !(*this < i); }
};

template < class D, class I, class Val, class Dist, class Ctg>
inline
_Polyhedron_iterator<I,Val,Dist,Ctg>
operator+( D n, _Polyhedron_iterator<I,Val,Dist,Ctg> i) {
    return i += Dist(n);
}

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
template < class I, class Val, class Dist, class Ctg>
inline  Ctg
iterator_category( const _Polyhedron_iterator<I,Val,Dist,Ctg>&){
    return Ctg();
}
template < class I, class Val, class Dist, class Ctg>
inline  Val*
value_type( const _Polyhedron_iterator<I,Val,Dist,Ctg>&) {
    return (Val*)(0);
}
template < class I, class Val, class Dist, class Ctg>
inline  Dist*
distance_type( const _Polyhedron_iterator<I,Val,Dist,Ctg>&) {
    return (Dist*)(0);
}
template < class I, class Val, class Dist, class Ctg>
inline  Iterator_tag
query_circulator_or_iterator(
    const _Polyhedron_iterator<I,Val,Dist,Ctg>&)
{
    return Iterator_tag();
}
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //


template < class I, class II, class Val, class Dist, class Ctg>
class _Polyhedron_const_iterator {
protected:
    I        nt;    // The internal iterator.
public:
    typedef  I  Iterator;
    typedef  _Polyhedron_const_iterator<I,II,Val,Dist,Ctg> Self;

    typedef  Ctg                          iterator_category;
    typedef  Val                          value_type;
    typedef  const value_type&            reference;
    typedef  const value_type*            pointer;
    typedef  Dist                         difference_type;

    typedef  _Polyhedron_iterator<II,Val,Dist,Ctg>  mutable_iterator;

// CREATION
// --------

    _Polyhedron_const_iterator() {}
    _Polyhedron_const_iterator( Iterator j) : nt(j) {}
    // _Polyhedron_const_iterator( Ptr p) : nt(p) {}
    _Polyhedron_const_iterator( mutable_iterator j) : nt( I(&*j)) {}

// OPERATIONS Forward Category
// ---------------------------

    Iterator  current_iterator() const { return nt;}
    pointer   ptr() const { return (pointer)(&(*nt)); }
    bool operator==( const Self& i) const { return ( nt == i.nt); }
    bool operator!=( const Self& i) const { return !(*this == i); }
    reference  operator*()  const { return *ptr(); }
    pointer    operator->() const { return ptr();  }
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
    Self& operator-=( difference_type n) { return operator+=( -n); }
    Self  operator-( difference_type n) const {
        Self tmp = *this;
        return tmp += -n;
    }
    difference_type  operator-( const Self& i) const { return nt - i.nt; }
    reference  operator[]( difference_type n) const {
        Self tmp = *this;
        tmp += n;
        return tmp.operator*();
    }
    bool operator< ( const Self& i) const { return ( nt < i.nt); }
    bool operator> ( const Self& i) const { return i < *this;    }
    bool operator<=( const Self& i) const { return !(i < *this); }
    bool operator>=( const Self& i) const { return !(*this < i); }
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

template < class D, class I, class II, class Val, class Dist, class Ctg>
inline
_Polyhedron_const_iterator<I,II,Val,Dist,Ctg>
operator+( D n, _Polyhedron_const_iterator<I,II,Val,Dist,Ctg> i) {
    return i += Dist(n);
}

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
template < class I, class II, class Val, class Dist, class Ctg>
inline  Ctg
iterator_category(
    const _Polyhedron_const_iterator<I,II,Val,Dist,Ctg>&){
    return Ctg();
}
template < class I, class II, class Val, class Dist, class Ctg>
inline  Val*
value_type( const _Polyhedron_const_iterator<I,II,Val,Dist,Ctg>&) {
    return (Val*)(0);
}
template < class I, class II, class Val, class Dist, class Ctg>
inline  Dist*
distance_type( const _Polyhedron_const_iterator<I,II,Val,Dist,Ctg>&) {
    return (Dist*)(0);
}
template < class I, class II, class Val, class Dist, class Ctg>
inline  Iterator_tag
query_circulator_or_iterator(
    const _Polyhedron_const_iterator<I,II,Val,Dist,Ctg>&)
{
    return Iterator_tag();
}
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
template < class It, class Val, class Dist, class Ctg>
class _Polyhedron_edge_iterator : public It {
protected:
    //It        nt;    // The internal iterator, inherited from It..
public:
    typedef  It  Base;
    typedef  _Polyhedron_edge_iterator<It,Val,Dist,Ctg> Self;

    typedef  Ctg                  iterator_category;
    typedef  Val                  value_type;
    typedef  Dist                 difference_type;

    typedef  value_type&          reference;
    typedef  value_type*          pointer;

// CREATION
// --------

    _Polyhedron_edge_iterator() {}
    _Polyhedron_edge_iterator( It j) : It(j) {}

// OPERATIONS Forward Category
// ---------------------------

    It  current_iterator() const { return It(*this);}

    Self& operator++() {
        It::operator++();
        It::operator++();
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
        It::operator--();
        It::operator--();
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
        It::operator+=( n << 1);
        return *this;
    }
    Self  operator+( difference_type n) const {
        Self tmp = *this;
        return tmp += n;
    }
    Self& operator-=( difference_type n) { return operator+=( -n); }
    Self  operator-( difference_type n) const {
        Self tmp = *this;
        return tmp += -n;
    }
    difference_type  operator-( const Self& i) const {
        return (It::operator-(i)) >> 1;
    }
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
    friend inline  iterator_category
    iterator_category( const Self&) { return iterator_category(); }
    friend inline  value_type*
    value_type( const Self&) { return (value_type*)(0); }
    friend inline  difference_type*
    distance_type( const Self&) { return (difference_type*)(0); }
    friend inline Iterator_tag
    query_circulator_or_iterator( const Self&) {
        return Iterator_tag();
    }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};
template < class D, class It, class Val, class Dist, class Ctg>
inline
_Polyhedron_edge_iterator<It,Val,Dist,Ctg>
operator+( D n, _Polyhedron_edge_iterator<It,Val,Dist,Ctg> i) {
    return i += Dist(n);
}


template < class It, class It2, class Val, class Dist, class Ctg>
class _Polyhedron_edge_const_iterator : public It {
protected:
    //It        nt;    // The internal iterator, inherited from It..
public:
    typedef  It  Base;
    typedef  _Polyhedron_edge_const_iterator<It,It2,Val,Dist,Ctg> Self;

    typedef  Ctg                  iterator_category;
    typedef  Val                  value_type;
    typedef  Dist                 difference_type;

    typedef  const value_type&    reference;
    typedef  const value_type*    pointer;
    typedef  _Polyhedron_edge_iterator<It2,Val,Dist,Ctg> Mutable;

// CREATION
// --------

    _Polyhedron_edge_const_iterator() {}
    _Polyhedron_edge_const_iterator( It j) : It(j) {}
    _Polyhedron_edge_const_iterator( Mutable j)
        : It( j.current_iterator()) {}

// OPERATIONS Forward Category
// ---------------------------

    It  current_iterator() const { return It(*this);}

    Self& operator++() {
        It::operator++();
        It::operator++();
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
        It::operator--();
        It::operator--();
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
        It::operator+=( n << 1);
        return *this;
    }
    Self  operator+( difference_type n) const {
        Self tmp = *this;
        return tmp += n;
    }
    Self& operator-=( difference_type n) { return operator+=( -n); }
    Self  operator-( difference_type n) const {
        Self tmp = *this;
        return tmp += -n;
    }
    difference_type  operator-( const Self& i) const {
        return (It::operator-(i)) >> 1;
    }
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
    friend inline  iterator_category
    iterator_category( const Self&) { return iterator_category(); }
    friend inline  value_type*
    value_type( const Self&) { return (value_type*)(0); }
    friend inline  difference_type*
    distance_type( const Self&) { return (difference_type*)(0); }
    friend inline Iterator_tag
    query_circulator_or_iterator( const Self&) {
        return Iterator_tag();
    }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};
template < class D, class It, class It2, class Val, class Dist, class Ctg>
inline
_Polyhedron_edge_const_iterator<It,It2,Val,Dist,Ctg>
operator+( D n,
           _Polyhedron_edge_const_iterator<It,It2,Val,Dist,Ctg> i) {
    return i += Dist(n);
}
template < class Node, class It, class Ctg>
class _Polyhedron_facet_circ : public It {
    // Ptr      nt;    // The internal node ptr inherited from It.
public:
    typedef  It  Base;
    typedef  _Polyhedron_facet_circ<Node,It,Ctg> Self;

    typedef  Ctg                iterator_category;
    typedef  Node               value_type;
    typedef  std::ptrdiff_t     difference_type;
    typedef  std::size_t        size_type;
    typedef  value_type&        reference;
    typedef  value_type*        pointer;


// CREATION
// --------

    _Polyhedron_facet_circ() : It(0) {}
    //_Polyhedron_facet_circ( pointer p) : It(p) {}
    _Polyhedron_facet_circ( It i) : It(i) {}

// OPERATIONS Forward Category
// ---------------------------

    pointer  ptr() const { return & It::operator*();}

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        return It::operator==( It(NULL));
    }
    bool operator!=( CGAL_NULL_TYPE p) const { return !(*this == p); }
    bool operator==( const Self& i) const { return  It::operator==(i); }
    bool operator!=( const Self& i) const { return !(*this == i); }

    Self& operator++() {
        nt = (*nt).next();
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
        nt = (*nt).prev();
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
    friend inline  iterator_category
    iterator_category( const Self&) { return iterator_category(); }
    friend inline  value_type*
    value_type( const Self&) { return (value_type*)(0); }
    friend inline  difference_type*
    distance_type( const Self&) { return (difference_type*)(0); }
    friend inline  Circulator_tag
    query_circulator_or_iterator( const Self&) {
        return Circulator_tag();
    }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};


template < class Node, class It, class Ctg>
class _Polyhedron_facet_const_circ : public It {
    // Ptr      nt;    // The internal node ptr inherited from It.
public:
    typedef  It  Base;
    typedef  _Polyhedron_facet_const_circ<Node,It,Ctg> Self;

    typedef  Ctg                iterator_category;
    typedef  Node               value_type;
    typedef  std::ptrdiff_t     difference_type;
    typedef  std::size_t        size_type;
    typedef  const value_type&  reference;
    typedef  const value_type*  pointer;

// CREATION
// --------

    _Polyhedron_facet_const_circ() : It(0) {}
    _Polyhedron_facet_const_circ( pointer p) : It(p) {}
    _Polyhedron_facet_const_circ( It i) : It(i) {}

// OPERATIONS Forward Category
// ---------------------------

    pointer  ptr() const { return & It::operator*();}

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        return It::operator==( It(NULL));
    }
    bool operator!=( CGAL_NULL_TYPE p) const { return !(*this == p); }
    bool operator==( const Self& i) const { return  It::operator==(i); }
    bool operator!=( const Self& i) const { return !(*this == i); }

    Self& operator++() {
        nt = (*nt).next();
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
        nt = (*nt).prev();
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
    friend inline  iterator_category
    iterator_category( const Self&) { return iterator_category(); }
    friend inline  value_type*
    value_type( const Self&) { return (value_type*)(0); }
    friend inline  difference_type*
    distance_type( const Self&) { return (difference_type*)(0); }
    friend inline  Circulator_tag
    query_circulator_or_iterator( const Self&) {
        return Circulator_tag();
    }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};
template < class Node, class It, class Ctg>
class _Polyhedron_vertex_circ : public It {
    // Ptr      nt;    // The internal node ptr inherited from It.
public:
    typedef  It  Base;
    typedef  _Polyhedron_vertex_circ<Node,It,Ctg> Self;

    typedef  Ctg                iterator_category;
    typedef  Node               value_type;
    typedef  std::ptrdiff_t     difference_type;
    typedef  std::size_t        size_type;
    typedef  value_type&        reference;
    typedef  value_type*        pointer;

// CREATION
// --------

    _Polyhedron_vertex_circ() : It(0) {}
    //_Polyhedron_vertex_circ( pointer p) : It(p) {}
    _Polyhedron_vertex_circ( It i) : It(i) {}

// OPERATIONS Forward Category
// ---------------------------

    pointer  ptr() const { return & It::operator*();}

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        return It::operator==( It(NULL));
    }
    bool operator!=( CGAL_NULL_TYPE p) const { return !(*this == p); }
    bool operator==( const Self& i) const { return  It::operator==(i); }
    bool operator!=( const Self& i) const { return !(*this == i); }

    Self& operator++() {
        nt = (*nt).next()->opposite();
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
        nt = (*nt).opposite()->prev();
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
    friend inline  iterator_category
    iterator_category( const Self&) { return iterator_category(); }
    friend inline  value_type*
    value_type( const Self&) { return (value_type*)(0); }
    friend inline  difference_type*
    distance_type( const Self&) { return (difference_type*)(0); }
    friend inline  Circulator_tag
    query_circulator_or_iterator( const Self&) {
        return Circulator_tag();
    }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};


template < class Node, class It, class Ctg>
class _Polyhedron_vertex_const_circ : public It {
    // Ptr      nt;    // The internal node ptr inherited from It.
public:
    typedef  It  Base;
    typedef  _Polyhedron_vertex_const_circ<Node,It,Ctg> Self;

    typedef  Ctg                iterator_category;
    typedef  Node               value_type;
    typedef  std::ptrdiff_t     difference_type;
    typedef  std::size_t        size_type;
    typedef  const value_type&  reference;
    typedef  const value_type*  pointer;

// CREATION
// --------

    _Polyhedron_vertex_const_circ() : It(0) {}
    _Polyhedron_vertex_const_circ( pointer p) : It(p) {}
    _Polyhedron_vertex_const_circ( It i) : It(i) {}

// OPERATIONS Forward Category
// ---------------------------

    pointer  ptr() const { return & It::operator*();}

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        return It::operator==( It(NULL));
    }
    bool operator!=( CGAL_NULL_TYPE p) const { return !(*this == p); }
    bool operator==( const Self& i) const { return  It::operator==(i); }
    bool operator!=( const Self& i) const { return !(*this == i); }

    Self& operator++() {
        nt = (*nt).next()->opposite();
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
        nt = (*nt).opposite()->prev();
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
    friend inline  iterator_category
    iterator_category( const Self&) { return iterator_category(); }
    friend inline  value_type*
    value_type( const Self&) { return (value_type*)(0); }
    friend inline  difference_type*
    distance_type( const Self&) { return (difference_type*)(0); }
    friend inline  Circulator_tag
    query_circulator_or_iterator( const Self&) {
        return Circulator_tag();
    }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};

template < class T>
struct Polyhedron_circulator_traits {};

CGAL_TEMPLATE_NULL
struct Polyhedron_circulator_traits<Tag_true> {
    typedef Bidirectional_circulator_tag  iterator_category;
};

CGAL_TEMPLATE_NULL
struct Polyhedron_circulator_traits<Tag_false> {
    typedef Forward_circulator_tag  iterator_category;
};

CGAL_END_NAMESPACE
#endif // CGAL_POLYHEDRON_ITERATOR_3_H //
// EOF //
