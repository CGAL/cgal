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
// file          : Counting_iterator.h
// package       : STL_Extension (2.7)
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
// source        : stl_extension.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// An iterator adaptor for the identity function.
// ======================================================================

#ifndef CGAL_COUNTING_ITERATOR_H
#define CGAL_COUNTING_ITERATOR_H 1
#ifndef CGAL_CIRCULATOR_H
#include <CGAL/circulator.h>
#endif

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
template < class I, class Val>
#else
template < class I,
           class Val = typename std::iterator_traits<I>::value_type>
#endif
class Counting_iterator {
private:
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
    friend inline  value_type*
    value_type( const Self&) { return (value_type*)(0); }
    friend inline  iterator_category
    iterator_category( const Self&){ return iterator_category(); }
    friend inline  difference_type*
    distance_type( const Self&) { return (difference_type*)(0); }
    friend inline  Iterator_tag
    query_circulator_or_iterator( const Self&) { return Iterator_tag(); }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};

CGAL_END_NAMESPACE
#endif // CGAL_COUNTING_ITERATOR_H //
// EOF //
