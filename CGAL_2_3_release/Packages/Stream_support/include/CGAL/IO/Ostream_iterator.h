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
// release_date  :
//
// file          : IO/Ostream_iterator.h
// package       : Stream_support (2.4)
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// A General Ostream_iterator
// ======================================================================

#ifndef CGAL_IO_OSTREAM_ITERATOR_H
#define CGAL_IO_OSTREAM_ITERATOR_H

#include <CGAL/circulator.h>

CGAL_BEGIN_NAMESPACE

// This proxy is for the Ostream_iterator.
template <class T, class Stream>
class Ostream_proxy {
    Stream& stream;
public:
    Ostream_proxy( Stream& s) : stream(s) {}
    Ostream_proxy<T,Stream>&  operator=( const T& t) {
        stream << t;
        return *this;
    }
};

template <class T, class Stream>
class Ostream_iterator {
    Stream& stream;
public:
    typedef  T                         value_type;
    typedef  T&                        reference;
    typedef  const T&                  const_reference;
    typedef  T*                        pointer;
    typedef  const T*                  const_pointer;
    typedef  std::ptrdiff_t            difference_type;
    typedef  std::output_iterator_tag  iterator_category;

    Ostream_iterator( Stream& s) : stream(s) {}
    Ostream_iterator<T,Stream>& operator++()      { return *this;}
    Ostream_iterator<T,Stream>  operator++(int)   { return *this;}
    Ostream_proxy<T,Stream>     operator*() const {
        return Ostream_proxy<T,Stream>(stream);
    }
};

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
template <class T, class Stream> inline
std::output_iterator_tag
iterator_category( const Ostream_iterator<T,Stream>&) {
    return std::output_iterator_tag();
}
template <class T, class Stream> inline
T*
value_type( const Ostream_iterator<T,Stream>&) {
    return (T*)0;
}
template <class T, class Stream> inline
Iterator_tag
query_circulator_or_iterator(
    const Ostream_iterator<T,Stream>&) {
    return Iterator_tag();
}
#endif

CGAL_END_NAMESPACE

#endif // CGAL_IO_OSTREAM_ITERATOR_H
