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
// file          : IO/Istream_iterator.h
// package       : Stream_support (2.4)
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// A General Istream_iterator
// ======================================================================

#ifndef CGAL_IO_ISTREAM_ITERATOR_H
#define CGAL_IO_ISTREAM_ITERATOR_H

#include <CGAL/circulator.h>

CGAL_BEGIN_NAMESPACE

template <class T, class Stream>
class Istream_iterator {
protected:
    Stream* stream;
    T value;
    void read() {
        if ( stream) {
            if ( *stream) {
                *stream >> value;
                if ( ! *stream)
                    stream = 0;
            } else
                stream = 0;
        }
    }
public:
    typedef  T                           value_type;
    typedef  const T&                    reference;
    typedef  const T&                    const_reference;
    typedef  const T*                    pointer;
    typedef  const T*                    const_pointer;
    typedef  std::size_t                 size_type;
    typedef  std::ptrdiff_t              difference_type;
    typedef  std::input_iterator_tag     iterator_category;
    typedef  Istream_iterator<T,Stream>  Self;

    Istream_iterator() : stream(0) {}
    Istream_iterator( Stream& s) : stream(&s) { read(); }
    bool      operator==( const Self& i) const {
                  return stream == i.stream;
    }
   bool      operator!=( const Self& i) const {
                  return stream != i.stream;
    }

    reference operator*()  const { return value; }
#ifdef  CGAL_ARROW_OPERATOR
    pointer   operator->() const { return &(operator*()); }
#endif
    Self&     operator++()      {
                  read();
                  return *this;
    }
    Self      operator++(int)   {
                  Self tmp = *this;
                  read();
                  return tmp;
    }
};

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
template <class T, class Stream> inline
std::input_iterator_tag
iterator_category( const Istream_iterator<T,Stream>&) {
    return std::input_iterator_tag();
}
template <class T, class Stream> inline
T*
value_type( const Istream_iterator<T,Stream>&) {
    return (T*)0;
}
template <class T, class Stream> inline
std::ptrdiff_t*
distance_type( const Istream_iterator<T,Stream>&) {
    return (std::ptrdiff_t*)0;
}
template <class T, class Stream> inline
Iterator_tag
query_circulator_or_iterator(
    const Istream_iterator<T,Stream>&) {
    return Iterator_tag();
}
#endif

CGAL_END_NAMESPACE

#endif // CGAL_IO_ISTREAM_ITERATOR_H
