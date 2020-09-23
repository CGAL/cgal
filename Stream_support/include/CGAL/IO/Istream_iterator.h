// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>

#ifndef CGAL_IO_ISTREAM_ITERATOR_H
#define CGAL_IO_ISTREAM_ITERATOR_H

#include <CGAL/circulator.h>

namespace CGAL {

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

} //namespace CGAL

#endif // CGAL_IO_ISTREAM_ITERATOR_H
