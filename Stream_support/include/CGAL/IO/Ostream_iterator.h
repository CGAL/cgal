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

#ifndef CGAL_IO_OSTREAM_ITERATOR_H
#define CGAL_IO_OSTREAM_ITERATOR_H

#include <CGAL/circulator.h>

namespace CGAL {

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
    Stream* stream;
public:
    typedef  T                         value_type;
    typedef  T&                        reference;
    typedef  const T&                  const_reference;
    typedef  T*                        pointer;
    typedef  const T*                  const_pointer;
    typedef  std::ptrdiff_t            difference_type;
    typedef  std::output_iterator_tag  iterator_category;

    Ostream_iterator( Stream& s) : stream(&s) {}
    Ostream_iterator<T,Stream>& operator++()      { return *this;}
    Ostream_iterator<T,Stream>  operator++(int)   { return *this;}
    Ostream_proxy<T,Stream>     operator*() const {
        return Ostream_proxy<T,Stream>(*stream);
    }
};

} //namespace CGAL

#endif // CGAL_IO_OSTREAM_ITERATOR_H
