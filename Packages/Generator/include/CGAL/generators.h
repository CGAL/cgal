// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>

#ifndef CGAL_GENERATORS_H
#define CGAL_GENERATORS_H 1

#include <CGAL/basic.h>
#include <cstddef>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <CGAL/function_objects.h>
#include <CGAL/Random.h>

CGAL_BEGIN_NAMESPACE
template < class T >
class Generator_base {
protected:
    T       d_item;
    double  d_range;
public:
    typedef std::input_iterator_tag iterator_category;
    typedef T                       value_type;
    typedef std::ptrdiff_t          difference_type;
    typedef const T*                pointer;
    typedef const T&                reference;
    typedef Generator_base<T>      This;

    Generator_base() {}
    Generator_base( double range) : d_range( range) {}
    Generator_base( const T& item, double range)
        : d_item(item), d_range(range) {}

    bool operator==( const This& base) const {
        return ( d_item == base.d_item);
    }
    bool operator!=( const This& base) const { return ! operator==(base);}
    double    range()      const { return d_range; }
    reference operator*()  const { return d_item; }
    pointer   operator->() const { return & operator*(); }
};

template < class T >
class Random_generator_base : public Generator_base<T> {
protected:
    Random& _rnd;
public:
    typedef  Random_generator_base<T> This;

    Random_generator_base() : _rnd( default_random) {}
    Random_generator_base( double range, Random& rnd)
        : Generator_base<T>( range), _rnd( rnd) {}
    Random_generator_base( const T& item, double range, Random& rnd)
        : Generator_base<T>( item, range), _rnd( rnd) {}
    bool operator==( const This& rb) const {
        return ( _rnd == rb._rnd && Generator_base<T>::operator==(rb));
    }
    bool operator!=( const This& rb) const { return ! operator==(rb); }
};
CGAL_END_NAMESPACE
#endif // CGAL_GENERATORS_H //
// EOF //
