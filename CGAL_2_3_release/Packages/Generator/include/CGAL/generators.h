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
// file          : generators.h
// chapter       : $CGAL_Chapter: Geometric Object Generators $
// package       : $CGAL_Package: Generator 2.12 (28 Jul 1999) $
// source        : generators.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// General Generator Support
// ============================================================================

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
