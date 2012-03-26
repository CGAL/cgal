// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
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

namespace CGAL {
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

class Random_double_in_interval : public Random_generator_base<double> {

 public:
  typedef Random_double_in_interval This;
  Random_double_in_interval(double a = 1, Random& rnd = default_random)
    // g is an input iterator creating points of type `P' uniformly
    // distributed in the half-open square with side length a,
    // centered around the origin, i.e. \forall p = `*g': -\frac{a}{2}
    // <= p.x() < \frac{a}{2} and -\frac{a}{2} <= p.y() < \frac{a}{2}
    // . Two random numbers are needed from `rnd' for each point.
    : Random_generator_base<double>( a, rnd)
    { 
      this->d_item = this->d_range * (2 * this->_rnd.get_double() - 1.0);
    }
  
  This& operator++() {
    this->d_item = this->d_range * (2 * this->_rnd.get_double() - 1.0);
    return *this;
  }
  This  operator++(int) {
    This tmp = *this;
    ++(*this);
    return tmp;
  }
};

} //namespace CGAL
#endif // CGAL_GENERATORS_H //
// EOF //
