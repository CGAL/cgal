// Copyright (c) 1997, 2007 
// GeometryFactory (France),
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//                 Fernando Cacciola <fernando.cacciola@geometryfactory.com>
// 

#ifndef CGAL_SAFE_CIRCULATOR_FROM_ITERATOR_H
#define CGAL_SAFE_CIRCULATOR_FROM_ITERATOR_H

#include <CGAL/circulator.h>

#include <boost/optional.hpp>

namespace CGAL {


// Note: Tt, Ss, and Dd are here for backwards compatibility, they are
// not used.
template < class  I, class Tt = int, class Ss = int, class Dd = int>
class Safe_circulator_from_iterator {
public:
// TYPES

    typedef Safe_circulator_from_iterator<I,Tt,Ss,Dd> Self;
    typedef Circulator_from_iterator<I,Tt,Ss,Dd>      Unsafe;
    typedef I                                         iterator;
    typedef std::iterator_traits<iterator>            Traits;

    typedef typename Traits::value_type               value_type;
    typedef std::size_t                               size_type;
    typedef typename Traits::difference_type          difference_type;
    typedef typename Traits::reference                reference;
    typedef typename Traits::pointer                  pointer;

    typedef typename Traits::iterator_category           Icategory;
    typedef I_Circulator_from_iterator_traits<Icategory> CTraits;
    typedef typename CTraits::iterator_category          iterator_category;

private:

    boost::optional<I> m_begin;
    boost::optional<I> m_end;
    boost::optional<I> m_current;
    bool m_empty;

public:
// CREATION

    Safe_circulator_from_iterator() : m_begin(),
                                 m_end(),
                                 m_current(),
				 m_empty( true)
  {}

    Safe_circulator_from_iterator( const I& bgn, const I& end)
        : m_begin(bgn), m_end(end), m_current(bgn), m_empty(bgn==end) {}

    Safe_circulator_from_iterator( const I& bgn, const I& end, const I& cur)
        : m_begin(bgn), m_end(end), m_current(cur), m_empty(bgn==end) {}

    Safe_circulator_from_iterator( const Self& c, const I& cur)
        : m_begin( c.m_begin), m_end( c.m_end), m_current(cur), m_empty(c.m_empty) {}


    template <class II, class A1, class A2, class A3>
    // Allow construction from Circulator_from_iterator with
    // assignment compatible iterator II:
    Safe_circulator_from_iterator(
        const Safe_circulator_from_iterator<II,A1,A2,A3>& ii)
    : m_begin( ii.m_begin), m_end( ii.m_end),
        m_current(ii.m_current), m_empty(ii.m_empty) {}

//
// OPERATIONS

    bool operator==( Nullptr_t p) const {
        CGAL_assertion( p == NULL);
        return m_empty;
    }
    bool operator!=( Nullptr_t p) const { return !(*this == p); }
    
    bool operator==( const Self& c) const
    {
      CGAL_assertion( is_not_singular() ) ;
      return current_iterator() == c.current_iterator();
    }
    
    bool operator!=( const Self& c) const { return !(*this == c); }
    
    reference  operator*() const {
        CGAL_assertion( is_not_singular() ) ;
        CGAL_assertion( current_iterator() != end() );
        return *current_iterator();
    }
    pointer  operator->() const {
        CGAL_assertion( is_not_singular() ) ;
        CGAL_assertion( current_iterator() != end() );
        return &(*current_iterator());
    }
    Self& operator++() {
        CGAL_assertion( is_not_singular() ) ;
        CGAL_assertion( current_iterator() != end() );
        ++current_iterator();
        if ( current_iterator() == end() )
            current_iterator() = begin();
        return *this;
    }
    Self  operator++(int) {
        Self tmp= *this;
        ++*this;
        return tmp;
    }
    Self& operator--() {
        CGAL_assertion( is_not_singular() ) ;
        CGAL_assertion( current_iterator() != end() );
        if ( current_iterator() == begin() )
            current_iterator() = end();
        --current_iterator();
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }
    Self& operator+=( difference_type n) {
        CGAL_assertion( is_not_singular() ) ;
        CGAL_assertion( current_iterator() != end() );
        difference_type i    = current_iterator() - begin();
        difference_type size = end()              - begin();
        CGAL_assertion( i    >= 0);
        CGAL_assertion( size >= 0);
        i = non_negative_mod( i + n, size);
        CGAL_assertion( i >= 0);
        CGAL_assertion( i < size);
        current_iterator() = begin() + i;
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
        CGAL_assertion((begin() == i.begin()) && (end() == i.end()));
        return current_iterator() - i.current_iterator();
    }
    reference  operator[](difference_type n) const {
        Self tmp = *this;
        tmp += n;
        return tmp.operator*();
    }
    
    iterator const& begin()            const { return *m_begin;}
    iterator const& end()              const { return *m_end;}
    iterator const& current_iterator() const { return *m_current;}
    Self            min_circulator()   const { return Self( m_begin, m_end, m_begin, m_empty); }
    
    Unsafe unsafe_circulator() const
    {
      CGAL_assertion( is_not_singular() ) ;
      
      return Unsafe(begin(),end(),current_iterator());
    }
    
private:

    bool is_not_singular() const { return !!m_begin && !!m_end && !!m_current ; }
    
    iterator& begin()            { return *m_begin;}
    iterator& end()              { return *m_end;}
    iterator& current_iterator() { return *m_current;}


};

//template < class I, class  T, class Size, class Dist>
//I Circulator_from_iterator< I, T, Size, Dist>::null_iterator = I();

template < class D, class I, class  T, class Size, class Dist> inline
Safe_circulator_from_iterator< I, T, Size, Dist>
operator+( D n, const
    Safe_circulator_from_iterator< I, T, Size, Dist>& circ) {
    Safe_circulator_from_iterator< I, T, Size, Dist>
        tmp = circ;
    return tmp += Dist(n);
}

} //namespace CGAL

#endif // CGAL_CIRCULATOR_H //
// EOF //
