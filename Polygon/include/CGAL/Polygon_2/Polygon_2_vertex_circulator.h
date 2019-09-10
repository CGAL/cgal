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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>

#ifndef CGAL_POLYGON_2_VERTEX_CIRCULATOR_H
#define CGAL_POLYGON_2_VERTEX_CIRCULATOR_H

namespace CGAL {

template < class  Ctnr>
class Polygon_circulator {
public:
// TYPES

    typedef Polygon_circulator<Ctnr> Self;
    typedef Circulator_from_container<Ctnr>       Mutable;
    typedef Ctnr                                  Container;
    typedef typename Ctnr::iterator               iterator;
    typedef typename Ctnr::const_iterator         const_iterator;
    typedef typename Ctnr::value_type             value_type;
    typedef typename Ctnr::const_reference        reference;
    typedef const value_type*                     pointer;
    typedef typename Ctnr::size_type              size_type;
    typedef typename Ctnr::difference_type        difference_type;

    typedef std::iterator_traits<const_iterator>  ITraits;
    typedef typename ITraits::iterator_category   Icategory;
    typedef I_Circulator_from_iterator_traits<Icategory> CTraits;
    typedef typename CTraits::iterator_category   iterator_category;

private:
    const Ctnr*    ctnr;
    iterator i;

public:
// CREATION

    Polygon_circulator() : ctnr(NULL) {}
    Polygon_circulator( const Ctnr* c)
        : ctnr(c), i(c->begin()) {}
    Polygon_circulator( const Ctnr* c, iterator j)
        : ctnr(c), i(j) {}
    Polygon_circulator( const Mutable& c)
        : ctnr( c.container()), i( c.current_iterator()) {}

// Gnu-bug workaround: define operator= explicitly.
    Self& operator=( const Self& c) {
        ctnr = c.ctnr;
        i    = c.i;
        return *this;
    }

// OPERATIONS

  bool operator==( Nullptr_t CGAL_assertion_code(p)) const {
        CGAL_assertion( p == NULL);
        return (ctnr == NULL) || (ctnr->begin() == ctnr->end());
    }
    bool operator!=( Nullptr_t p) const { return !(*this == p); }
    bool operator==( const Self& c) const { return i == c.i; }
    bool operator!=( const Self& c) const { return !(*this == c); }
    reference  operator*() const {
        CGAL_assertion( ctnr != NULL);
        CGAL_assertion( current_iterator() != ctnr->end());
        return *i;
    }

private:
// For cases where iterator is a pointer.
    template < typename T >
    static pointer deref(const T& t) { return t.operator->(); }
    template < typename T >
    static pointer deref(T* t) { return t; }

public:

    pointer  operator->() const {
        CGAL_assertion( ctnr != NULL);
        CGAL_assertion( current_iterator() != ctnr->end());
        return deref(i);
    }
    Self& operator++() {
        CGAL_assertion( ctnr != NULL);
        CGAL_assertion( current_iterator() != ctnr->end());
        ++i;
        if ( current_iterator() == ctnr->end())
            i = const_cast<Container*>(ctnr)->begin();
        return *this;
    }
    Self operator++(int) {
        Self tmp= *this;
        ++*this;
        return tmp;
    }
    Self& operator--() {
        CGAL_assertion( ctnr != NULL);
        CGAL_assertion( current_iterator() != ctnr->end());
        if ( current_iterator() == ctnr->begin())
            i = const_cast<Container*>(ctnr)->end();
        --i;
        return *this;
    }
    Self operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }
    Self& operator+=( difference_type n) {
        CGAL_assertion( ctnr != NULL);
        CGAL_assertion( current_iterator() != ctnr->end());
        typename Ctnr::difference_type j = current_iterator() - ctnr->begin();
        typename Ctnr::difference_type size = ctnr->size();
        CGAL_assertion( j    >= 0);
        CGAL_assertion( size >= 0);
        j = non_negative_mod( j + n, size);
        CGAL_assertion( j >= 0);
        CGAL_assertion( j < size);
        i = const_cast<Container*>(ctnr)->begin() + j;
        return *this;
    }
    Self operator+( difference_type n) const {
        Self tmp = *this;
        return tmp += n;
    }
    Self& operator-=( difference_type n) { return operator+=( -n); }
    Self operator-( difference_type n) const {
        Self tmp = *this;
        return tmp += -n;
    }
    difference_type operator-( const Self& c) const {
        CGAL_assertion( ctnr != NULL);
        CGAL_assertion( c.ctnr != NULL);
        return i - c.i;
    }
    reference  operator[]( difference_type n) const {
        Self tmp = *this;
        tmp += n;
        return *tmp;
    }
    const_iterator current_iterator() const { return i;}
    iterator       mod_iterator()     const { return i;}
    Self           min_circulator()   const { return Self(ctnr); }
    const Ctnr*    container()        const { return ctnr; }
};

template <class Ctnr>
inline
Polygon_circulator<Ctnr>
operator+( typename Polygon_circulator<Ctnr>::
               difference_type n,
           const Polygon_circulator<Ctnr>& c) {
    Polygon_circulator<Ctnr> tmp = c;
    return tmp += n;
}

}  // end of namespace CGAL

#endif  // CGAL_POLYGON_2_VERTEX_CIRCULATOR_H
