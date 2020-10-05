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

    Polygon_circulator() : ctnr(nullptr) {}
    Polygon_circulator( const Ctnr* c)
        : ctnr(c), i(c->begin()) {}
    Polygon_circulator( const Ctnr* c, iterator j)
        : ctnr(c), i(j) {}
    Polygon_circulator( const Mutable& c)
        : ctnr( c.container()), i( c.current_iterator()) {}



// OPERATIONS

  bool operator==( std::nullptr_t CGAL_assertion_code(p)) const {
        CGAL_assertion( p == nullptr);
        return (ctnr == nullptr) || (ctnr->begin() == ctnr->end());
    }
    bool operator!=( std::nullptr_t p) const { return !(*this == p); }
    bool operator==( const Self& c) const { return i == c.i; }
    bool operator!=( const Self& c) const { return !(*this == c); }
    reference  operator*() const {
        CGAL_assertion( ctnr != nullptr);
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
        CGAL_assertion( ctnr != nullptr);
        CGAL_assertion( current_iterator() != ctnr->end());
        return deref(i);
    }
    Self& operator++() {
        CGAL_assertion( ctnr != nullptr);
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
        CGAL_assertion( ctnr != nullptr);
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
        CGAL_assertion( ctnr != nullptr);
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
        CGAL_assertion( ctnr != nullptr);
        CGAL_assertion( c.ctnr != nullptr);
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
