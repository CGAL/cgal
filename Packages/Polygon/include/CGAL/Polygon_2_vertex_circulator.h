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
// release       :
// release_date  :
//
// file          : include/CGAL/Polygon_2_vertex_circulator.h
// source        : 
// author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>
//
// coordinator   : Utrecht University
//
// ============================================================================

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

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == CGAL_CIRC_NULL);
        return (ctnr == NULL) || (ctnr->begin() == ctnr->end());
    }
    bool operator!=( CGAL_NULL_TYPE p) const { return !(*this == p); }
    bool operator==( const Self& c) const { return i == c.i; }
    bool operator!=( const Self& c) const { return !(*this == c); }
    reference  operator*() const {
        CGAL_assertion( ctnr != NULL);
        CGAL_assertion( current_iterator() != ctnr->end());
        return *i;
    }
    pointer  operator->() const {
        CGAL_assertion( ctnr != NULL);
        CGAL_assertion( current_iterator() != ctnr->end());
        return i.operator->();
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
#ifndef _MSC_VER
    Self& operator+=( difference_type n);
#else
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
#endif
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

#ifndef _MSC_VER
template <class Ctnr>
Polygon_circulator<Ctnr>&
Polygon_circulator<Ctnr>::
operator+=( typename Ctnr::difference_type n) {
    CGAL_assertion( ctnr != NULL);
    CGAL_assertion( current_iterator() != ctnr->end());
    typename Ctnr::difference_type j    = current_iterator() - ctnr->begin();
    typename Ctnr::difference_type size = ctnr->size();
    CGAL_assertion( j    >= 0);
    CGAL_assertion( size >= 0);
    j = non_negative_mod( j + n, size);
    CGAL_assertion( j >= 0);
    CGAL_assertion( j < size);
    i = const_cast<Container*>(ctnr)->begin() + j;
    return *this;
}
#endif // _MSC_VER


}  // end of namespace CGAL

#endif  // CGAL_POLYGON_2_VERTEX_CIRCULATOR_H
