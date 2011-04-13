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
// release_date  : 1999, July 28
//
// file          : Circulator_project.h
// package       : STL_Extension (2.7)
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
// source        : stl_extension.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// An circulator adaptor performing a projection on the value type.
// ======================================================================

#ifndef CGAL_CIRCULATOR_PROJECT_H
#define CGAL_CIRCULATOR_PROJECT_H 1
#ifndef CGAL_CIRCULATOR_H
#include <CGAL/circulator.h>
#endif

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_CFG_NO_DEFAULT_PREVIOUS_TEMPLATE_ARGUMENTS
template < class C, class Fct, class Ref, class Ptr>
#else
template < class C,
           class Fct,
           class Ref = typename C::reference,
           class Ptr = typename C::pointer>
#endif
class Circulator_project {
protected:
    C        nt;    // The internal circulator.
public:
    typedef C  Circulator;
    typedef Circulator_project<C,Fct,Ref,Ptr> Self;

    typedef typename  C::iterator_category   iterator_category;
    typedef typename  Fct::result_type       value_type;
    typedef typename  C::difference_type     difference_type;
    typedef typename  C::size_type           size_type;
    typedef           Ref                    reference;
    typedef           Ptr                    pointer;

// CREATION
// --------

    Circulator_project() {}
    Circulator_project( Circulator j) : nt(j) {}

// OPERATIONS Forward Category
// ---------------------------

    Circulator  current_circulator() const { return nt;}
    Ptr ptr() const {
        Fct fct;
        return &(fct(*nt));
    }

    bool operator==( CGAL_NULL_TYPE p) const {
        CGAL_assertion( p == NULL);
        return ( nt == NULL);
    }
    bool  operator!=( CGAL_NULL_TYPE p) const { return !(*this == p); }
    bool  operator==( const Self& i) const { return ( nt == i.nt); }
    bool  operator!=( const Self& i) const { return !(*this == i); }
    Ref   operator*()  const { return *ptr(); }
    Ptr   operator->() const { return ptr(); }
    Self& operator++() {
        ++nt;
        return *this;
    }
    Self  operator++(int) {
        Self tmp = *this;
        ++*this;
        return tmp;
    }

// OPERATIONS Bidirectional Category
// ---------------------------------

    Self& operator--() {
        --nt;
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }

// OPERATIONS Random Access Category
// ---------------------------------

    Self  min_circulator() const {
        return Self( nt.min_circulator());
    }
    Self& operator+=( difference_type n) {
        nt += n;
        return *this;
    }
    Self  operator+( difference_type n) const {
        Self tmp = *this;
        return tmp += n;
    }
    Self& operator-=( difference_type n) {
        return operator+=( -n);
    }
    Self  operator-( difference_type n) const {
        Self tmp = *this;
        return tmp += -n;
    }
    difference_type  operator-( const Self& i) const {
        return nt - i.nt;
    }
    Ref  operator[]( difference_type n) const {
        Self tmp = *this;
        tmp += n;
        return tmp.operator*();
    }
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
    friend inline  iterator_category
    iterator_category( const Self&) { return iterator_category(); }
    friend inline  value_type*
    value_type( const Self&) { return (value_type*)(0); }
    friend inline  difference_type*
    distance_type( const Self&) { return (difference_type*)(0); }
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //
};

template < class Dist, class Fct, class C, class Ref, class Ptr>
inline
Circulator_project<C,Fct,Ref,Ptr>
operator+( Dist n, Circulator_project<C,Fct,Ref,Ptr> i) { return i += n; }

#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
template < class C, class Fct, class Ref, class Ptr>
inline  Circulator_tag
query_circulator_or_iterator(
    const Circulator_project<C,Fct,Ref,Ptr>&) {
    return Circulator_tag();
}
#endif // CGAL_CFG_NO_ITERATOR_TRAITS //

CGAL_END_NAMESPACE
#endif // CGAL_CIRCULATOR_PROJECT_H //
// EOF //
