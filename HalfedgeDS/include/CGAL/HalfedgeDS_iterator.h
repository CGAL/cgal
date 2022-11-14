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
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>)

#ifndef CGAL_HALFEDGEDS_ITERATOR_H
#define CGAL_HALFEDGEDS_ITERATOR_H 1

#include <CGAL/circulator.h>
#include <CGAL/tags.h>
#include <CGAL/use.h>

namespace CGAL {



template < class It, class Ctg>
class I_HalfedgeDS_facet_circ : public It {
public:
    typedef  It                                 Iterator;
    typedef  Ctg                                iterator_category;
    typedef  I_HalfedgeDS_facet_circ<It,Ctg>    Self;
    typedef  std::iterator_traits<It>           Traits;
    typedef  typename Traits::value_type        value_type;
    typedef  typename Traits::difference_type   difference_type;
    typedef  std::size_t                        size_type;
    typedef  typename Traits::reference         reference;
    typedef  typename Traits::pointer           pointer;

// CREATION
// --------

    I_HalfedgeDS_facet_circ() {}
    explicit I_HalfedgeDS_facet_circ( It i) : It(i) {}
    template <class It2, class Ctg2>
    I_HalfedgeDS_facet_circ( const I_HalfedgeDS_facet_circ<It2,Ctg2> &c)
        : It((const It2&)(c)) {}

// OPERATIONS Forward Category
// ---------------------------

    // pointer  ptr() const { return & It::operator*();}

    bool operator==( std::nullptr_t CGAL_assertion_code(p)) const {
        CGAL_assertion( p == 0);
        return It::operator==( It());
    }
    bool operator!=( std::nullptr_t p) const { return !(*this == p); }
    bool operator==( const Self& i)    const { return  It::operator==(i); }
    bool operator!=( const Self& i)    const { return !(*this == i); }

    // operator* and operator-> are inherited.

    Self& operator++() {
        *((Iterator*)this) = (*this)->next();
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
        *((Iterator*)this) = (*this)->prev();
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }
};
template < class It, class Ctg>
class I_HalfedgeDS_vertex_circ : public It {
public:
    typedef  It                                  Iterator;
    typedef  Ctg                                 iterator_category;
    typedef  I_HalfedgeDS_vertex_circ<It,Ctg>    Self;
    typedef  std::iterator_traits<It>            Traits;
    typedef  typename Traits::value_type         value_type;
    typedef  typename Traits::difference_type    difference_type;
    typedef  std::size_t                         size_type;
    typedef  typename Traits::reference          reference;
    typedef  typename Traits::pointer            pointer;

// CREATION
// --------

    I_HalfedgeDS_vertex_circ() {}
    explicit I_HalfedgeDS_vertex_circ( It i) : It(i) {}
    template <class It2, class Ctg2>
    I_HalfedgeDS_vertex_circ( const I_HalfedgeDS_vertex_circ<It2,Ctg2> &c)
        : It((const It2&)(c)) {}

// OPERATIONS Forward Category
// ---------------------------

    // pointer  ptr() const { return & It::operator*();}

    bool operator==( std::nullptr_t CGAL_assertion_code(p)) const {
        CGAL_assertion( p == 0);
        return It::operator==( It());
    }
    bool operator!=( std::nullptr_t p) const { return !(*this == p); }
    bool operator==( const Self& i)    const { return  It::operator==(i); }
    bool operator!=( const Self& i)    const { return !(*this == i); }

    // operator* and operator-> are inherited.

    Self& operator++() {
        *((Iterator*)this) = (*this)->next()->opposite();
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
        *((Iterator*)this) = (*this)->opposite()->prev();
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }
};

// Determine the circulator category: If prev() is supported,
// its bidirectional, otherwise its forward.

template < class T> struct HalfedgeDS_circulator_traits {};

template <>
struct HalfedgeDS_circulator_traits<Tag_true> {
    typedef Bidirectional_circulator_tag  iterator_category;
};

template <>
struct HalfedgeDS_circulator_traits<Tag_false> {
    typedef Forward_circulator_tag  iterator_category;
};

// portability with other code using the old HalfedgeDS circulators

template < class Node, class It, class Ctg>
class _HalfedgeDS_facet_circ : public It {
    // Ptr      nt;    // The internal node ptr inherited from It.
public:
    typedef  It  Base;
    typedef  _HalfedgeDS_facet_circ<Node,It,Ctg> Self;

    typedef  Ctg                iterator_category;
    typedef  Node               value_type;
    typedef  std::ptrdiff_t     difference_type;
    typedef  std::size_t        size_type;
    typedef  value_type&        reference;
    typedef  value_type*        pointer;


// CREATION
// --------

    _HalfedgeDS_facet_circ() : It(nullptr) {}
    //_HalfedgeDS_facet_circ( pointer p) : It(p) {}
    _HalfedgeDS_facet_circ( It i) : It(i) {}

// OPERATIONS Forward Category
// ---------------------------

    pointer  ptr() const { return & It::operator*();}

    bool operator==( std::nullptr_t p) const {
        CGAL_USE(p);
        CGAL_assertion( p == nullptr);
        return It::operator==( It(nullptr));
    }
    bool operator!=( std::nullptr_t p) const { return !(*this == p); }
    bool operator==( const Self& i) const { return  It::operator==(i); }
    bool operator!=( const Self& i) const { return !(*this == i); }

    Self& operator++() {
        this->nt = typename It::Iterator((*this->nt).next());
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
        this->nt = typename It::Iterator((*this->nt).prev());
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }
};


template < class Node, class It, class Ctg>
class _HalfedgeDS_facet_const_circ : public It {
    // Ptr      nt;    // The internal node ptr inherited from It.
public:
    typedef  It  Base;
    typedef  _HalfedgeDS_facet_const_circ<Node,It,Ctg> Self;

    typedef  Ctg                iterator_category;
    typedef  Node               value_type;
    typedef  std::ptrdiff_t     difference_type;
    typedef  std::size_t        size_type;
    typedef  const value_type&  reference;
    typedef  const value_type*  pointer;

// CREATION
// --------

    _HalfedgeDS_facet_const_circ() : It(nullptr) {}
    _HalfedgeDS_facet_const_circ( pointer p) : It(p) {}
    _HalfedgeDS_facet_const_circ( It i) : It(i) {}

    template <class NN, class II, class CTG>
    _HalfedgeDS_facet_const_circ( const _HalfedgeDS_facet_circ<NN,II,CTG>& c)
        : It(c.ptr()) {}

// OPERATIONS Forward Category
// ---------------------------

    pointer  ptr() const { return & It::operator*();}

    bool operator==( std::nullptr_t p) const {
        CGAL_USE(p);
        CGAL_assertion( p == nullptr);
        return It::operator==( It(nullptr));
    }
    bool operator!=( std::nullptr_t p) const { return !(*this == p); }
    bool operator==( const Self& i) const { return  It::operator==(i); }
    bool operator!=( const Self& i) const { return !(*this == i); }
    bool operator==( const It& i) const { return  It::operator==(i); }
    bool operator!=( const It& i) const { return !(*this == i); }

    Self& operator++() {
        this->nt = typename It::Iterator((*this->nt).next());
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
        this->nt = typename It::Iterator((*this->nt).prev());
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }
};

template < class Node, class It, class Ctg>
class _HalfedgeDS_vertex_circ : public It {
    // Ptr      nt;    // The internal node ptr inherited from It.
public:
    typedef  It  Base;
    typedef  _HalfedgeDS_vertex_circ<Node,It,Ctg> Self;

    typedef  Ctg                iterator_category;
    typedef  Node               value_type;
    typedef  std::ptrdiff_t     difference_type;
    typedef  std::size_t        size_type;
    typedef  value_type&        reference;
    typedef  value_type*        pointer;

// CREATION
// --------

    _HalfedgeDS_vertex_circ() : It(nullptr) {}
    //_HalfedgeDS_vertex_circ( pointer p) : It(p) {}
    _HalfedgeDS_vertex_circ( It i) : It(i) {}

// OPERATIONS Forward Category
// ---------------------------

    pointer  ptr() const { return & It::operator*();}

    bool operator==( std::nullptr_t p) const {
        CGAL_USE(p);
        CGAL_assertion( p == nullptr);
        return It::operator==( It(nullptr));
    }
    bool operator!=( std::nullptr_t p) const { return !(*this == p); }
    bool operator==( const It& i) const { return  It::operator==(i); }
    bool operator==( const Self& i) const { return  It::operator==(i); }
    bool operator!=( const Self& i) const { return !(*this == i); }

    Self& operator++() {
        this->nt = typename It::Iterator((*this->nt).next()->opposite());
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
        this->nt = typename It::Iterator((*this->nt).opposite()->prev());
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }
};


template < class Node, class It, class Ctg>
class _HalfedgeDS_vertex_const_circ : public It {
    // Ptr      nt;    // The internal node ptr inherited from It.
public:
    typedef  It  Base;
    typedef  _HalfedgeDS_vertex_const_circ<Node,It,Ctg> Self;

    typedef  Ctg                iterator_category;
    typedef  Node               value_type;
    typedef  std::ptrdiff_t     difference_type;
    typedef  std::size_t        size_type;
    typedef  const value_type&  reference;
    typedef  const value_type*  pointer;

// CREATION
// --------

    _HalfedgeDS_vertex_const_circ() : It(nullptr) {}
    _HalfedgeDS_vertex_const_circ( pointer p) : It(p) {}
    _HalfedgeDS_vertex_const_circ( It i) : It(i) {}

    template <class NN, class II, class CTG>
    _HalfedgeDS_vertex_const_circ( const _HalfedgeDS_vertex_circ<NN,II,CTG>& c)
        : It(c.ptr()) {}

// OPERATIONS Forward Category
// ---------------------------

    pointer  ptr() const { return & It::operator*();}

    bool operator==( std::nullptr_t p) const {
        CGAL_USE(p);
        CGAL_assertion( p == nullptr);
        return It::operator==( It(nullptr));
    }
    bool operator!=( std::nullptr_t p) const { return !(*this == p); }
    bool operator==( const Self& i) const { return  It::operator==(i); }
    bool operator!=( const Self& i) const { return !(*this == i); }

    Self& operator++() {
      this->nt = typename It::Iterator((*this->nt).next()->opposite());
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
        this->nt = typename It::Iterator((*this->nt).opposite()->prev());
        return *this;
    }
    Self  operator--(int) {
        Self tmp = *this;
        --*this;
        return tmp;
    }
};

} //namespace CGAL
#endif // CGAL_HALFEDGEDS_ITERATOR_H //
// EOF //
