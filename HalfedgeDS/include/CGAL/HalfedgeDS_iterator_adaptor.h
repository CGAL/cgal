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
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_HALFEDGEDS_ITERATOR_ADAPTOR_H
#define CGAL_HALFEDGEDS_ITERATOR_ADAPTOR_H 1

#include <iterator>

namespace CGAL {

//  The iterator identity adaptor will be used for the HDS implementations
//  that are based on STL (or other) container classes which do not
//  guarantee that the default construction of its iterator gives always
//  the same singular value (i.e. something like nullptr). This adaptor
//  assumes that iterator traits are fully supported. It works for all
//  kinds of iterators, from input iterators to random access iterators.
//  It does no longer require that these iterators have a constructor
//  accepting zero as an argument (e.g. a pointer). This approach has
//  failed with the KCC compiler (Kai's C++).

//  Instead, we rely now on a static local variable. Static variables are
//  first of all zero-initialized (Section 3.6.2), which guarantees that
//  pointers and such are set to zero even if the construtor does not
//  initialize them (Section 8.5). With static variables, the order of
//  initialization could be critical, if the initialization of one
//  requires another one to be initialized already (I have not seen such a
//  design for iterator yet, but who knows, maybe if reference counting
//  come into play). If global static variables are in different
//  compilation units, there order of initialization is unspecified.
//  However, if the static variable is a local variable, the
//  initialization happens before the surrounding function is called if
//  the type is a POD, otherwise the initialization happens when the
//  control flow passes through the declaration of the static variable the
//  first time (Section 6.7), which costs performance but makes it safe
//  for weird static initialization situations. Usually the std::vector
//  class uses a plain C-pointer as iterator, which would be a POD and
//  thus efficient. However, the std::list iterators might not be POD's if
//  they define their own copy contructor. This is the case for
//  std::list::iterator of the current SGI STL, but not for the
//  std::list::const_iterator, which is a funny side-effect of having
//  only a single class for both and a constructor that allows iterator to
//  const_iterator conversions. The bottom line is, std::vector is fine,
//  std::list could cost a bit performance and In_place_list might be
//  preferable.


template < class I>
class HalfedgeDS_iterator_adaptor {
private:
    I        nt;    // The internal iterator.

    // keep the local static variable for the null iterator
    // in a small inline function, optimizers might like that
    // better than a static function.
    const I& null_iterator() const {
        static const I it = I(); // zero-initialized and a default constructor
        return it;
    }
public:
    typedef I                                  Iterator;
    typedef HalfedgeDS_iterator_adaptor<I>     Self;
    typedef std::iterator_traits<I>            Traits;
    typedef typename Traits::reference         reference;
    typedef typename Traits::pointer           pointer;
    typedef typename Traits::value_type        value_type;
    typedef typename Traits::difference_type   difference_type;
    typedef typename Traits::iterator_category iterator_category;

// CREATION
// --------
                                              // explicitly set to 0
    HalfedgeDS_iterator_adaptor() : nt( null_iterator()) {}
    HalfedgeDS_iterator_adaptor( Iterator j) : nt(j) {} // down cast

    template < class J>
    HalfedgeDS_iterator_adaptor( const HalfedgeDS_iterator_adaptor<J>& j)
        : nt(j.iterator()) {}

// OPERATIONS Forward Category
// ---------------------------

    Iterator  iterator() const                 { return nt;}

    bool      operator==( const Self& i) const { return ( nt == i.nt); }
    bool      operator!=( const Self& i) const { return !(*this == i); }
    reference operator*()  const               { return *nt; }
    pointer   operator->() const               { return &*nt; }
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

    Self& operator+=( difference_type n) {
        nt += n;
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
    difference_type  operator-( const Self& i) const { return nt - i.nt; }
    reference  operator[]( difference_type n) const {
        Self tmp = *this;
        tmp += n;
        return tmp.operator*();
    }
    bool operator< ( const Self& i) const { return ( nt < i.nt); }
    bool operator> ( const Self& i) const { return i < *this;    }
    bool operator<=( const Self& i) const { return !(i < *this); }
    bool operator>=( const Self& i) const { return !(*this < i); }
};

} //namespace CGAL

// we don't need Koenig lookup here
template < class D, class I>
inline
CGAL::HalfedgeDS_iterator_adaptor<I>
operator+( D n, CGAL::HalfedgeDS_iterator_adaptor<I> i) {
    return i += Dist(n);
}
#endif // CGAL_HALFEDGEDS_ITERATOR_ADAPTOR_H //
// EOF //
