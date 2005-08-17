// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Laurent Saboret, Pierre Alliez


#ifndef CGAL_CONVERTIBLE_FILTER_ITERATOR_H
#define CGAL_CONVERTIBLE_FILTER_ITERATOR_H

#include <CGAL/iterator.h>

CGAL_BEGIN_NAMESPACE


// This class inherits from Filter_iterator<> +
// adds a conversion to handle/const handle.
// See Filter_iterator<> documentation.

template<class Iterator,        // internal iterator
         class Predicate,       // predicate to filter out elements
         class ConstHandle,     // const-handle type to convert to
         class Handle = void*>  // non-const-handle type to convert to
                                // (void* means none)
class Convertible_filter_iterator
    : public Filter_iterator<Iterator, Predicate>
{
    typedef Filter_iterator<Iterator, Predicate>    Base;
    typedef Convertible_filter_iterator             Self;

public:

  // CREATION
  // --------

    Convertible_filter_iterator() {}
    Convertible_filter_iterator(Base base)
        : Base(base) {}
    Convertible_filter_iterator(Iterator e, const Predicate& p)
        : Base(e,p) {}
    Convertible_filter_iterator(Iterator e, const Predicate& p, Iterator c)
        : Base(e,p,c) {}

    Convertible_filter_iterator(const Self& it) : Base(it) {}
    Self& operator=(const Self& it) { Base::operator=(it); return *this; }

  // OPERATIONS Forward Category
  // ---------------------------

    bool  operator==(CGAL_NULL_TYPE ptr) const { return Base::operator==(ptr); }
    bool  operator!=(CGAL_NULL_TYPE ptr) const { return ! (*this == ptr); }
    bool  operator==(const Self& it) const { return (Base&)*this == (Base&)it; }
    bool  operator!=(const Self& it) const { return ! (*this == it); }

    Self& operator++()     { Base::operator++(); return *this; }
    Self  operator++(int)  { Self tmp(*this); ++(*this); return tmp; }

  // OPERATIONS Bidirectional Category
  // ---------------------------------

    Self& operator--()     { Base::operator--(); return *this; }
    Self  operator--(int)  { Self tmp(*this); --(*this); return tmp; }

  // EXTRA CASTS
  // ---------------------------

    operator Handle()               { return Base::operator->(); }
    operator ConstHandle() const    { return Base::operator->(); }
};


CGAL_END_NAMESPACE

#endif //CGAL_CONVERTIBLE_FILTER_ITERATOR_H
