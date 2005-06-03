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


#ifndef CGAL_CONVERTIBLE_ITERATOR_PROJECT_H
#define CGAL_CONVERTIBLE_ITERATOR_PROJECT_H

#include <CGAL/Iterator_project.h>

CGAL_BEGIN_NAMESPACE


// This class inherits from Iterator_project<> + 
// adds a conversion to handle/const handle.
// See Iterator_project<> documentation.

template<class I,               // internal iterator
         class Fct,             // conversion functor
         class ConstHandle,     // const-handle type to convert to
         class Handle = void*>  // non-const-handle type to convert to
                                // (void* means none)
class Convertible_iterator_project 
    : public Iterator_project<I, Fct>
{
    typedef Iterator_project<I, Fct>        Base;
    typedef Convertible_iterator_project    Self;

public:

  // CREATION
  // --------

    Convertible_iterator_project() {}
    Convertible_iterator_project(Base base) : Base(base) {}

    Convertible_iterator_project(const Self& it) : Base(it) {}
    Self& operator=(const Self& it) { Base::operator=(it); return *this; }

  // OPERATIONS Forward Category
  // ---------------------------

    bool  operator==(CGAL_NULL_TYPE ptr) const { return Base::operator==(ptr); }
    bool  operator!=(CGAL_NULL_TYPE ptr) const { return ! (*this == ptr); }
    bool  operator==(const Self& it) const { return Base::operator==(it); }
    bool  operator!=(const Self& it) const { return ! (*this == it); }

    Self& operator++()     { Base::operator++(); return *this; }
    Self  operator++(int)  { Self tmp(*this); ++(*this); return tmp; }

  // OPERATIONS Bidirectional Category
  // ---------------------------------

    Self& operator--()     { Base::operator--(); return *this; }
    Self  operator--(int)  { Self tmp(*this); --(*this); return tmp; }

  // EXTRA CASTS
  // ---------------------------

    operator Handle()               { return operator->(); }
    operator ConstHandle() const    { return operator->(); }
};


CGAL_END_NAMESPACE

#endif //CGAL_CONVERTIBLE_ITERATOR_PROJECT_H
