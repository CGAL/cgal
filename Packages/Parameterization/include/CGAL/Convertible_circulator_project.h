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


#ifndef CGAL_CONVERTIBLE_CIRCULATOR_PROJECT_H
#define CGAL_CONVERTIBLE_CIRCULATOR_PROJECT_H

#include <CGAL/Circulator_project.h>

CGAL_BEGIN_NAMESPACE


// This class inherits from Circulator_project<> + 
// adds a conversion to handle/const handle.
// See Circulator_project<> documentation.

template<class C,               // internal circulator
         class Fct,             // conversion functor
         class Ref,
         class Ptr,
         class ConstHandle,     // const-handle type to convert to
         class Handle = void*>  // non-const-handle type to convert to
                                // (void* means none)
class Convertible_circulator_project 
    : public Circulator_project<C, Fct, Ref, Ptr>
{
    typedef Circulator_project<C, Fct, Ref, Ptr>    Base;
    typedef Convertible_circulator_project          Self;

public:

  // CREATION
  // --------

    Convertible_circulator_project() {}
    Convertible_circulator_project(Base base) : Base(base) {}

    Convertible_circulator_project(const Self& cir) : Base(cir) {}
    Self& operator=(const Self& cir) { Base::operator=(cir); return *this; }

  // OPERATIONS Forward Category
  // ---------------------------

    bool  operator==(CGAL_NULL_TYPE ptr) const { return Base::operator==(ptr); }
    bool  operator!=(CGAL_NULL_TYPE ptr) const { return ! (*this == ptr); }
    bool  operator==(const Self& cir)    const { return Base::operator==(cir); }
    bool  operator!=(const Self& cir)    const { return ! (*this == cir); }

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

#endif //CGAL_CONVERTIBLE_CIRCULATOR_PROJECT_H
