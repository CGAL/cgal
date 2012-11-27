// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy


#ifndef CGAL_CONVERTIBLE_ITERATOR_PROJECT_H
#define CGAL_CONVERTIBLE_ITERATOR_PROJECT_H

#include <CGAL/Iterator_project.h>

namespace CGAL {

/// \cond SKIP_IN_MANUAL

/// This class inherits from Iterator_project<> +
/// adds a conversion to handle/const handle.
/// See Iterator_project<> documentation.

template<class I,               ///< Internal iterator.
         class Fct,             ///< Conversion functor.
         class ConstHandle,     ///< Const-handle type to convert to.
         class Handle = void*   ///< Non-const-handle type to convert to (void*=none).
>
class Convertible_iterator_project
    : public Iterator_project<I, Fct>
{
    typedef Iterator_project<I, Fct>        Base;
    typedef Convertible_iterator_project    Self;

public:

  /// \name Creation
  /// @{
    Convertible_iterator_project() {}
    Convertible_iterator_project(Base base) : Base(base) {}

    Convertible_iterator_project(const Self& it) : Base(it) {}
    Self& operator=(const Self& it) { Base::operator=(it); return *this; }

  /// @}

  /// \name Forward Category
  /// @{

    bool  operator==(Nullptr_t ptr) const { return (const Base&)*this == ptr; }
    bool  operator!=(Nullptr_t ptr) const { return ! (*this == ptr); }
    bool  operator==(const Self& it) const { return (const Base&)*this == it; }
    bool  operator!=(const Self& it) const { return ! (*this == it); }

    Self& operator++()     { Base::operator++(); return *this; }
    Self  operator++(int)  { Self tmp(*this); ++(*this); return tmp; }

  /// @}

  /// \name Bidirectional Category
  /// @{

    Self& operator--()     { Base::operator--(); return *this; }
    Self  operator--(int)  { Self tmp(*this); --(*this); return tmp; }

  /// @}

  /// \name Conversion
  /// @{
    operator Handle()               { return Base::operator->(); }
    operator ConstHandle() const    { return Base::operator->(); }
  /// @}

};

/// \endcond

} //namespace CGAL

#endif //CGAL_CONVERTIBLE_ITERATOR_PROJECT_H
