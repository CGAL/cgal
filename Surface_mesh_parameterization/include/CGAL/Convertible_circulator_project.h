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


#ifndef CGAL_CONVERTIBLE_CIRCULATOR_PROJECT_H
#define CGAL_CONVERTIBLE_CIRCULATOR_PROJECT_H

#include <CGAL/Circulator_project.h>

namespace CGAL {

/// \cond SKIP_IN_MANUAL

/// This class inherits from Circulator_project<> +
/// adds a conversion to handle/const handle.
/// See Circulator_project<> documentation.
template<class C,               ///< Internal circulator.
         class Fct,             ///< Conversion functor.
         class Ref,
         class Ptr,
         class ConstHandle,     ///< Const-handle type to convert to.
         class Handle = void*   ///< Non-const-handle type to convert to (void*=none).
>
class Convertible_circulator_project
    : public Circulator_project<C, Fct, Ref, Ptr>
{
    typedef Circulator_project<C, Fct, Ref, Ptr>    Base;
    typedef Convertible_circulator_project          Self;

public:

  /// \name Creation
  /// @{

    Convertible_circulator_project() {}
    Convertible_circulator_project(Base base) : Base(base) {}

    Convertible_circulator_project(const Self& cir) : Base(cir) {}
    Self& operator=(const Self& cir) { Base::operator=(cir); return *this; }

  /// @}

  /// \name Forward Category
  ///@{

    bool  operator==(Nullptr_t ptr) const { return (const Base&)*this == ptr; }
    bool  operator!=(Nullptr_t ptr) const { return ! (*this == ptr); }
    bool  operator==(const Self& cir)    const { return (const Base&)*this == cir; }
    bool  operator!=(const Self& cir)    const { return ! (*this == cir); }

    Self& operator++()     { Base::operator++(); return *this; }
    Self  operator++(int)  { Self tmp(*this); ++(*this); return tmp; }
  ///@}

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

#endif //CGAL_CONVERTIBLE_CIRCULATOR_PROJECT_H
