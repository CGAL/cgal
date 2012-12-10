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


#ifndef CGAL_CONVERTIBLE_FILTER_ITERATOR_H
#define CGAL_CONVERTIBLE_FILTER_ITERATOR_H

#include <CGAL/iterator.h>

namespace CGAL {

/// \cond SKIP_IN_MANUAL
  
/// This class inherits from Filter_iterator<> +
/// adds a conversion to handle/const handle.
/// See Filter_iterator<> documentation.

template<class Iterator,        ///< Internal iterator.
         class Predicate,       ///< Predicate to filter out elements.
         class ConstHandle,     ///< Const-handle type to convert to.
         class Handle = void*   ///< Non-const-handle type to convert to (void*=none).
>
class Convertible_filter_iterator
    : public Filter_iterator<Iterator, Predicate>
{
    typedef Filter_iterator<Iterator, Predicate>    Base;
    typedef Convertible_filter_iterator             Self;

public:

  /// \name Creation
  /// @{
    Convertible_filter_iterator() {}
    Convertible_filter_iterator(Base base)
        : Base(base) {}
    Convertible_filter_iterator(Iterator e, const Predicate& p)
        : Base(e,p) {}
    Convertible_filter_iterator(Iterator e, const Predicate& p, Iterator c)
        : Base(e,p,c) {}

    Convertible_filter_iterator(const Self& it) : Base(it) {}
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

#endif //CGAL_CONVERTIBLE_FILTER_ITERATOR_H
