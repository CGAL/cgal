// Copyright (c) 2007-2020  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jackson Campolattaro, Cédric Portaneri, Tong Zhao

#ifndef CGAL_ORTHTREE_TRAVERSAL_ITERATOR_H
#define CGAL_ORTHTREE_TRAVERSAL_ITERATOR_H

#include <CGAL/license/Orthtree.h>

#include <boost/function.hpp>
#include <boost/iterator/iterator_facade.hpp>

/// \cond SKIP_IN_MANUAL

namespace CGAL {

/*!
 * \ingroup PkgOrthtreeClasses
 *
 * \brief
 *
 * \todo
 *
 * \tparam Value
 */
template<class Value>
class Traversal_iterator :
        public boost::iterator_facade<Traversal_iterator<Value>, Value, boost::forward_traversal_tag> {

public:

  /// \name Types
  /// @{

  /*!
   * \brief
   *
   * \todo
   */
  typedef std::function<Value(Value)> Traversal_function;

  /// @}

public:

  /// \name Creation
  /// @{

  /*!
   * \brief
   *
   * \todo
   */
  Traversal_iterator() : m_value(), m_next() {}

  /*!
   * \brief
   *
   * \todo
   *
   * \param first
   * \param next
   */
  Traversal_iterator(Value first, const Traversal_function &next) : m_value(first), m_next(next) {}

  /// @}

private:
  friend class boost::iterator_core_access;

  bool equal(Traversal_iterator<Value> const &other) const {
    return m_value == other.m_value;
  }

  void increment() {
    m_value = m_next(m_value);
  }

  Value &dereference() const {
    return const_cast<Value&>(m_value);
  }

private:

  Value m_value;
  Traversal_function m_next;
};
}

/// \endcond

#endif //CGAL_ORTHTREE_TRAVERSAL_ITERATOR_H
