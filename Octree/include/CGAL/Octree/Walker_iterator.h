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

#ifndef CGAL_OCTREE_WALKER_ITERATOR_H
#define CGAL_OCTREE_WALKER_ITERATOR_H

#include <CGAL/license/Octree.h>

#include <boost/function.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include <CGAL/Octree/Node.h>
#include "Node.h"

namespace CGAL {

/*!
 * \ingroup PkgOctreeClasses
 *
 * \brief
 *
 * \todo
 *
 * \tparam Value
 */
template<class Value>
class Walker_iterator :
        public boost::iterator_facade<Walker_iterator<Value>, Value, boost::forward_traversal_tag> {

public:

  /// \name Types
  /// @{

  /*!
   * \brief
   *
   * \todo
   */
  typedef std::function<Value *(Value *)> Walker_function;

  /// @}

public:

  /// \name Creation
  /// @{

  /*!
   * \brief
   *
   * \todo
   */
  Walker_iterator() : m_value(nullptr), m_next() {}

  /*!
   * \brief
   *
   * \todo
   *
   * \param first
   * \param next
   */
  Walker_iterator(Value *first, const Walker_function &next) : m_value(first), m_next(next) {}

  /// @}

private:
  friend class boost::iterator_core_access;

  bool equal(Walker_iterator<Value> const &other) const {
    return m_value == other.m_value;
  }

  void increment() {
    m_value = m_next(m_value);
  }

  Value &dereference() const {
    return *m_value;
  }

private:

  Value *m_value;
  Walker_function m_next;
};
}

#endif //CGAL_OCTREE_WALKER_ITERATOR_H
