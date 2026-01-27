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

#include <optional>

#include <boost/function.hpp>
#include <boost/optional.hpp>
#include <boost/iterator/iterator_facade.hpp>

/// \cond SKIP_IN_MANUAL

namespace CGAL {

/*!
 * \ingroup PkgOrthtreeTraversal
 *
 * \brief Wraps a traversal definition to produce an iterator which traverses the tree when incremented.
 *
 * \todo
 *
 * \tparam Tree The orthtree type to iterate over
 */
template <class Tree>
class Index_traversal_iterator : public boost::iterator_facade<
  Index_traversal_iterator<Tree>,
  const typename Tree::Node_index,
  boost::forward_traversal_tag,
  const typename Tree::Node_index
> {
public:

  /// \name Types
  /// @{

  /*!
   * \brief
   *
   * \todo
   */
  using Traversal_function = std::function<std::optional<std::size_t>(const Tree&, std::size_t)>;

  using Node_index = typename Tree::Node_index;

  /// @}

public:

  /// \name Creation
  /// @{

  /*!
   * \brief Default constructor, creates an end sentinel
   *
   * \todo
   */
  Index_traversal_iterator() : m_next() {}

  /*!
   * \brief
   *
   * \todo
   *
   * \param tree
   * \param first
   * \param next
   */
  Index_traversal_iterator(const Tree& tree, Node_index first, const Traversal_function& next) :
    m_tree(&tree), m_index(first), m_next(next)  {}

  /// @}

private:

  friend class boost::iterator_core_access;

  bool equal(Index_traversal_iterator<Tree> const& other) const {
    return m_index == other.m_index;
  }

  void increment() {
    // invoking increment on the sentinel is undefined behavior
    m_index = m_next(*m_tree, *m_index);
  }

  Node_index dereference() const {
    return *m_index;
  }

private:

  const Tree* m_tree = nullptr;
  std::optional<std::size_t> m_index;
  Traversal_function m_next;

};

}

/// \endcond

#endif //CGAL_ORTHTREE_TRAVERSAL_ITERATOR_H
