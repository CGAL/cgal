// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_SURFACE_SWEEP_2_RANDOM_ACCESS_OUTPUT_ITERATOR_H
#define CGAL_SURFACE_SWEEP_2_RANDOM_ACCESS_OUTPUT_ITERATOR_H

#include <CGAL/license/Surface_sweep_2.h>

#include <CGAL/utility.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class
 * A random-access iterator that can automatically resize its container.
 */
template <typename Container_>
class Random_access_output_iterator {
public:
  typedef Container_                                Container;
  typedef typename Container::value_type            value_type;
  typedef Random_access_output_iterator<Container>  Self;

private:
  // Data members:
  Container* m_container;               // The container.
  unsigned int m_index;                 // The current index.

public:
  Random_access_output_iterator() {}

  Random_access_output_iterator(Container& _container, unsigned int _index = 0) :
    m_container(&_container),
    m_index(_index)
  {}

  value_type& operator*()
  {
    if(m_index >= m_container->capacity()) {
      m_container->reserve(2 * m_index + 1);
      m_container->resize(m_index+1);
    }
    else if(m_index >= m_container->size())
      m_container->resize(m_index+1);
    return (*m_container)[m_index];
  }

  Self& operator++()
  {
    ++m_index;
    return (*this);
  }

  Self operator++ (int)
  {
    Self temp = *this;
    ++m_index;
    return (temp);
  }

  Self& operator--()
  {
    --m_index;
    return (*this);
  }

  Self operator--(int)
  {
    Self temp = *this;
    --m_index;
    return (temp);
  }

  bool operator==(const Self& other) const
  {
    CGAL_precondition(m_container == other.m_container);
    return (m_index == other.m_index);
  }

  bool operator!=(const Self& other) const
  {
    CGAL_precondition(m_container == other.m_container);
    return !(*this == other);
  }

  unsigned int operator-(const Self& other) const
  {
    CGAL_precondition(m_container == other.m_container);
    return (m_index - other.m_index);
  }
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
