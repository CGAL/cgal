// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:  $
// $Id:  $
//
// Author(s)     : Clement Jamin

#ifndef CGAL_COMPACT_CONTAINER_STRATEGIES_H
#define CGAL_COMPACT_CONTAINER_STRATEGIES_H

#include <CGAL/tags.h>

namespace CGAL {

// A basic "do nothing" strategy
// One can inheritate from it for partial specialisation
class Compact_container_strategy_base {
public:
  typedef Tag_false Uses_erase_counter;

  // Do nothing
  template <typename Element>
  static unsigned int get_erase_counter(const Element &) { return 0; }
  template <typename Element>
  static void set_erase_counter(Element &, unsigned int) {}
  template <typename Element>
  static void increment_erase_counter(Element &) {}
};


// A strategy managing an internal counter
class Compact_container_strategy_with_counter
{
public:
  typedef Tag_true Uses_erase_counter;

  template <typename Element>
  static unsigned int get_erase_counter(const Element &e)
  {
    return e.get_erase_counter();
  }

  template <typename Element>
  static void set_erase_counter(Element &e, unsigned int c)
  {
    e.set_erase_counter(c);
  }

  template <typename Element>
  static void increment_erase_counter(Element &e)
  {
    e.increment_erase_counter();
  }
};

// CC_safe_handle is a helper that store a CC handle and its erase 
// counter value (value when the CC_safe_handle instance was created).
// The is_zombie() function allows to know if the pointee was erased since.
template <typename CC_iterator>
class CC_safe_handle
{
  typedef typename CC_iterator::Strategy Strategy;

public:
  CC_safe_handle(CC_iterator handle)
    : m_handle(handle)
    , m_erase_counter_value(Strategy::get_erase_counter(*handle))
  {
    CGAL_static_assertion(
      (boost::is_same<Strategy::Uses_erase_counter, Tag_true>::value));
  }

  bool is_zombie() const
  {
    return Strategy::get_erase_counter(*m_handle) != m_erase_counter_value;
  }

  CC_iterator get_cc_handle() const
  {
    return m_handle;
  }

protected:
  CC_iterator       m_handle;
  unsigned int    m_erase_counter_value;
};

template <typename CC_iterator>
CC_safe_handle<CC_iterator> make_cc_safe_handle(CC_iterator handle)
{
  return CC_safe_handle<CC_iterator>(handle);
}

} //namespace CGAL

#endif // CGAL_COMPACT_CONTAINER_STRATEGIES_H