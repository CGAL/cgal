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

namespace CGAL {

// A basic "do nothing" strategy
// One can inheritate from it for partial specialisation
class Compact_container_strategy_base {
public:
  static const bool Uses_erase_counter = false;

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
  static const bool Uses_erase_counter = true;

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

} //namespace CGAL

#endif // CGAL_COMPACT_CONTAINER_STRATEGIES_H