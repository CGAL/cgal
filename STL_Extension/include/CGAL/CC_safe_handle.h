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
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Clement Jamin

#ifndef CGAL_CC_SAFE_HANDLE_H
#define CGAL_CC_SAFE_HANDLE_H

namespace CGAL {

// CC_safe_handle is a helper that store a CC handle and its erase 
// counter value (value when the CC_safe_handle instance was created).
// The is_zombie() function allows to know if the pointee was erased since.
template <typename CC_iterator>
class CC_safe_handle
{
public:
  CC_safe_handle(CC_iterator iterator)
    : m_handle(iterator)
    , m_erase_counter_value(iterator->erase_counter())
  {
  }

  bool is_zombie() const
  {
    return m_handle->erase_counter() != m_erase_counter_value;
  }

  CC_iterator cc_iterator() const
  {
    return m_handle;
  }

protected:
  CC_iterator     m_handle;
  unsigned int    m_erase_counter_value;
};

template <typename CC_iterator>
CC_safe_handle<CC_iterator> make_cc_safe_handle(CC_iterator iterator)
{
  return CC_safe_handle<CC_iterator>(iterator);
}

} //namespace CGAL

#endif // CGAL_CC_SAFE_HANDLE_H
