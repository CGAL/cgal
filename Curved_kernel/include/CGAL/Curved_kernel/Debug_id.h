// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Monique Teillaud, Sylvain Pion

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CURVED_KERNEL_DEBUG_ID_H
#define CGAL_CURVED_KERNEL_DEBUG_ID_H

namespace CGAL {
namespace CGALi {

// Debug_id is a base class storing an identifier, to help recognizing
// the geometric objects in a debugging dump.
// It is empty, unless CGAL_CURVED_KERNEL_DEBUG is defined.

#ifndef CGAL_CURVED_KERNEL_DEBUG
  template <typename T = void>
  struct Debug_id {
    const Debug_id & id() const { return *this; }
  };

  template < typename T >
  std::ostream &
  operator<< (std::ostream & os, const Debug_id<T> &)
  {
    return os;
  }
#else
  // I make it a template in order to avoid linking problems
  // with the static data member.
  template <typename T = void>
  class Debug_id {
    static int cnt;
  public:
    int _id;
    Debug_id() : _id(cnt++) {}
    const Debug_id & id() const { return *this; }
  };

  template <typename T>
  int Debug_id<T>::cnt = 0;

  template < typename T >
  std::ostream &
  operator<< (std::ostream & os, const Debug_id<T> &a)
  {
    return os << " [[[ DEBUG_ID = " << a._id << " ]]] ";
  }
#endif

} // namespace CGALi
} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_DEBUG_ID_H
