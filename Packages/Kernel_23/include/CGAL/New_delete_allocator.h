// Copyright (c) 2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_NEW_DELETE_ALLOCATOR_H
#define CGAL_NEW_DELETE_ALLOCATOR_H

#include <cstddef> 

CGAL_BEGIN_NAMESPACE

template <class T>
class New_delete_allocator
{
 public:
  typedef std::size_t      size_type;
  typedef std::ptrdiff_t   difference_type;
  typedef T                value_type;
  typedef T*               pointer;
  typedef const T*         const_pointer;
  typedef T&               reference;
  typedef const T&         const_reference;

  template <class U> struct rebind { typedef New_delete_allocator<U> other; };

  New_delete_allocator() {}

  // Uses the default constructor of T.
  pointer
  allocate(size_type n, const_pointer = 0) const
  { return new T[n]; }

  // p should have an appropriate type
  // That's why we need the parameterization
  // T should have a virtual destructor
  void
  deallocate(pointer p, size_type) const
  { delete[] p; }

  pointer
  address(reference) const
  { return static_cast<pointer>(0); }

  const_pointer
  address(const_reference) const
  { return static_cast<const_pointer>(0); }

  void
  construct(pointer ptr, const_reference ref) const
  { *ptr = ref; }

  // We can't do anything here.
  void
  destroy(pointer) const
  { }

  size_type
  max_size() const
  { return 0; }
};

CGAL_END_NAMESPACE

#endif // CGAL_NEW_DELETE_ALLOCATOR_H
