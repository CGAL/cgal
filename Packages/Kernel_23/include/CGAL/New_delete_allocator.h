// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : New_delete_allocator.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 
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
